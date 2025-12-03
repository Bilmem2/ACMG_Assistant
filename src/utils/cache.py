"""
Result Cache Module
==================

Provides a strict, validated caching layer for predictor and population data.
All cached entries are validated before use; invalid entries are rejected.

Key Design Principles:
- Cache is an optimization layer, not a data source of truth
- All cached entries are validated before use
- Invalid/corrupted entries are rejected (treated as cache miss)
- CacheKey is constructed from category, variant_id, source, and version
- Thread-safe with file locking

Author: Can Sevilmiş  
License: MIT License
"""

import json
import os
import hashlib
import threading
from dataclasses import dataclass, asdict
from typing import Optional, Any
from datetime import datetime, timedelta
from pathlib import Path


@dataclass
class CacheKey:
    """
    Cache key for predictor and population data.
    
    Attributes:
        category: Type of data ('predictor' or 'population')
        source: API/database source (e.g., 'dbNSFP', 'gnomAD_GraphQL')
        variant_id: Normalized variant identifier (e.g., 'GRCh38:17-7674234-G-A')
        version: API/database version (optional)
    """
    category: str  # "predictor" or "population"
    source: str    # "CADD_API", "dbNSFP", "gnomAD", etc.
    variant_id: str  # Normalized representation, e.g. "GRCh38:17-7674234-G-A"
    version: Optional[str] = None  # API/DB version
    
    def to_hash(self) -> str:
        """Generate a unique hash for this cache key."""
        key_string = f"{self.category}:{self.source}:{self.variant_id}:{self.version or ''}"
        return hashlib.sha256(key_string.encode()).hexdigest()[:32]
    
    def to_dict(self) -> dict:
        """Convert to dictionary for serialization."""
        return asdict(self)
    
    @classmethod
    def from_dict(cls, data: dict) -> 'CacheKey':
        """Create CacheKey from dictionary."""
        return cls(
            category=data['category'],
            source=data['source'],
            variant_id=data['variant_id'],
            version=data.get('version')
        )


@dataclass
class CacheEntry:
    """
    Cache entry with metadata and validation.
    
    Attributes:
        key: CacheKey for this entry
        value: Cached data (dictionary)
        timestamp: When entry was created (ISO format)
        valid_until: When entry expires (ISO format)
        validated: Whether entry passed validation
    """
    key: CacheKey
    value: dict
    timestamp: str
    valid_until: str
    validated: bool = True
    
    def is_expired(self) -> bool:
        """Check if cache entry has expired."""
        try:
            valid_until = datetime.fromisoformat(self.valid_until)
            return datetime.now() > valid_until
        except (ValueError, TypeError):
            return True  # Treat parsing errors as expired
    
    def to_dict(self) -> dict:
        """Convert to dictionary for serialization."""
        return {
            'key': self.key.to_dict(),
            'value': self.value,
            'timestamp': self.timestamp,
            'valid_until': self.valid_until,
            'validated': self.validated
        }
    
    @classmethod
    def from_dict(cls, data: dict) -> 'CacheEntry':
        """Create CacheEntry from dictionary."""
        return cls(
            key=CacheKey.from_dict(data['key']),
            value=data['value'],
            timestamp=data['timestamp'],
            valid_until=data['valid_until'],
            validated=data.get('validated', True)
        )


class ResultCache:
    """
    Strict, validated cache for predictor and population data.
    
    This cache:
    - Validates all entries before use
    - Rejects corrupted/invalid entries
    - Provides thread-safe access
    - Supports TTL-based expiration
    
    Usage:
        cache = ResultCache(cache_dir='./cache')
        key = CacheKey(category='predictor', source='dbNSFP', 
                       variant_id='GRCh38:17-7674234-G-A')
        
        # Get cached value (returns None if not found or invalid)
        value = cache.get(key)
        
        # Set cache value
        cache.set(key, {'revel': 0.85, 'cadd_phred': 25.0})
        
        # Invalidate entry
        cache.invalidate(key)
    """
    
    # Default TTL: 7 days for predictor data, 30 days for population data
    DEFAULT_TTL = {
        'predictor': timedelta(days=7),
        'population': timedelta(days=30),
    }
    
    def __init__(
        self,
        cache_dir: Optional[str] = None,
        ttl: Optional[timedelta] = None,
        enabled: bool = True
    ):
        """
        Initialize the result cache.
        
        Args:
            cache_dir: Directory to store cache files (default: ./api_cache)
            ttl: Time-to-live for cache entries (default: category-specific)
            enabled: Whether caching is enabled
        """
        self.enabled = enabled
        self.ttl = ttl
        self._lock = threading.RLock()
        
        # Set up cache directory
        if cache_dir:
            self.cache_dir = Path(cache_dir)
        else:
            # Use src/api_cache by default
            self.cache_dir = Path(__file__).parent.parent / 'api_cache'
        
        # Create cache directory if it doesn't exist
        if self.enabled:
            self.cache_dir.mkdir(parents=True, exist_ok=True)
    
    def _get_cache_file(self, key: CacheKey) -> Path:
        """Get the cache file path for a given key."""
        # Organize by category and source
        category_dir = self.cache_dir / key.category / key.source
        category_dir.mkdir(parents=True, exist_ok=True)
        return category_dir / f"{key.to_hash()}.json"
    
    def _get_ttl(self, category: str) -> timedelta:
        """Get TTL for a category."""
        if self.ttl:
            return self.ttl
        return self.DEFAULT_TTL.get(category, timedelta(days=7))
    
    def get(self, key: CacheKey) -> Optional[dict]:
        """
        Get cached value if valid and not expired.
        
        Args:
            key: CacheKey to look up
            
        Returns:
            Cached value dict or None if not found/invalid/expired
        """
        if not self.enabled:
            return None
        
        with self._lock:
            cache_file = self._get_cache_file(key)
            
            if not cache_file.exists():
                return None
            
            try:
                with open(cache_file, 'r', encoding='utf-8') as f:
                    data = json.load(f)
                
                entry = CacheEntry.from_dict(data)
                
                # Check expiration
                if entry.is_expired():
                    self._remove_file(cache_file)
                    return None
                
                # Check if entry was validated
                if not entry.validated:
                    self._remove_file(cache_file)
                    return None
                
                # Verify key matches (sanity check)
                if entry.key.to_hash() != key.to_hash():
                    self._remove_file(cache_file)
                    return None
                
                return entry.value
                
            except (json.JSONDecodeError, KeyError, TypeError, ValueError) as e:
                # Corrupted cache file - remove and treat as miss
                self._remove_file(cache_file)
                return None
    
    def set(
        self,
        key: CacheKey,
        value: dict,
        validated: bool = True
    ) -> None:
        """
        Set a cache entry.
        
        Args:
            key: CacheKey for this entry
            value: Dictionary value to cache
            validated: Whether the value passed validation
        """
        if not self.enabled:
            return
        
        if not validated:
            # Do not cache invalid values
            return
        
        with self._lock:
            cache_file = self._get_cache_file(key)
            
            now = datetime.now()
            ttl = self._get_ttl(key.category)
            
            entry = CacheEntry(
                key=key,
                value=value,
                timestamp=now.isoformat(),
                valid_until=(now + ttl).isoformat(),
                validated=validated
            )
            
            try:
                with open(cache_file, 'w', encoding='utf-8') as f:
                    json.dump(entry.to_dict(), f, indent=2)
            except (OSError, IOError) as e:
                # Silently fail on write errors - cache is optional
                pass
    
    def invalidate(self, key: CacheKey) -> None:
        """
        Invalidate (remove) a cache entry.
        
        Args:
            key: CacheKey to invalidate
        """
        if not self.enabled:
            return
        
        with self._lock:
            cache_file = self._get_cache_file(key)
            self._remove_file(cache_file)
    
    def invalidate_all(self, category: Optional[str] = None) -> int:
        """
        Invalidate all cache entries, optionally filtered by category.
        
        Args:
            category: Optional category to filter ('predictor' or 'population')
            
        Returns:
            Number of entries invalidated
        """
        if not self.enabled:
            return 0
        
        count = 0
        with self._lock:
            if category:
                category_dir = self.cache_dir / category
                if category_dir.exists():
                    count = self._remove_directory_contents(category_dir)
            else:
                count = self._remove_directory_contents(self.cache_dir)
        
        return count
    
    def _remove_file(self, path: Path) -> None:
        """Safely remove a file."""
        try:
            if path.exists():
                path.unlink()
        except OSError:
            pass
    
    def _remove_directory_contents(self, path: Path) -> int:
        """Recursively remove directory contents and count files removed."""
        count = 0
        try:
            for item in path.rglob('*.json'):
                item.unlink()
                count += 1
        except OSError:
            pass
        return count
    
    def get_stats(self) -> dict:
        """Get cache statistics."""
        if not self.enabled:
            return {'enabled': False}
        
        stats = {
            'enabled': True,
            'cache_dir': str(self.cache_dir),
            'categories': {}
        }
        
        with self._lock:
            for category in ['predictor', 'population']:
                category_dir = self.cache_dir / category
                if category_dir.exists():
                    files = list(category_dir.rglob('*.json'))
                    stats['categories'][category] = {
                        'entries': len(files),
                        'sources': list(set(f.parent.name for f in files))
                    }
        
        return stats


# =============================================================================
# Variant ID Normalization Helpers
# =============================================================================

def normalize_variant_id(
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    genome_build: str = 'GRCh38'
) -> str:
    """
    Normalize variant to a canonical string identifier.
    
    Args:
        chrom: Chromosome (e.g., '17', 'chr17')
        pos: Genomic position
        ref: Reference allele
        alt: Alternate allele
        genome_build: Genome build (default: GRCh38)
        
    Returns:
        Normalized variant ID (e.g., 'GRCh38:17-7674234-G-A')
    """
    # Normalize chromosome (remove 'chr' prefix)
    norm_chrom = str(chrom).upper().replace('CHR', '')
    
    # Uppercase alleles
    norm_ref = str(ref).upper()
    norm_alt = str(alt).upper()
    
    return f"{genome_build}:{norm_chrom}-{pos}-{norm_ref}-{norm_alt}"


def normalize_variant_id_from_hgvs(
    gene: str,
    transcript: Optional[str],
    hgvs_c: str,
    genome_build: str = 'GRCh38'
) -> str:
    """
    Normalize variant to a canonical string from HGVS notation.
    
    Args:
        gene: Gene symbol
        transcript: Transcript ID (e.g., 'NM_000546.6')
        hgvs_c: cDNA HGVS notation (e.g., 'c.1528C>T')
        genome_build: Genome build (default: GRCh38)
        
    Returns:
        Normalized variant ID (e.g., 'GRCh38:TP53:NM_000546.6:c.1528C>T')
    """
    gene = gene.upper() if gene else 'UNKNOWN'
    transcript = transcript or 'UNKNOWN'
    
    # Normalize HGVS notation
    norm_hgvs = hgvs_c.strip()
    
    return f"{genome_build}:{gene}:{transcript}:{norm_hgvs}"


def build_predictor_cache_key(
    source: str,
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    version: Optional[str] = None,
    genome_build: str = 'GRCh38'
) -> CacheKey:
    """
    Build a CacheKey for predictor data.
    
    Args:
        source: API source name (e.g., 'dbNSFP', 'CADD_API')
        chrom: Chromosome
        pos: Genomic position
        ref: Reference allele
        alt: Alternate allele
        version: API/database version
        genome_build: Genome build
        
    Returns:
        CacheKey for predictor data
    """
    variant_id = normalize_variant_id(chrom, pos, ref, alt, genome_build)
    return CacheKey(
        category='predictor',
        source=source,
        variant_id=variant_id,
        version=version
    )


def build_population_cache_key(
    source: str,
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    version: Optional[str] = None,
    genome_build: str = 'GRCh38'
) -> CacheKey:
    """
    Build a CacheKey for population data.
    
    Args:
        source: API source name (e.g., 'gnomAD_GraphQL')
        chrom: Chromosome
        pos: Genomic position
        ref: Reference allele
        alt: Alternate allele
        version: Database version
        genome_build: Genome build
        
    Returns:
        CacheKey for population data
    """
    variant_id = normalize_variant_id(chrom, pos, ref, alt, genome_build)
    return CacheKey(
        category='population',
        source=source,
        variant_id=variant_id,
        version=version
    )

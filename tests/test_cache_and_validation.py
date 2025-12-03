"""
Test Suite for Cache and Validation Layer
==========================================

Tests for the strict, validated caching layer for predictor and population data.

Key test scenarios:
- Valid cache hit uses cached data
- Invalid cache entry (e.g., REVEL=5.0) rejected + invalidated + API fallback
- Population cache with invalid AF/AC/AN rejected
- Multi-source priority with one invalid, one valid source
- Cache expiration handling
- Corrupted cache file recovery

Author: Can Sevilmiş
License: MIT License
"""

import pytest
import json
import tempfile
import shutil
from pathlib import Path
from datetime import datetime, timedelta
from unittest.mock import patch, MagicMock

# Import modules under test
import sys
sys.path.insert(0, str(Path(__file__).parent.parent / 'src'))

from utils.cache import (
    CacheKey, CacheEntry, ResultCache,
    normalize_variant_id,
    build_predictor_cache_key, build_population_cache_key
)
from config.predictors import (
    PredictorScore, PopulationStats,
    validate_predictor_score, validate_population_stats,
    validate_predictor_score_object, validate_population_stats_object,
    validate_cached_predictor_data, validate_cached_population_data
)


# =============================================================================
# CacheKey Tests
# =============================================================================

class TestCacheKey:
    """Tests for CacheKey dataclass."""
    
    def test_cache_key_creation(self):
        """Test basic CacheKey creation."""
        key = CacheKey(
            category='predictor',
            source='dbNSFP',
            variant_id='GRCh38:17-7674234-G-A'
        )
        assert key.category == 'predictor'
        assert key.source == 'dbNSFP'
        assert key.variant_id == 'GRCh38:17-7674234-G-A'
        assert key.version is None
    
    def test_cache_key_with_version(self):
        """Test CacheKey with version."""
        key = CacheKey(
            category='population',
            source='gnomAD_GraphQL',
            variant_id='GRCh38:17-7674234-G-A',
            version='v4.0'
        )
        assert key.version == 'v4.0'
    
    def test_cache_key_hash_uniqueness(self):
        """Test that different keys have different hashes."""
        key1 = CacheKey('predictor', 'dbNSFP', 'GRCh38:17-7674234-G-A')
        key2 = CacheKey('predictor', 'dbNSFP', 'GRCh38:17-7674234-G-T')
        key3 = CacheKey('predictor', 'CADD_API', 'GRCh38:17-7674234-G-A')
        
        assert key1.to_hash() != key2.to_hash()  # Different variant
        assert key1.to_hash() != key3.to_hash()  # Different source
    
    def test_cache_key_hash_consistency(self):
        """Test that same key always produces same hash."""
        key1 = CacheKey('predictor', 'dbNSFP', 'GRCh38:17-7674234-G-A')
        key2 = CacheKey('predictor', 'dbNSFP', 'GRCh38:17-7674234-G-A')
        
        assert key1.to_hash() == key2.to_hash()
    
    def test_cache_key_to_dict(self):
        """Test CacheKey serialization."""
        key = CacheKey('predictor', 'dbNSFP', 'GRCh38:17-7674234-G-A', 'v4.0')
        data = key.to_dict()
        
        assert data == {
            'category': 'predictor',
            'source': 'dbNSFP',
            'variant_id': 'GRCh38:17-7674234-G-A',
            'version': 'v4.0'
        }
    
    def test_cache_key_from_dict(self):
        """Test CacheKey deserialization."""
        data = {
            'category': 'population',
            'source': 'gnomAD_GraphQL',
            'variant_id': 'GRCh38:17-7674234-G-A',
            'version': 'v4.0'
        }
        key = CacheKey.from_dict(data)
        
        assert key.category == 'population'
        assert key.source == 'gnomAD_GraphQL'
        assert key.version == 'v4.0'


# =============================================================================
# CacheEntry Tests
# =============================================================================

class TestCacheEntry:
    """Tests for CacheEntry dataclass."""
    
    def test_cache_entry_creation(self):
        """Test basic CacheEntry creation."""
        key = CacheKey('predictor', 'dbNSFP', 'GRCh38:17-7674234-G-A')
        now = datetime.now()
        entry = CacheEntry(
            key=key,
            value={'revel': 0.85},
            timestamp=now.isoformat(),
            valid_until=(now + timedelta(days=7)).isoformat()
        )
        
        assert entry.key == key
        assert entry.value == {'revel': 0.85}
        assert entry.validated is True
    
    def test_cache_entry_not_expired(self):
        """Test that fresh entry is not expired."""
        key = CacheKey('predictor', 'dbNSFP', 'test-variant')
        now = datetime.now()
        entry = CacheEntry(
            key=key,
            value={'revel': 0.85},
            timestamp=now.isoformat(),
            valid_until=(now + timedelta(days=7)).isoformat()
        )
        
        assert not entry.is_expired()
    
    def test_cache_entry_expired(self):
        """Test that old entry is expired."""
        key = CacheKey('predictor', 'dbNSFP', 'test-variant')
        past = datetime.now() - timedelta(days=10)
        entry = CacheEntry(
            key=key,
            value={'revel': 0.85},
            timestamp=past.isoformat(),
            valid_until=(past + timedelta(days=7)).isoformat()  # Expired 3 days ago
        )
        
        assert entry.is_expired()
    
    def test_cache_entry_serialization_roundtrip(self):
        """Test CacheEntry serialization and deserialization."""
        key = CacheKey('predictor', 'dbNSFP', 'GRCh38:17-7674234-G-A')
        now = datetime.now()
        entry = CacheEntry(
            key=key,
            value={'revel': 0.85, 'cadd_phred': 25.0},
            timestamp=now.isoformat(),
            valid_until=(now + timedelta(days=7)).isoformat()
        )
        
        data = entry.to_dict()
        restored = CacheEntry.from_dict(data)
        
        assert restored.key.to_hash() == entry.key.to_hash()
        assert restored.value == entry.value
        assert restored.validated == entry.validated


# =============================================================================
# ResultCache Tests
# =============================================================================

class TestResultCache:
    """Tests for ResultCache class."""
    
    @pytest.fixture
    def temp_cache_dir(self):
        """Create a temporary directory for cache tests."""
        temp_dir = tempfile.mkdtemp(prefix='acmg_cache_test_')
        yield temp_dir
        shutil.rmtree(temp_dir, ignore_errors=True)
    
    def test_cache_initialization(self, temp_cache_dir):
        """Test cache initialization."""
        cache = ResultCache(cache_dir=temp_cache_dir)
        
        assert cache.enabled is True
        assert cache.cache_dir == Path(temp_cache_dir)
    
    def test_cache_disabled(self):
        """Test that disabled cache returns None."""
        cache = ResultCache(enabled=False)
        key = CacheKey('predictor', 'dbNSFP', 'test-variant')
        
        cache.set(key, {'revel': 0.85})
        assert cache.get(key) is None
    
    def test_cache_set_and_get(self, temp_cache_dir):
        """Test basic set and get operations."""
        cache = ResultCache(cache_dir=temp_cache_dir)
        key = CacheKey('predictor', 'dbNSFP', 'GRCh38:17-7674234-G-A')
        value = {'revel': 0.85, 'cadd_phred': 25.0}
        
        cache.set(key, value)
        result = cache.get(key)
        
        assert result == value
    
    def test_cache_miss_returns_none(self, temp_cache_dir):
        """Test that cache miss returns None."""
        cache = ResultCache(cache_dir=temp_cache_dir)
        key = CacheKey('predictor', 'dbNSFP', 'nonexistent-variant')
        
        assert cache.get(key) is None
    
    def test_cache_invalidate(self, temp_cache_dir):
        """Test cache invalidation."""
        cache = ResultCache(cache_dir=temp_cache_dir)
        key = CacheKey('predictor', 'dbNSFP', 'test-variant')
        
        cache.set(key, {'revel': 0.85})
        assert cache.get(key) is not None
        
        cache.invalidate(key)
        assert cache.get(key) is None
    
    def test_cache_invalidate_all(self, temp_cache_dir):
        """Test invalidating all cache entries."""
        cache = ResultCache(cache_dir=temp_cache_dir)
        
        # Add multiple entries
        for i in range(5):
            key = CacheKey('predictor', 'dbNSFP', f'variant-{i}')
            cache.set(key, {'revel': 0.5 + i * 0.1})
        
        count = cache.invalidate_all()
        assert count == 5
        
        # Verify all are gone
        for i in range(5):
            key = CacheKey('predictor', 'dbNSFP', f'variant-{i}')
            assert cache.get(key) is None
    
    def test_cache_rejects_invalid_data(self, temp_cache_dir):
        """Test that cache rejects invalid (unvalidated) data."""
        cache = ResultCache(cache_dir=temp_cache_dir)
        key = CacheKey('predictor', 'dbNSFP', 'test-variant')
        
        # Try to set with validated=False
        cache.set(key, {'revel': 999.0}, validated=False)
        
        # Should not be cached
        assert cache.get(key) is None
    
    def test_cache_handles_corrupted_file(self, temp_cache_dir):
        """Test that cache handles corrupted JSON gracefully."""
        cache = ResultCache(cache_dir=temp_cache_dir)
        key = CacheKey('predictor', 'dbNSFP', 'test-variant')
        
        # Set valid data first
        cache.set(key, {'revel': 0.85})
        
        # Corrupt the file
        cache_file = cache._get_cache_file(key)
        with open(cache_file, 'w') as f:
            f.write('not valid json {{{')
        
        # Should return None (and remove corrupted file)
        assert cache.get(key) is None
    
    def test_cache_expired_entry_returns_none(self, temp_cache_dir):
        """Test that expired entries are not returned."""
        cache = ResultCache(
            cache_dir=temp_cache_dir,
            ttl=timedelta(seconds=-1)  # Already expired
        )
        key = CacheKey('predictor', 'dbNSFP', 'test-variant')
        
        cache.set(key, {'revel': 0.85})
        
        # Entry should be expired immediately
        assert cache.get(key) is None
    
    def test_cache_stats(self, temp_cache_dir):
        """Test cache statistics."""
        cache = ResultCache(cache_dir=temp_cache_dir)
        
        # Add entries
        key1 = CacheKey('predictor', 'dbNSFP', 'variant-1')
        key2 = CacheKey('predictor', 'CADD_API', 'variant-2')
        key3 = CacheKey('population', 'gnomAD_GraphQL', 'variant-1')
        
        cache.set(key1, {'revel': 0.85})
        cache.set(key2, {'cadd_phred': 25.0})
        cache.set(key3, {'af': 0.001})
        
        stats = cache.get_stats()
        
        assert stats['enabled'] is True
        assert 'predictor' in stats['categories']
        assert 'population' in stats['categories']


# =============================================================================
# Predictor Score Validation Tests
# =============================================================================

class TestPredictorScoreValidation:
    """Tests for predictor score validation."""
    
    def test_validate_revel_valid(self):
        """Test valid REVEL scores."""
        assert validate_predictor_score('revel', 0.0) is True
        assert validate_predictor_score('revel', 0.5) is True
        assert validate_predictor_score('revel', 1.0) is True
        assert validate_predictor_score('revel', 0.85) is True
    
    def test_validate_revel_invalid_out_of_range(self):
        """Test invalid REVEL scores (out of range)."""
        assert validate_predictor_score('revel', -0.1) is False
        assert validate_predictor_score('revel', 1.1) is False
        assert validate_predictor_score('revel', 5.0) is False
    
    def test_validate_cadd_phred_valid(self):
        """Test valid CADD PHRED scores."""
        assert validate_predictor_score('cadd_phred', 0.0) is True
        assert validate_predictor_score('cadd_phred', 25.0) is True
        assert validate_predictor_score('cadd_phred', 60.0) is True
    
    def test_validate_cadd_phred_invalid(self):
        """Test invalid CADD PHRED scores."""
        assert validate_predictor_score('cadd_phred', -5.0) is False
        assert validate_predictor_score('cadd_phred', 100.0) is False
    
    def test_validate_none_value(self):
        """Test that None is valid (means no data)."""
        assert validate_predictor_score('revel', None) is True
        assert validate_predictor_score('cadd_phred', None) is True
    
    def test_validate_none_strict_mode(self):
        """Test that None is invalid in strict mode."""
        assert validate_predictor_score('revel', None, strict=True) is False
    
    def test_validate_nan_inf(self):
        """Test that NaN and Inf are invalid."""
        import math
        assert validate_predictor_score('revel', math.nan) is False
        assert validate_predictor_score('revel', math.inf) is False
        assert validate_predictor_score('revel', -math.inf) is False
    
    def test_validate_wrong_type(self):
        """Test that non-numeric types are invalid."""
        assert validate_predictor_score('revel', 'high') is False
        assert validate_predictor_score('revel', [0.5]) is False
    
    def test_validate_predictor_score_object(self):
        """Test validation of PredictorScore objects."""
        valid_score = PredictorScore(
            predictor='revel',
            value=0.85,
            source='dbNSFP'
        )
        assert validate_predictor_score_object(valid_score) is True
        
        invalid_score = PredictorScore(
            predictor='revel',
            value=5.0,  # Out of range
            source='dbNSFP'
        )
        assert validate_predictor_score_object(invalid_score) is False


# =============================================================================
# Population Stats Validation Tests
# =============================================================================

class TestPopulationStatsValidation:
    """Tests for population statistics validation."""
    
    def test_validate_af_valid(self):
        """Test valid allele frequencies."""
        assert validate_population_stats(af=0.0) is True
        assert validate_population_stats(af=0.001) is True
        assert validate_population_stats(af=0.5) is True
        assert validate_population_stats(af=1.0) is True
    
    def test_validate_af_invalid(self):
        """Test invalid allele frequencies."""
        assert validate_population_stats(af=-0.001) is False
        assert validate_population_stats(af=1.1) is False
        assert validate_population_stats(af=5.0) is False
    
    def test_validate_an_ac_valid(self):
        """Test valid allele number and count."""
        assert validate_population_stats(an=10000, ac=10) is True
        assert validate_population_stats(an=100, ac=0) is True
        assert validate_population_stats(an=1000, ac=1000) is True  # 100% AF
    
    def test_validate_ac_exceeds_an_invalid(self):
        """Test that AC > AN is invalid."""
        assert validate_population_stats(an=100, ac=200) is False
    
    def test_validate_an_zero_with_ac_invalid(self):
        """Test that AN=0 with AC present is invalid."""
        assert validate_population_stats(an=0, ac=5) is False
    
    def test_validate_all_none_valid(self):
        """Test that all None values are valid (no data)."""
        assert validate_population_stats() is True
    
    def test_validate_strict_mode(self):
        """Test strict mode requires at least AF or (AC, AN)."""
        assert validate_population_stats(strict=True) is False
        assert validate_population_stats(af=0.001, strict=True) is True
        assert validate_population_stats(an=1000, ac=10, strict=True) is True
    
    def test_validate_population_stats_object(self):
        """Test validation of PopulationStats objects."""
        valid_stats = PopulationStats(
            population='gnomad_v4',
            af=0.001,
            an=100000,
            ac=100
        )
        assert validate_population_stats_object(valid_stats) is True
        
        invalid_stats = PopulationStats(
            population='gnomad_v4',
            af=5.0  # Invalid
        )
        assert validate_population_stats_object(invalid_stats) is False


# =============================================================================
# Cached Data Validation Tests
# =============================================================================

class TestCachedDataValidation:
    """Tests for cached data dictionary validation."""
    
    def test_validate_cached_predictor_data_valid(self):
        """Test valid cached predictor data."""
        data = {
            'revel': 0.85,
            'cadd_phred': 25.0,
            'sift': 0.01
        }
        is_valid, errors = validate_cached_predictor_data(data)
        assert is_valid is True
        assert len(errors) == 0
    
    def test_validate_cached_predictor_data_invalid_score(self):
        """Test cached predictor data with invalid score."""
        data = {
            'revel': 5.0,  # Out of range
            'cadd_phred': 25.0
        }
        is_valid, errors = validate_cached_predictor_data(data)
        assert is_valid is False
        assert 'revel' in errors[0]
    
    def test_validate_cached_predictor_data_partial_valid(self):
        """Test that valid subset is accepted."""
        data = {
            'revel': 0.85,
            'unknown_predictor': 999.0  # Unknown predictor, any value OK
        }
        is_valid, errors = validate_cached_predictor_data(data)
        assert is_valid is True
    
    def test_validate_cached_population_data_valid(self):
        """Test valid cached population data."""
        data = {
            'af': 0.001,
            'an': 100000,
            'ac': 100
        }
        is_valid, errors = validate_cached_population_data(data)
        assert is_valid is True
    
    def test_validate_cached_population_data_invalid(self):
        """Test invalid cached population data."""
        data = {
            'af': 2.0,  # AF > 1
            'an': 100000
        }
        is_valid, errors = validate_cached_population_data(data)
        assert is_valid is False


# =============================================================================
# Variant ID Normalization Tests
# =============================================================================

class TestVariantIdNormalization:
    """Tests for variant ID normalization."""
    
    def test_normalize_variant_id_basic(self):
        """Test basic variant ID normalization."""
        result = normalize_variant_id('17', 7674234, 'G', 'A')
        assert result == 'GRCh38:17-7674234-G-A'
    
    def test_normalize_variant_id_chr_prefix(self):
        """Test chromosome prefix stripping."""
        result = normalize_variant_id('chr17', 7674234, 'G', 'A')
        assert result == 'GRCh38:17-7674234-G-A'
    
    def test_normalize_variant_id_lowercase(self):
        """Test allele case normalization."""
        result = normalize_variant_id('17', 7674234, 'g', 'a')
        assert result == 'GRCh38:17-7674234-G-A'
    
    def test_build_predictor_cache_key(self):
        """Test predictor cache key building."""
        key = build_predictor_cache_key(
            source='dbNSFP',
            chrom='17',
            pos=7674234,
            ref='G',
            alt='A'
        )
        assert key.category == 'predictor'
        assert key.source == 'dbNSFP'
        assert 'GRCh38:17-7674234-G-A' in key.variant_id
    
    def test_build_population_cache_key(self):
        """Test population cache key building."""
        key = build_population_cache_key(
            source='gnomAD_GraphQL',
            chrom='17',
            pos=7674234,
            ref='G',
            alt='A',
            version='v4.0'
        )
        assert key.category == 'population'
        assert key.source == 'gnomAD_GraphQL'
        assert key.version == 'v4.0'


# =============================================================================
# Integration Tests with API Clients
# =============================================================================

class TestCacheIntegrationWithAPIClients:
    """Integration tests for cache with API clients."""
    
    @pytest.fixture
    def temp_cache_dir(self):
        """Create a temporary directory for cache tests."""
        temp_dir = tempfile.mkdtemp(prefix='acmg_cache_test_')
        yield temp_dir
        shutil.rmtree(temp_dir, ignore_errors=True)
    
    def test_predictor_client_uses_validated_cache(self, temp_cache_dir):
        """Test that PredictorAPIClient uses validated cache."""
        from utils.predictor_api_client import PredictorAPIClient
        
        cache = ResultCache(cache_dir=temp_cache_dir)
        client = PredictorAPIClient(
            api_enabled=False,  # Disable API for test
            result_cache=cache
        )
        
        # Pre-populate cache with valid data
        key = build_predictor_cache_key('dbNSFP', '17', 7674234, 'G', 'A')
        cache.set(key, {'revel': 0.85, 'cadd_phred': 25.0})
        
        # Client should use validated cache
        assert client._use_validated_cache is True
        cached = client._get_cached_predictor_data('dbNSFP', '17', 7674234, 'G', 'A')
        assert cached is not None
        assert cached['revel'] == 0.85
    
    def test_predictor_client_rejects_invalid_cache(self, temp_cache_dir):
        """Test that PredictorAPIClient rejects invalid cached data."""
        from utils.predictor_api_client import PredictorAPIClient
        
        cache = ResultCache(cache_dir=temp_cache_dir)
        
        # Manually write invalid data to cache file
        key = build_predictor_cache_key('dbNSFP', '17', 7674234, 'G', 'A')
        cache.set(key, {'revel': 5.0})  # Invalid REVEL score
        
        # Client should reject this
        client = PredictorAPIClient(
            api_enabled=False,
            result_cache=cache
        )
        
        cached = client._get_cached_predictor_data('dbNSFP', '17', 7674234, 'G', 'A')
        assert cached is None  # Invalid data rejected
    
    def test_population_client_uses_validated_cache(self, temp_cache_dir):
        """Test that PopulationAPIClient uses validated cache."""
        from utils.predictor_api_client import PopulationAPIClient
        
        cache = ResultCache(cache_dir=temp_cache_dir)
        client = PopulationAPIClient(
            api_enabled=False,
            result_cache=cache
        )
        
        # Pre-populate cache with valid data
        key = build_population_cache_key('gnomAD_GraphQL', '17', 7674234, 'G', 'A')
        cache.set(key, {'af': 0.001, 'an': 100000, 'ac': 100})
        
        # Client should use validated cache
        assert client._use_validated_cache is True
        cached = client._get_cached_population_data('gnomAD_GraphQL', '17', 7674234, 'G', 'A')
        assert cached is not None
        assert cached['af'] == 0.001
    
    def test_population_client_rejects_invalid_af(self, temp_cache_dir):
        """Test that PopulationAPIClient rejects invalid AF in cache."""
        from utils.predictor_api_client import PopulationAPIClient
        
        cache = ResultCache(cache_dir=temp_cache_dir)
        
        # Manually write invalid data to cache
        key = build_population_cache_key('gnomAD_GraphQL', '17', 7674234, 'G', 'A')
        cache.set(key, {'af': 2.0})  # Invalid AF > 1
        
        client = PopulationAPIClient(
            api_enabled=False,
            result_cache=cache
        )
        
        cached = client._get_cached_population_data('gnomAD_GraphQL', '17', 7674234, 'G', 'A')
        assert cached is None  # Invalid data rejected


if __name__ == '__main__':
    pytest.main([__file__, '-v'])

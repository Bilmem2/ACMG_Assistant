"""
HGVS Parser Module
==================

This module provides utilities for parsing HGVS-formatted variant strings.
Supports multiple HGVS formats including full RefSeq notation.

Supported formats:
- Full HGVS: NM_000546.6:c.1528C>T
- cDNA only: c.1528C>T
- Position only: 1528C>T
"""

import re
from typing import Dict, Optional, Tuple


class HGVSParser:
    """Parser for HGVS variant notation."""
    
    # Comprehensive HGVS patterns
    PATTERNS = {
        # Full format with RefSeq: NM_000546.6:c.1528C>T
        'full': r'^(?P<refseq>NM_\d+\.\d+):(?P<notation>[cgmrn])\.(?P<position>[\d\+\-\*]+)(?P<ref>[A-Z]+)>(?P<alt>[A-Z]+)$',
        
        # cDNA format: c.1528C>T
        'cdna': r'^(?P<notation>c)\.(?P<position>[\d\+\-\*]+)(?P<ref>[A-Z]+)>(?P<alt>[A-Z]+)$',
        
        # Simple position format: 1528C>T
        'simple': r'^(?P<position>[\d\+\-\*]+)(?P<ref>[A-Z]+)>(?P<alt>[A-Z]+)$',
        
        # Deletion format: c.1528del or NM_000546.6:c.1528del
        'deletion': r'^((?P<refseq>NM_\d+\.\d+):)?(?P<notation>[cgmrn])\.(?P<position>[\d\+\-\*]+(_[\d\+\-\*]+)?)del(?P<ref>[A-Z]+)?$',
        
        # Insertion format: c.1528_1529insA or NM_000546.6:c.1528_1529insA
        'insertion': r'^((?P<refseq>NM_\d+\.\d+):)?(?P<notation>[cgmrn])\.(?P<position>[\d\+\-\*]+_[\d\+\-\*]+)ins(?P<alt>[A-Z]+)$',
        
        # Duplication format: c.1528dup or NM_000546.6:c.1528dupC
        'duplication': r'^((?P<refseq>NM_\d+\.\d+):)?(?P<notation>[cgmrn])\.(?P<position>[\d\+\-\*]+(_[\d\+\-\*]+)?)dup(?P<ref>[A-Z]+)?$'
    }
    
    @classmethod
    def parse(cls, variant_string: str) -> Optional[Dict[str, str]]:
        """
        Parse a variant string in HGVS format.
        
        Args:
            variant_string: The variant string to parse
            
        Returns:
            Dict with parsed components or None if parsing fails:
            {
                'refseq_id': 'NM_000546.6' (optional),
                'notation_type': 'c',
                'position': '1528',
                'ref_base': 'C',
                'alt_base': 'T',
                'variant_type': 'substitution' | 'deletion' | 'insertion' | 'duplication',
                'original': original input string
            }
        """
        if not variant_string:
            return None
        
        variant_string = variant_string.strip()
        
        # Try each pattern in order
        for variant_type, pattern in cls.PATTERNS.items():
            match = re.match(pattern, variant_string, re.IGNORECASE)
            if match:
                result = match.groupdict()
                
                # Add metadata
                result['original'] = variant_string
                
                # Determine variant type
                if variant_type == 'full' or variant_type == 'cdna' or variant_type == 'simple':
                    result['variant_type'] = 'substitution'
                else:
                    result['variant_type'] = variant_type
                
                # Set default notation type if not present
                if 'notation' not in result or result['notation'] is None:
                    result['notation'] = 'c'  # Default to coding DNA
                
                # Clean up None values
                result = {k: v for k, v in result.items() if v is not None}
                
                # Standardize key names
                if 'refseq' in result:
                    result['refseq_id'] = result.pop('refseq')
                if 'notation' in result:
                    result['notation_type'] = result.pop('notation')
                if 'ref' in result:
                    result['ref_base'] = result.pop('ref')
                if 'alt' in result:
                    result['alt_base'] = result.pop('alt')
                
                return result
        
        return None
    
    @classmethod
    def validate(cls, variant_string: str) -> bool:
        """
        Validate if a string matches any HGVS pattern.
        
        Args:
            variant_string: The variant string to validate
            
        Returns:
            True if valid HGVS format, False otherwise
        """
        return cls.parse(variant_string) is not None
    
    @classmethod
    def extract_position_and_bases(cls, variant_string: str) -> Optional[Tuple[str, str, str]]:
        """
        Extract position and bases from variant string.
        
        Args:
            variant_string: The variant string to parse
            
        Returns:
            Tuple of (position, ref_base, alt_base) or None
        """
        parsed = cls.parse(variant_string)
        if not parsed:
            return None
        
        position = parsed.get('position', '')
        ref_base = parsed.get('ref_base', '')
        alt_base = parsed.get('alt_base', '')
        
        return (position, ref_base, alt_base)
    
    @classmethod
    def to_simple_format(cls, variant_string: str) -> Optional[str]:
        """
        Convert any HGVS format to simple c.notation format.
        
        Args:
            variant_string: The variant string to convert
            
        Returns:
            Simple format string (e.g., 'c.1528C>T') or None if parsing fails
        """
        parsed = cls.parse(variant_string)
        if not parsed:
            return None
        
        notation = parsed.get('notation_type', 'c')
        position = parsed.get('position', '')
        
        if parsed.get('variant_type') == 'substitution':
            ref_base = parsed.get('ref_base', '')
            alt_base = parsed.get('alt_base', '')
            return f"{notation}.{position}{ref_base}>{alt_base}"
        elif parsed.get('variant_type') == 'deletion':
            ref_base = parsed.get('ref_base', '')
            return f"{notation}.{position}del{ref_base}" if ref_base else f"{notation}.{position}del"
        elif parsed.get('variant_type') == 'insertion':
            alt_base = parsed.get('alt_base', '')
            return f"{notation}.{position}ins{alt_base}"
        elif parsed.get('variant_type') == 'duplication':
            ref_base = parsed.get('ref_base', '')
            return f"{notation}.{position}dup{ref_base}" if ref_base else f"{notation}.{position}dup"
        
        return None
    
    @classmethod
    def to_full_format(cls, variant_string: str, refseq_id: Optional[str] = None) -> Optional[str]:
        """
        Convert any HGVS format to full format with RefSeq ID.
        
        Args:
            variant_string: The variant string to convert
            refseq_id: Optional RefSeq ID to use (e.g., 'NM_000546.6')
            
        Returns:
            Full format string (e.g., 'NM_000546.6:c.1528C>T') or None if parsing fails
        """
        parsed = cls.parse(variant_string)
        if not parsed:
            return None
        
        # Use existing RefSeq ID if present, otherwise use provided one
        refseq = parsed.get('refseq_id', refseq_id)
        if not refseq:
            return None
        
        simple = cls.to_simple_format(variant_string)
        if not simple:
            return None
        
        return f"{refseq}:{simple}"


def parse_hgvs_variant(variant_string: str) -> Optional[Dict[str, str]]:
    """
    Convenience function to parse HGVS variant string.
    
    Args:
        variant_string: The variant string to parse
        
    Returns:
        Dict with parsed components or None if parsing fails
    """
    return HGVSParser.parse(variant_string)


def validate_hgvs_variant(variant_string: str) -> bool:
    """
    Convenience function to validate HGVS variant string.
    
    Args:
        variant_string: The variant string to validate
        
    Returns:
        True if valid HGVS format, False otherwise
    """
    return HGVSParser.validate(variant_string)

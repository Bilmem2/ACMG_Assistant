"""
Version and metadata information for ACMG Variant Classification Assistant
=========================================================================

Author: Can Sevilmiş
License: MIT License
"""

# Version and metadata information
VERSION_INFO = {
    'version': '4.0.0',
    'author': 'Can Sevilmiş',
    'license': 'MIT License',
    'last_updated': 'December 3, 2025',
    'guidelines': 'ACMG/AMP 2015 & 2023',
    'description': 'ACMG Variant Classification Assistant - Multi-Source API Integration with Validated Caching',
    'major_features': [
        'Dual-mode ACMG 2015/2023 guidelines support',
        'Complete 28 ACMG/AMP criteria implementation',
        'Multi-source predictor system (myvariant.info, AlphaMissense, CADD)',
        'Multi-source population data (gnomAD GraphQL, ExAC, TOPMed)',
        'Strict validated caching with TTL expiration',
        'PM1 via CancerHotspots API + UniProt domains',
        'HPO-based phenotype matching (PP4/BP5)',
        'Interactive evidence collection (PS3/BS3, PS4, PP1/BS4, PS1/PM5, PP5/BP6)',
        'Enhanced computational metascore for PP3/BP4',
        'Dynamic evidence weighting with confidence tracking',
        '🔬 gnomAD v4 constraint API (LOF intolerance)',
        '🔬 ClinGen eRepo API (disease-specific LOF mechanism)',
        '🔬 Multi-source LOF validation (population + clinical)',
        '🔬 Transparent data provenance tracking',
        '✨ Gene-specific BA1/BS1 thresholds',
        '✨ Automated Fisher\'s exact test (PS4)',
        '✨ Automated LOD scoring (PP1/BS4)',
        '✨ Strict PS2 2023 upgrade rules',
        '✨ Enhanced PP5/BP6 source validation',
        '✨ HGVS format support with auto-extraction',
        '🆕 ResultCache with strict validation',
        '🆕 198 comprehensive tests passing'
    ]
}

# Convenience access
__version__ = VERSION_INFO['version']
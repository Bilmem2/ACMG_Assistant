"""
User Interface Configuration for ACMG Variant Classification Assistant
=====================================================================

Contains colors, messages, and formatting configurations for the CLI interface.

Author: Can Sevilmiş
License: MIT License
"""

# Colorama color codes for terminal output
COLORAMA_COLORS = {
    'RED': '\033[91m',
    'GREEN': '\033[92m',
    'YELLOW': '\033[93m',
    'BLUE': '\033[94m',
    'MAGENTA': '\033[95m',
    'CYAN': '\033[96m',
    'WHITE': '\033[97m',
    'ENDC': '\033[0m',  # End color
    'RESET': '\033[0m',  # Reset color (alias for ENDC)
    'BOLD': '\033[1m',
    'UNDERLINE': '\033[4m'
}

# Section headers for organized output
SECTION_HEADERS = {
    'variant_info': f"{COLORAMA_COLORS['CYAN']}🧬 VARIANT INFORMATION{COLORAMA_COLORS['ENDC']}",
    'pathogenic_evidence': f"{COLORAMA_COLORS['RED']}🔴 PATHOGENIC EVIDENCE{COLORAMA_COLORS['ENDC']}",
    'benign_evidence': f"{COLORAMA_COLORS['GREEN']}🟢 BENIGN EVIDENCE{COLORAMA_COLORS['ENDC']}",
    'final_classification': f"{COLORAMA_COLORS['BOLD']}📊 FINAL CLASSIFICATION{COLORAMA_COLORS['ENDC']}",
    'recommendations': f"{COLORAMA_COLORS['YELLOW']}💡 RECOMMENDATIONS{COLORAMA_COLORS['ENDC']}",
    'report_generation': f"{COLORAMA_COLORS['BLUE']}📋 REPORT GENERATION{COLORAMA_COLORS['ENDC']}"
}

# Classification-specific colors
CLASSIFICATION_COLORS = {
    'Pathogenic': f"{COLORAMA_COLORS['RED']}{COLORAMA_COLORS['BOLD']}",
    'Likely Pathogenic': f"{COLORAMA_COLORS['RED']}",
    'Uncertain Significance': f"{COLORAMA_COLORS['YELLOW']}",
    'Likely Benign': f"{COLORAMA_COLORS['GREEN']}",
    'Benign': f"{COLORAMA_COLORS['GREEN']}{COLORAMA_COLORS['BOLD']}"
}

# Error and message constants
ERROR_MESSAGES = {
    'invalid_gene': 'Invalid gene symbol provided',
    'invalid_chromosome': 'Invalid chromosome format',
    'invalid_position': 'Invalid genomic position',
    'invalid_allele': 'Invalid allele format',
    'missing_data': 'Required data is missing',
    'api_error': 'API request failed',
    'file_not_found': 'File not found',
    'invalid_format': 'Invalid input format'
}

SUCCESS_MESSAGES = {
    'data_loaded': 'Data loaded successfully',
    'analysis_complete': 'Analysis completed successfully',
    'report_generated': 'Report generated successfully'
}

WARNING_MESSAGES = {
    'missing_optional_data': 'Optional data missing, proceeding with available data',
    'low_confidence': 'Low confidence result due to limited data'
}

INFO_MESSAGES = {
    'interactive_mode': 'Interactive mode enabled - user input required',
    'test_mode': 'Test mode enabled - using default responses'
}
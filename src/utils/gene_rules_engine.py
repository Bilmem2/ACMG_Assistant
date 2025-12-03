# =============================================================================
# Gene-Specific Rules Engine (Utils)
# =============================================================================
# Changes in this documentation pass (Dec 2024):
# - Added clarification about relationship with core/gene_rules_engine.py
# - Enhanced function docstrings with detailed Args/Returns sections
# - Added type hints throughout
#
# NOTE: This is the FUNCTION-BASED implementation for ClinGen threshold calculations.
# A CLASS-BASED implementation also exists in core/gene_rules_engine.py which
# handles rule-based evidence modifications.
# TODO: Consider merging these two modules in a future refactor.
# =============================================================================
"""
Gene-Specific Rules Engine
===========================

Dynamically generates gene-specific thresholds based on ClinGen Dosage Sensitivity data.
Integrates haploinsufficiency (HI) scores with ACMG classification thresholds.

This module provides FUNCTION-BASED utilities for:
- Calculating gene-specific BA1/BS1/PM2 thresholds based on HI scores
- Building comprehensive threshold dictionaries from ClinGen TSV data
- Generating human-readable threshold explanations

See also: core/gene_rules_engine.py for class-based rule application.
"""

from typing import Dict, Any, Optional
from config.constants import GENE_SPECIFIC_THRESHOLDS

def get_gene_specific_thresholds(gene: str, hi_score: Optional[int] = None, 
                                   clingen_dosage_data: Optional[Dict] = None) -> Dict[str, float]:
    """
    Get gene-specific thresholds with ClinGen Dosage Sensitivity integration.
    
    **Logic:**
    - HI=3 (Sufficient evidence): Stricter thresholds (lower BA1/BS1, higher PM2 sensitivity)
    - HI=2 (Some evidence): Moderately strict thresholds
    - HI=1 (Little evidence): Relaxed thresholds
    - HI=0/30/40 (No/unlikely/AR): Default thresholds
    
    Args:
        gene (str): Gene symbol (e.g., 'BRCA1')
        hi_score (int, optional): Haploinsufficiency score (0-3, 30, 40)
        clingen_dosage_data (dict, optional): Full ClinGen dosage data
    
    Returns:
        Dict with BA1, BS1, PM2 thresholds
    """
    # Start with manual gene-specific thresholds if available
    if gene in GENE_SPECIFIC_THRESHOLDS:
        thresholds = GENE_SPECIFIC_THRESHOLDS[gene].copy()
    else:
        thresholds = GENE_SPECIFIC_THRESHOLDS['default'].copy()
    
    # Apply ClinGen Dosage-based adjustments
    if hi_score is not None:
        if hi_score == 3:
            # Sufficient evidence for HI: STRICTER thresholds
            # Rationale: Haploinsufficient genes have LOF pathogenic variants at lower frequencies
            thresholds['BA1'] = min(thresholds.get('BA1', 0.05), 0.03)  # Reduce to 3% from 5%
            thresholds['BS1'] = min(thresholds.get('BS1', 0.01), 0.005)  # Reduce to 0.5% from 1%
            thresholds['PM2'] = min(thresholds.get('PM2', 0.0001), 0.00005)  # More sensitive absence threshold
            thresholds['rationale'] = f"{gene} HI Score=3 (Sufficient): Strict thresholds for haploinsufficient gene"
        
        elif hi_score == 2:
            # Some evidence for HI: MODERATELY STRICT thresholds
            thresholds['BA1'] = min(thresholds.get('BA1', 0.05), 0.04)  # Reduce to 4%
            thresholds['BS1'] = min(thresholds.get('BS1', 0.01), 0.0075)  # Reduce to 0.75%
            thresholds['PM2'] = min(thresholds.get('PM2', 0.0001), 0.00008)
            thresholds['rationale'] = f"{gene} HI Score=2 (Some): Moderately strict thresholds"
        
        elif hi_score == 1:
            # Little evidence for HI: RELAXED thresholds
            thresholds['BA1'] = max(thresholds.get('BA1', 0.05), 0.06)  # Increase to 6%
            thresholds['BS1'] = max(thresholds.get('BS1', 0.01), 0.015)  # Increase to 1.5%
            thresholds['PM2'] = max(thresholds.get('PM2', 0.0001), 0.0002)
            thresholds['rationale'] = f"{gene} HI Score=1 (Little): Relaxed thresholds"
        
        elif hi_score == 0:
            # No evidence for HI: DEFAULT or slightly relaxed
            thresholds['BA1'] = max(thresholds.get('BA1', 0.05), 0.07)  # Increase to 7%
            thresholds['BS1'] = max(thresholds.get('BS1', 0.01), 0.02)  # Increase to 2%
            thresholds['PM2'] = max(thresholds.get('PM2', 0.0001), 0.0003)
            thresholds['rationale'] = f"{gene} HI Score=0 (No evidence): Default/relaxed thresholds"
        
        elif hi_score == 30:
            # AR phenotype: DEFAULT thresholds (not haploinsufficient)
            thresholds['rationale'] = f"{gene} HI Score=30 (AR phenotype): Default thresholds"
        
        elif hi_score == 40:
            # Unlikely HI: VERY RELAXED thresholds
            thresholds['BA1'] = max(thresholds.get('BA1', 0.05), 0.1)  # Increase to 10%
            thresholds['BS1'] = max(thresholds.get('BS1', 0.01), 0.03)  # Increase to 3%
            thresholds['PM2'] = max(thresholds.get('PM2', 0.0001), 0.0005)
            thresholds['rationale'] = f"{gene} HI Score=40 (Unlikely HI): Very relaxed thresholds"
    
    # Add ClinGen metadata if available
    if clingen_dosage_data:
        thresholds['clingen_data'] = {
            'hi_score': clingen_dosage_data.get('haploinsufficiency_score'),
            'hi_description': clingen_dosage_data.get('haploinsufficiency_description'),
            'ts_score': clingen_dosage_data.get('triplosensitivity_score'),
            'source': clingen_dosage_data.get('source'),
            'confidence': clingen_dosage_data.get('confidence')
        }
    
    return thresholds


def build_gene_specific_thresholds_from_clingen(clingen_tsv_path: str = None) -> Dict[str, Dict[str, float]]:
    """
    Build comprehensive gene-specific thresholds from ClinGen TSV file.
    
    This function can be run once to generate a comprehensive GENE_SPECIFIC_THRESHOLDS
    dictionary based on all ClinGen dosage sensitivity data.
    
    Args:
        clingen_tsv_path (str, optional): Path to ClinGen TSV file
    
    Returns:
        Dict mapping gene symbols to threshold dictionaries
    """
    import os
    
    if clingen_tsv_path is None:
        # Try to find local TSV file
        possible_paths = [
            'ClinGen_gene_curation_list_GRCh38.tsv',
            '../ClinGen_gene_curation_list_GRCh38.tsv',
            '../../ClinGen_gene_curation_list_GRCh38.tsv'
        ]
        for path in possible_paths:
            if os.path.exists(path):
                clingen_tsv_path = path
                break
    
    if clingen_tsv_path is None or not os.path.exists(clingen_tsv_path):
        print(f"⚠️  ClinGen TSV file not found. Using default thresholds.")
        return GENE_SPECIFIC_THRESHOLDS.copy()
    
    # Parse ClinGen TSV
    gene_thresholds = {}
    
    try:
        with open(clingen_tsv_path, 'r', encoding='utf-8') as f:
            lines = f.readlines()
        
        # Skip header
        for line in lines[1:]:
            parts = line.strip().split('\t')
            if len(parts) < 5:
                continue
            
            gene_symbol = parts[0].strip()
            hi_score_str = parts[4].strip()
            
            # Parse HI score
            try:
                hi_score = int(hi_score_str) if hi_score_str.isdigit() else None
            except:
                hi_score = None
            
            if hi_score is not None and gene_symbol:
                # Generate thresholds based on HI score
                thresholds = get_gene_specific_thresholds(gene_symbol, hi_score=hi_score)
                gene_thresholds[gene_symbol] = thresholds
        
        print(f"✅ Generated thresholds for {len(gene_thresholds)} genes from ClinGen data")
        
        # Add default
        gene_thresholds['default'] = GENE_SPECIFIC_THRESHOLDS['default'].copy()
        
        return gene_thresholds
    
    except Exception as e:
        print(f"❌ Error parsing ClinGen TSV: {e}")
        return GENE_SPECIFIC_THRESHOLDS.copy()


def get_threshold_explanation(gene: str, hi_score: Optional[int] = None) -> str:
    """
    Get human-readable explanation of why specific thresholds are applied.
    
    Args:
        gene (str): Gene symbol
        hi_score (int, optional): Haploinsufficiency score
    
    Returns:
        str: Explanation text
    """
    thresholds = get_gene_specific_thresholds(gene, hi_score=hi_score)
    
    explanation = f"\n📊 **Gene-Specific Thresholds for {gene}**:\n"
    explanation += f"  BA1 (Benign Stand-alone): AF > {thresholds['BA1']:.4f} ({thresholds['BA1']*100:.2f}%)\n"
    explanation += f"  BS1 (Benign Strong): AF > {thresholds['BS1']:.4f} ({thresholds['BS1']*100:.2f}%)\n"
    explanation += f"  PM2 (Absent/Rare): AF < {thresholds['PM2']:.6f} ({thresholds['PM2']*100:.4f}%)\n"
    
    if 'rationale' in thresholds:
        explanation += f"\n  📝 Rationale: {thresholds['rationale']}\n"
    
    if 'clingen_data' in thresholds:
        clingen = thresholds['clingen_data']
        explanation += f"\n  🔬 ClinGen Data:\n"
        explanation += f"    - HI Score: {clingen.get('hi_score')}\n"
        explanation += f"    - HI Description: {clingen.get('hi_description')}\n"
        explanation += f"    - TS Score: {clingen.get('ts_score')}\n"
        explanation += f"    - Confidence: {clingen.get('confidence')}\n"
    
    return explanation


# Example usage:
if __name__ == "__main__":
    from colorama import init, Fore, Style
    init()
    
    print(f"\n{Fore.CYAN}{'='*80}{Style.RESET_ALL}")
    print(f"{Fore.CYAN}GENE-SPECIFIC RULES ENGINE TEST{Style.RESET_ALL}")
    print(f"{Fore.CYAN}{'='*80}{Style.RESET_ALL}\n")
    
    # Test with different HI scores
    test_cases = [
        ("BRCA1", 3, "Sufficient HI evidence"),
        ("TP53", 3, "Sufficient HI evidence"),
        ("TTN", 1, "Little HI evidence"),
        ("MUC16", 0, "No HI evidence"),
        ("CFTR", 30, "AR phenotype"),
        ("UNKNOWN_GENE", None, "No ClinGen data")
    ]
    
    for gene, hi_score, description in test_cases:
        print(f"{Fore.YELLOW}{gene} (HI={hi_score}): {description}{Style.RESET_ALL}")
        explanation = get_threshold_explanation(gene, hi_score)
        print(explanation)
    
    # Test building comprehensive thresholds
    print(f"\n{Fore.GREEN}{'='*80}{Style.RESET_ALL}")
    print(f"{Fore.GREEN}Building Comprehensive Thresholds from ClinGen TSV{Style.RESET_ALL}")
    print(f"{Fore.GREEN}{'='*80}{Style.RESET_ALL}\n")
    
    all_thresholds = build_gene_specific_thresholds_from_clingen()
    print(f"\nTotal genes with custom thresholds: {len(all_thresholds)}")
    print(f"Sample genes: {list(all_thresholds.keys())[:10]}")

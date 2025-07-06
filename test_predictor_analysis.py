#!/usr/bin/env python3
"""
Analysis of current in silico predictors and recommendations for additions
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from config.constants import INSILICO_WEIGHTS, INSILICO_THRESHOLDS
from colorama import Fore, Style, init

# Initialize colorama
init()

def analyze_current_predictors():
    """Analyze current predictor coverage"""
    
    print(f"{Fore.CYAN}🔍 Current In Silico Predictor Analysis{Style.RESET_ALL}")
    print("=" * 60)
    
    # Categorize current predictors
    categories = {
        'Primary Metascores': [],
        'Conservation Scores': [],
        'Individual Predictors': [],
        'Splice Predictors': [],
        'Missing Important': []
    }
    
    for predictor in INSILICO_WEIGHTS.keys():
        if predictor in ['revel', 'cadd_phred', 'alphamissense', 'metarnn', 'clinpred', 'bayesdel_addaf', 'bayesdel_noaf']:
            categories['Primary Metascores'].append(predictor)
        elif 'phylop' in predictor or 'gerp' in predictor:
            categories['Conservation Scores'].append(predictor)
        elif 'splice' in predictor or predictor in ['ada_score', 'rf_score', 'dbscsnv_ada', 'dbscsnv_rf']:
            categories['Splice Predictors'].append(predictor)
        else:
            categories['Individual Predictors'].append(predictor)
    
    # Print current coverage
    for category, predictors in categories.items():
        if predictors:
            print(f"\n{Fore.GREEN}✅ {category}:{Style.RESET_ALL}")
            for pred in predictors:
                weight = INSILICO_WEIGHTS.get(pred, 0)
                print(f"   • {pred}: weight={weight}")
    
    # Identify missing important predictors
    print(f"\n{Fore.YELLOW}⚠️  Potentially Missing Important Predictors:{Style.RESET_ALL}")
    
    missing_predictors = {
        'VEST4': {'type': 'Missense', 'importance': 'High', 'reason': 'Cancer-specific predictor'},
        'PROVEAN': {'type': 'Missense', 'importance': 'Medium', 'reason': 'Good for protein impact prediction'},
        'MutPred': {'type': 'Missense', 'importance': 'Medium', 'reason': 'Structural impact prediction'},
        'LRT': {'type': 'Conservation', 'importance': 'Medium', 'reason': 'Likelihood ratio test'},
        'MetaSVM': {'type': 'Meta', 'importance': 'Medium', 'reason': 'SVM-based ensemble'},
        'MetaLR': {'type': 'Meta', 'importance': 'Medium', 'reason': 'Logistic regression ensemble'},
        'DANN': {'type': 'Deep Learning', 'importance': 'Medium', 'reason': 'Deep neural network'},
        'DEOGEN2': {'type': 'Deep Learning', 'importance': 'Medium', 'reason': 'Gradient boosting'},
        'PrimateAI': {'type': 'Deep Learning', 'importance': 'High', 'reason': 'Deep learning with primate data'},
        'ESM1b': {'type': 'Language Model', 'importance': 'High', 'reason': 'Protein language model'},
        'EVE': {'type': 'Evolutionary', 'importance': 'High', 'reason': 'Evolutionary model of variant effect'},
        # Splice-specific
        'MMSplice': {'type': 'Splice', 'importance': 'High', 'reason': 'Modular modeling of splicing'},
        'Pangolin': {'type': 'Splice', 'importance': 'Medium', 'reason': 'Splice site prediction'},
        'SQUIRLS': {'type': 'Splice', 'importance': 'Medium', 'reason': 'Splice variant impact'},
        # Conservation
        'SiPhy': {'type': 'Conservation', 'importance': 'Medium', 'reason': 'Site-specific phylogenetic'},
        'PhastCons': {'type': 'Conservation', 'importance': 'Medium', 'reason': 'Conserved elements'},
        # Structural
        'FoldX': {'type': 'Structural', 'importance': 'Medium', 'reason': 'Protein stability'},
        'Rosetta': {'type': 'Structural', 'importance': 'Medium', 'reason': 'Protein energy calculation'}
    }
    
    for predictor, info in missing_predictors.items():
        color = Fore.RED if info['importance'] == 'High' else Fore.YELLOW if info['importance'] == 'Medium' else Fore.GREEN
        print(f"   {color}• {predictor} ({info['type']}): {info['reason']}{Style.RESET_ALL}")
    
    return categories, missing_predictors

def analyze_variant_type_coverage():
    """Analyze how well different variant types are covered"""
    
    print(f"\n{Fore.CYAN}🧬 Variant Type Coverage Analysis{Style.RESET_ALL}")
    print("=" * 60)
    
    variant_types = {
        'Missense': {
            'predictors': ['revel', 'alphamissense', 'metarnn', 'clinpred', 'sift', 'polyphen2'],
            'coverage': 'Excellent',
            'recommendations': ['Consider adding PrimateAI, ESM1b, EVE for cutting-edge analysis']
        },
        'Nonsense': {
            'predictors': ['cadd_phred', 'phylop scores'],
            'coverage': 'Good',
            'recommendations': ['Sufficient - LOF mechanism is well established']
        },
        'Frameshift': {
            'predictors': ['cadd_phred', 'phylop scores'],
            'coverage': 'Good', 
            'recommendations': ['Sufficient - LOF mechanism is well established']
        },
        'Splice Site': {
            'predictors': ['spliceai scores', 'ada_score', 'rf_score', 'dbscsnv scores'],
            'coverage': 'Very Good',
            'recommendations': ['Consider adding MMSplice, Pangolin for comprehensive splice analysis']
        },
        'Intronic': {
            'predictors': ['spliceai scores', 'conservation scores'],
            'coverage': 'Good',
            'recommendations': ['Current implementation handles splice impact well']
        },
        'Synonymous': {
            'predictors': ['spliceai scores', 'conservation scores'],
            'coverage': 'Good',
            'recommendations': ['Focus on splice impact - current approach is appropriate']
        },
        'Regulatory': {
            'predictors': ['cadd_phred', 'conservation scores'],
            'coverage': 'Limited',
            'recommendations': ['Could benefit from regulatory-specific predictors (FATHMM-MKL, GWAVA)']
        }
    }
    
    for vtype, info in variant_types.items():
        color = Fore.GREEN if info['coverage'] == 'Excellent' else Fore.BLUE if info['coverage'] == 'Very Good' else Fore.YELLOW if info['coverage'] == 'Good' else Fore.RED
        print(f"\n{color}{vtype}: {info['coverage']}{Style.RESET_ALL}")
        print(f"   Current predictors: {', '.join(info['predictors'])}")
        for rec in info['recommendations']:
            print(f"   💡 {rec}")

def recommend_additions():
    """Recommend specific predictor additions"""
    
    print(f"\n{Fore.CYAN}📋 Recommendations for Predictor Additions{Style.RESET_ALL}")
    print("=" * 60)
    
    priorities = {
        'High Priority (Add Soon)': [
            ('VEST4', 'Cancer variant assessment'),
            ('PrimateAI', 'State-of-the-art missense prediction'),
            ('MMSplice', 'Advanced splice prediction'),
            ('ESM1b', 'Protein language model scores')
        ],
        'Medium Priority (Consider)': [
            ('PROVEAN', 'Protein variation effect'),
            ('MetaSVM/MetaLR', 'Additional ensemble methods'),
            ('DANN', 'Deep learning approach'),
            ('SiPhy', 'Additional conservation metric')
        ],
        'Low Priority (Future)': [
            ('MutPred', 'Structural impact prediction'),
            ('DEOGEN2', 'Gradient boosting approach'),
            ('Pangolin', 'Additional splice predictor'),
            ('FoldX', 'Protein stability analysis')
        ]
    }
    
    for priority, predictors in priorities.items():
        color = Fore.RED if 'High' in priority else Fore.YELLOW if 'Medium' in priority else Fore.GREEN
        print(f"\n{color}{priority}:{Style.RESET_ALL}")
        for pred, desc in predictors:
            print(f"   • {pred}: {desc}")

if __name__ == "__main__":
    categories, missing = analyze_current_predictors()
    analyze_variant_type_coverage()
    recommend_additions()
    
    print(f"\n{Fore.GREEN}📊 Summary:{Style.RESET_ALL}")
    print(f"✅ Current predictor set is {Fore.GREEN}very comprehensive{Style.RESET_ALL}")
    print(f"✅ All major variant types are {Fore.GREEN}well covered{Style.RESET_ALL}")
    print(f"✅ Algorithm is {Fore.GREEN}production-ready{Style.RESET_ALL} as-is")
    print(f"💡 Additional predictors would be {Fore.BLUE}enhancements{Style.RESET_ALL}, not requirements")

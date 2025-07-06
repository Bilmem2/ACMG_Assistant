#!/usr/bin/env python3
"""
ACMG Variant Classification Assistant
====================================

A comprehensive tool for classifying genetic variants according to ACMG/AMP 2015 and 2023 guidelines.

Author: Can Sevilmiş
License: MIT
Last Updated: July 6, 2025

Features:
- ACMG/AMP 2015 and 2023 guidelines implementation
- Support for all major variant types
- Integration with multiple in silico predictors
- Population frequency analysis
- Statistical frameworks (Fisher's Exact Test, VAMPP-score-like metascore)
- Interactive command-line interface
- Automatic ClinVar and Ensembl API integration
- Comprehensive reporting and validation
"""

import sys
import os
import argparse
import json
from typing import Dict, List, Optional, Tuple, Any
from datetime import datetime
from colorama import init, Fore, Back, Style

# Initialize colorama for cross-platform colored output
init(autoreset=True)

# Import custom modules
from core.variant_data import VariantData
from core.acmg_classifier import ACMGClassifier
from core.evidence_evaluator import EvidenceEvaluator
from utils.input_handler import InputHandler
from utils.api_client import APIClient
from utils.validators import Validators
from utils.report_generator import ReportGenerator
from config.constants import (
    VERSION_INFO, COLORAMA_COLORS, SECTION_HEADERS, 
    CLASSIFICATION_COLORS, get_colored_message
)
from utils.report_generator import ReportGenerator
from config.constants import *

class ACMGAssistant:
    """Main ACMG Classification Assistant class."""
    
    def __init__(self, use_2023_guidelines: bool = False, test_mode: bool = False):
        """
        Initialize the ACMG Assistant.
        
        Args:
            use_2023_guidelines (bool): Whether to use ACMG 2023 guidelines
            test_mode (bool): Whether to run in test mode
        """
        self.use_2023_guidelines = use_2023_guidelines
        self.test_mode = test_mode
        self.api_client = APIClient()
        self.input_handler = InputHandler(test_mode=test_mode, api_client=self.api_client)
        self.validators = Validators()
        self.evidence_evaluator = EvidenceEvaluator(use_2023_guidelines)
        self.classifier = ACMGClassifier(use_2023_guidelines)
        self.report_generator = ReportGenerator()
        
        # Print welcome message
        self._print_welcome_message()
    
    def _print_welcome_message(self):
        """Print welcome message and guidelines information."""
        print(f"\n{Fore.CYAN}{'='*80}")
        print(f"{Fore.CYAN}{Style.BRIGHT}🧬 ACMG Variant Classification Assistant")
        print(f"{Fore.CYAN}{'='*80}")
        print(f"{Fore.YELLOW}Guidelines: {Fore.WHITE}{Style.BRIGHT}ACMG/AMP {'2023' if self.use_2023_guidelines else '2015'}")
        print(f"{Fore.YELLOW}Mode: {Fore.WHITE}{Style.BRIGHT}{'Test Mode' if self.test_mode else 'Standard Mode'}")
        print(f"{Fore.YELLOW}Date: {Fore.WHITE}{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print(f"{Fore.CYAN}{'='*80}")
        
        if not self.test_mode:
            print(f"\n{Fore.RED}{Style.BRIGHT}⚠️  DISCLAIMER:")
            print(f"{Fore.RED}This tool is for research and educational purposes only.")
            print(f"{Fore.RED}Clinical variant interpretations must be validated by certified professionals.")
            print(f"{Fore.RED}Always cross-check results with primary sources (ClinVar, Varsome).")
            print(f"{Fore.CYAN}{'='*80}")
    
    def run(self):
        """Main execution method."""
        while True:
            try:
                # Collect variant data
                print(f"\n{SECTION_HEADERS['basic_info']}")
                print(f"{Fore.BLUE}{'-'*40}{Style.RESET_ALL}")
                variant_data = self._collect_variant_data()
                
                # Evaluate evidence
                print(f"\n{COLORAMA_COLORS['CYAN']}{COLORAMA_COLORS['BOLD']}🔍 EVIDENCE EVALUATION{COLORAMA_COLORS['RESET']}")
                print(f"{Fore.BLUE}{'-'*40}{Style.RESET_ALL}")
                evidence_results = self._evaluate_evidence(variant_data)
                
                # Classify variant
                print(f"\n{COLORAMA_COLORS['MAGENTA']}{COLORAMA_COLORS['BOLD']}⚖️  VARIANT CLASSIFICATION{COLORAMA_COLORS['RESET']}")
                print(f"{Fore.BLUE}{'-'*40}{Style.RESET_ALL}")
                classification_result = self._classify_variant(evidence_results)
                
                # Generate report
                print(f"\n{SECTION_HEADERS['report_generation']}")
                print(f"{Fore.BLUE}{'-'*40}{Style.RESET_ALL}")
                report_path = self._generate_report(variant_data, evidence_results, classification_result)
                
                # Display results
                self._display_results(classification_result, report_path, variant_data)
                
                # Ask if user wants to analyze another variant
                if not self.test_mode:
                    self._ask_for_restart()
                else:
                    break
                    
            except KeyboardInterrupt:
                print(f"\n\n{Fore.RED}❌ Operation cancelled by user.{Style.RESET_ALL}")
                self._pause_before_exit()
                sys.exit(1)
            except Exception as e:
                print(f"\n❌ An error occurred: {str(e)}")
                print(f"Error type: {type(e).__name__}")
                import traceback
                print("Full traceback:")
                traceback.print_exc()
                
                if not self.test_mode:
                    print("\nPlease check your inputs and try again.")
                    self._ask_for_restart()
                else:
                    self._pause_before_exit()
                    sys.exit(1)
    
    def _ask_for_restart(self):
        """Ask user if they want to analyze another variant."""
        print(f"\n{COLORAMA_COLORS['CYAN']}{'='*80}{COLORAMA_COLORS['RESET']}")
        print(f"{COLORAMA_COLORS['GREEN']}{COLORAMA_COLORS['BOLD']}🔄 ANALYSIS COMPLETE{COLORAMA_COLORS['RESET']}")
        print(f"{COLORAMA_COLORS['CYAN']}{'='*80}{COLORAMA_COLORS['RESET']}")
        
        while True:
            try:
                choice = input(f"\n{COLORAMA_COLORS['YELLOW']}Would you like to analyze another variant? (y/n) or press Enter to exit: {COLORAMA_COLORS['RESET']}").strip().lower()
                
                if choice == '' or choice == 'n' or choice == 'no':
                    print(f"\n{COLORAMA_COLORS['GREEN']}👋 Thank you for using ACMG Variant Classification Assistant!{COLORAMA_COLORS['RESET']}")
                    self._pause_before_exit()
                    sys.exit(0)
                elif choice == 'y' or choice == 'yes':
                    print(f"\n{COLORAMA_COLORS['GREEN']}🔄 Starting new analysis...{COLORAMA_COLORS['RESET']}")
                    return
                else:
                    print(f"{COLORAMA_COLORS['RED']}❌ Please enter 'y' for yes, 'n' for no, or press Enter to exit.{COLORAMA_COLORS['RESET']}")
                    
            except KeyboardInterrupt:
                print(f"\n\n{COLORAMA_COLORS['RED']}❌ Operation cancelled by user.{COLORAMA_COLORS['RESET']}")
                self._pause_before_exit()
                sys.exit(1)
    
    def _pause_before_exit(self):
        """Pause before exiting to allow user to see any error messages."""
        try:
            print(f"{COLORAMA_COLORS['CYAN']}Press Enter to exit...{COLORAMA_COLORS['RESET']}")
            input()
        except KeyboardInterrupt:
            pass
    
    def _collect_variant_data(self) -> VariantData:
        """Collect all variant data from user input."""
        print("Please provide the following information:")
        print("(All in silico scores must be entered manually)")
        
        # Basic variant information
        basic_info = self.input_handler.collect_basic_info()
        
        # Population data
        population_data = self.input_handler.collect_population_data()
        
        # In silico scores
        insilico_data = self.input_handler.collect_insilico_data(basic_info['variant_type'])
        
        # Genetic data
        genetic_data = self.input_handler.collect_genetic_data()
        
        # Functional data
        functional_data = self.input_handler.collect_functional_data()
        
        # Create variant data object
        variant_data = VariantData(
            basic_info=basic_info,
            population_data=population_data,
            insilico_data=insilico_data,
            genetic_data=genetic_data,
            functional_data=functional_data
        )
        
        # Enrich with API data
        if not self.test_mode:
            self._enrich_with_api_data(variant_data)
        
        return variant_data
    
    def _enrich_with_api_data(self, variant_data: VariantData):
        """Enrich variant data with API information."""
        print("\n🌐 Fetching additional data from APIs...")
        
        # Get ClinVar status
        print("- Checking ClinVar status...")
        clinvar_data = self.api_client.get_clinvar_status(
            variant_data.basic_info['chromosome'],
            variant_data.basic_info['position'],
            variant_data.basic_info['ref_allele'],
            variant_data.basic_info['alt_allele']
        )
        variant_data.clinvar_data = clinvar_data
    
    def _evaluate_evidence(self, variant_data: VariantData) -> Dict[str, Any]:
        """Evaluate ACMG evidence criteria."""
        return self.evidence_evaluator.evaluate_all_criteria(variant_data)
    
    def _classify_variant(self, evidence_results: Dict[str, Any]) -> Dict[str, Any]:
        """Classify variant based on evidence."""
        return self.classifier.classify(evidence_results)
    
    def _generate_report(self, variant_data: VariantData, evidence_results: Dict[str, Any], 
                        classification_result: Dict[str, Any]) -> str:
        """Generate detailed classification report."""
        return self.report_generator.generate_report(
            variant_data, evidence_results, classification_result
        )
    
    def _display_results(self, classification_result: Dict[str, Any], report_path: str, variant_data=None):
        """Display classification results to user."""
        print("\n" + "="*80)
        print("🎯 CLASSIFICATION RESULTS")
        print("="*80)
        
        # Main classification
        classification = classification_result['classification']
        confidence = classification_result.get('confidence', 'N/A')
        
        # Display classification result with color
        classification = classification_result['classification']
        confidence = classification_result['confidence']
        
        # Get color for classification
        if classification in CLASSIFICATION_COLORS:
            color_info = CLASSIFICATION_COLORS[classification]
            color = COLORAMA_COLORS[color_info[0]]
            style = COLORAMA_COLORS[color_info[1]] if color_info[1] else ""
            print(f"\n{color}{style}📊 FINAL CLASSIFICATION: {classification}{COLORAMA_COLORS['RESET']}")
        else:
            print(f"\n{COLORAMA_COLORS['MAGENTA']}📊 FINAL CLASSIFICATION: {classification}{COLORAMA_COLORS['RESET']}")
        
        print(f"{COLORAMA_COLORS['CYAN']}🔍 CONFIDENCE LEVEL: {confidence}{COLORAMA_COLORS['RESET']}")
        
        # Applied criteria
        print(f"\n{COLORAMA_COLORS['BLUE']}{COLORAMA_COLORS['BOLD']}📋 APPLIED CRITERIA:{COLORAMA_COLORS['RESET']}")
        for category, criteria in classification_result['applied_criteria'].items():
            if criteria:
                print(f"  {COLORAMA_COLORS['WHITE']}{category}: {', '.join(criteria)}{COLORAMA_COLORS['RESET']}")
        
        # Suggestions
        if classification_result.get('suggestions'):
            print(f"\n{COLORAMA_COLORS['YELLOW']}{COLORAMA_COLORS['BOLD']}💡 SUGGESTIONS:{COLORAMA_COLORS['RESET']}")
            for suggestion in classification_result['suggestions']:
                print(f"  {COLORAMA_COLORS['YELLOW']}• {suggestion}{COLORAMA_COLORS['RESET']}")
        
        # Report location
        print(f"\n{COLORAMA_COLORS['GREEN']}📄 DETAILED REPORT: {report_path}{COLORAMA_COLORS['RESET']}")
        
        # Varsome URL
        if variant_data:
            varsome_url = self.api_client.generate_varsome_url(
                variant_data.basic_info['chromosome'],
                variant_data.basic_info['position'],
                variant_data.basic_info['ref_allele'],
                variant_data.basic_info['alt_allele']
            )
            if varsome_url:
                print(f"{COLORAMA_COLORS['CYAN']}🔗 VARSOME URL: {varsome_url}{COLORAMA_COLORS['RESET']}")
        
        print(f"{COLORAMA_COLORS['CYAN']}{'='*80}{COLORAMA_COLORS['RESET']}")


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="ACMG Variant Classification Assistant",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python acmg_assistant.py                    # Standard mode with ACMG 2015 guidelines
  python acmg_assistant.py --acmg-2023        # Use ACMG 2023 guidelines
  python acmg_assistant.py --test             # Run in test mode
  python acmg_assistant.py --test --acmg-2023 # Test mode with 2023 guidelines
        """
    )
    
    parser.add_argument(
        '--acmg-2023',
        action='store_true',
        help='Use ACMG/AMP 2023 guidelines instead of 2015'
    )
    
    parser.add_argument(
        '--test',
        action='store_true',
        help='Run in test mode with sample data'
    )
    
    parser.add_argument(
        '--version',
        action='version',
        version='ACMG Assistant v1.0.0'
    )
    
    args = parser.parse_args()
    
    # Create and run assistant
    assistant = ACMGAssistant(
        use_2023_guidelines=args.acmg_2023,
        test_mode=args.test
    )
    
    assistant.run()

if __name__ == "__main__":
    main()

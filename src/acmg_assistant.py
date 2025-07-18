#!/usr/bin/env python3
"""
ACMG Variant Classification Assistant
====================================

A comprehensive tool for classifying genetic variants according to ACMG/AMP 2015 and 2023 guidelines.

Author: Can Sevilmiş
License: MIT
Last Updated: July 10, 2025

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
        self.evidence_evaluator = EvidenceEvaluator(use_2023_guidelines, test_mode=test_mode)
        self.classifier = ACMGClassifier(use_2023_guidelines)
        
        # Get the correct output directory for reports
        if getattr(sys, 'frozen', False):
            # Running as executable - use the directory where the executable is located
            output_dir = os.path.dirname(sys.executable)
        else:
            # Running as script - use the project root directory (parent of src)
            current_dir = os.path.dirname(os.path.abspath(__file__))
            output_dir = os.path.dirname(current_dir)
        
        self.report_generator = ReportGenerator(output_dir=output_dir)
        
        # Print welcome message (show disclaimer)
        self._print_welcome_message(show_disclaimer=True)
    
    def _print_welcome_message(self, show_disclaimer=False):
        """Print English welcome message and interactive mode selection info. Optionally show disclaimer."""
        print(f"\n{COLORAMA_COLORS['CYAN']}{'='*80}")
        print(f"{COLORAMA_COLORS['CYAN']}{COLORAMA_COLORS['BOLD']}🧬 ACMG Variant Classification Assistant v{VERSION_INFO['version']}")
        print(f"{COLORAMA_COLORS['CYAN']}{'='*80}")
        print(f"{COLORAMA_COLORS['YELLOW']}Available Modes:{COLORAMA_COLORS['RESET']}")
        print(f"{COLORAMA_COLORS['WHITE']}  • ACMG/AMP 2015 (default)        →  Starts in standard mode")
        print(f"{COLORAMA_COLORS['WHITE']}  • ACMG/AMP 2023                  →  Start with '--acmg-2023' or select in exe")
        print(f"{COLORAMA_COLORS['WHITE']}  • Batch Mode (Multiple Variants)  →  Start with '--batch file.csv'")
        print(f"{COLORAMA_COLORS['YELLOW']}Active Mode: {COLORAMA_COLORS['WHITE']}{COLORAMA_COLORS['BOLD']}ACMG/AMP {'2023' if self.use_2023_guidelines else '2015'}{COLORAMA_COLORS['RESET']}")
        print(f"{COLORAMA_COLORS['YELLOW']}Run Type: {COLORAMA_COLORS['WHITE']}{COLORAMA_COLORS['BOLD']}{'Test Mode' if self.test_mode else 'Standard Mode'}{COLORAMA_COLORS['RESET']}")
        print(f"{COLORAMA_COLORS['YELLOW']}Date: {COLORAMA_COLORS['WHITE']}{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}{COLORAMA_COLORS['RESET']}")
        print(f"{COLORAMA_COLORS['CYAN']}{'='*80}")
        print(f"{COLORAMA_COLORS['GREEN']}Created by: {VERSION_INFO['author']}{COLORAMA_COLORS['RESET']}")

        print(f"\n{COLORAMA_COLORS['YELLOW']}Mode Selection Guide:{COLORAMA_COLORS['RESET']}")
        print(f"{COLORAMA_COLORS['WHITE']}  • Press Enter to start in default mode (ACMG/AMP 2015).{COLORAMA_COLORS['RESET']}")
        print(f"{COLORAMA_COLORS['WHITE']}  • Or type one of the following commands to switch mode:{COLORAMA_COLORS['RESET']}")
        print(f"{COLORAMA_COLORS['WHITE']}      --acmg-2023   → ACMG/AMP 2023 mode{COLORAMA_COLORS['RESET']}")
        print(f"{COLORAMA_COLORS['WHITE']}      --batch file.csv → Batch mode (CSV file required){COLORAMA_COLORS['RESET']}")
        print(f"{COLORAMA_COLORS['WHITE']}  • In exe, you can select mode from the main screen.{COLORAMA_COLORS['RESET']}")


        # Show disclaimer only on initial welcome page
        if show_disclaimer:
            print(f"\n{COLORAMA_COLORS['RED']}{COLORAMA_COLORS['BOLD']}⚠️  DISCLAIMER:{COLORAMA_COLORS['RESET']}")
            print(f"{COLORAMA_COLORS['RED']}This tool is for research and educational purposes only.{COLORAMA_COLORS['RESET']}")
            print(f"{COLORAMA_COLORS['RED']}Clinical variant interpretations must be validated by certified professionals.{COLORAMA_COLORS['RESET']}")
            print(f"{COLORAMA_COLORS['RED']}Always cross-check results with primary sources (ClinVar, Varsome).{COLORAMA_COLORS['RESET']}")
            print(f"\n{COLORAMA_COLORS['YELLOW']}{COLORAMA_COLORS['BOLD']}📝 IMPORTANT NOTE:{COLORAMA_COLORS['RESET']}")
            print(f"{COLORAMA_COLORS['YELLOW']}All in silico scores must be entered manually - no automatic retrieval.{COLORAMA_COLORS['RESET']}")
            print(f"{COLORAMA_COLORS['CYAN']}{'='*80}{COLORAMA_COLORS['RESET']}")


        # Interactive mode selection (only if not test or batch mode)
        if not self.test_mode and not getattr(self, 'batch_mode', False):
            print(f"\n{COLORAMA_COLORS['CYAN']}To continue, press Enter for default mode or type a command to switch mode.{COLORAMA_COLORS['RESET']}")
            while True:
                valid_modes = {
                    'default': None,
                    '--acmg-2023': 'acmg2023',
                    '--batch': 'batch'
                }
                user_input = input(f"{COLORAMA_COLORS['YELLOW']}Mode selection (Enter for default): {COLORAMA_COLORS['RESET']}" ).strip()
                mode_message = None
                selected_mode = 'default'
                batch_file = None
                if user_input:
                    if user_input.startswith('--batch'):
                        import re
                        match = re.match(r'--batch\s+(\S+)', user_input)
                        if match:
                            candidate_file = match.group(1)
                            import os
                            # Dosya var mı kontrol et
                            if not os.path.isfile(candidate_file):
                                candidate_file_path = os.path.join(os.getcwd(), candidate_file)
                                if not os.path.isfile(candidate_file_path):
                                    print(f"{COLORAMA_COLORS['RED']}Batch file not found: {candidate_file}{COLORAMA_COLORS['RESET']}")
                                    print(f"{COLORAMA_COLORS['YELLOW']}Please provide the file name (e.g. batch.csv) if the file is in the current directory, or the full path (e.g. C:/path/to/batch.csv) if it is elsewhere.{COLORAMA_COLORS['RESET']}")
                                    continue  # Dosya yoksa tekrar mod seçimine dön
                                else:
                                    batch_file = candidate_file_path
                            else:
                                batch_file = candidate_file
                            selected_mode = 'batch'
                        else:
                            print(f"{COLORAMA_COLORS['RED']}Invalid batch mode command. Usage: --batch file.csv{COLORAMA_COLORS['RESET']}")
                            continue  # Return to mode selection
                    elif user_input == '--acmg-2023':
                        selected_mode = 'acmg2023'
                    else:
                        print(f"{COLORAMA_COLORS['RED']}Invalid mode selection. Please try again.{COLORAMA_COLORS['RESET']}")
                        continue  # Return to mode selection

                # Apply selected mode
                if selected_mode == 'acmg2023':
                    self.use_2023_guidelines = True
                    mode_message = f"{COLORAMA_COLORS['CYAN']}{COLORAMA_COLORS['BOLD']}You are now in ACMG/AMP 2023 mode.{COLORAMA_COLORS['RESET']}"
                elif selected_mode == 'batch' and batch_file:
                    self.batch_mode = True
                    self.batch_file = batch_file
                    mode_message = f"{COLORAMA_COLORS['YELLOW']}{COLORAMA_COLORS['BOLD']}You are now in Batch Mode. File: {self.batch_file}{COLORAMA_COLORS['RESET']}"
                else:
                    mode_message = f"{COLORAMA_COLORS['CYAN']}{COLORAMA_COLORS['BOLD']}You are now in ACMG/AMP 2015 (default) mode.{COLORAMA_COLORS['RESET']}"
                print(f"{COLORAMA_COLORS['GREEN']}Mode switched according to your selection.{COLORAMA_COLORS['RESET']}")
                if mode_message:
                    print(f"\n{mode_message}")
                break  # Exit loop after valid selection

        # Show mode message for direct mode selection (non-interactive)
        elif self.test_mode:
            print(f"\n{COLORAMA_COLORS['MAGENTA']}{COLORAMA_COLORS['BOLD']}You are now in Test Mode.{COLORAMA_COLORS['RESET']}")
        elif getattr(self, 'batch_mode', False):
            print(f"\n{COLORAMA_COLORS['YELLOW']}{COLORAMA_COLORS['BOLD']}You are now in Batch Mode. File: {getattr(self, 'batch_file', 'N/A')}{COLORAMA_COLORS['RESET']}")
        elif self.use_2023_guidelines:
            print(f"\n{COLORAMA_COLORS['CYAN']}{COLORAMA_COLORS['BOLD']}You are now in ACMG/AMP 2023 mode.{COLORAMA_COLORS['RESET']}")
        else:
            print(f"\n{COLORAMA_COLORS['CYAN']}{COLORAMA_COLORS['BOLD']}You are now in ACMG/AMP 2015 (default) mode.{COLORAMA_COLORS['RESET']}")

    
    def run(self):
        """Main execution method."""
        while True:
            try:
                # Collect variant data
                variant_data = self._collect_variant_data()
                
                # Evaluate evidence
                print(f"\n{COLORAMA_COLORS['CYAN']}{COLORAMA_COLORS['BOLD']}🔍 EVIDENCE EVALUATION{COLORAMA_COLORS['RESET']}")
                print(f"{COLORAMA_COLORS['BLUE']}{'-'*40}{COLORAMA_COLORS['RESET']}")
                evidence_results = self._evaluate_evidence(variant_data)
                
                # Classify variant
                print(f"\n{COLORAMA_COLORS['MAGENTA']}{COLORAMA_COLORS['BOLD']}⚖️  VARIANT CLASSIFICATION{COLORAMA_COLORS['RESET']}")
                print(f"{COLORAMA_COLORS['BLUE']}{'-'*40}{COLORAMA_COLORS['RESET']}")
                classification_result = self._classify_variant(evidence_results)
                
                # Generate report
                print(f"\n{SECTION_HEADERS['report_generation']}")
                print(f"{COLORAMA_COLORS['BLUE']}{'-'*40}{COLORAMA_COLORS['RESET']}")
                report_path = self._generate_report(variant_data, evidence_results, classification_result)
                
                # Display results
                self._display_results(classification_result, report_path, variant_data)
                
                # Ask if user wants to analyze another variant
                if not self.test_mode:
                    self._ask_for_restart()
                else:
                    break
                    
            except KeyboardInterrupt:
                print(f"\n\n{COLORAMA_COLORS['RED']}❌ Operation cancelled by user.{COLORAMA_COLORS['RESET']}")
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
            color_code = CLASSIFICATION_COLORS[classification]
            print(f"\n{color_code}📊 FINAL CLASSIFICATION: {classification}{COLORAMA_COLORS['RESET']}")
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
  python acmg_assistant.py --batch variants.csv # Batch mode (CSV)

Batch dosya formatı örneği:
chromosome,position,ref_allele,alt_allele,variant_type,consequence,gene
17,43093464,C,T,missense,missense_variant,TP53
13,32906641,G,A,frameshift,frameshift_variant,BRCA2
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
        '--batch',
        metavar='BATCH_FILE',
        type=str,
        help='Batch mode: analyze multiple variants from a CSV file'
    )

    parser.add_argument(
        '--version',
        action='version',
        version=f'ACMG Assistant v{VERSION_INFO["version"]}'
    )

    args = parser.parse_args()

    # Prevent test mode in exe
    test_mode_arg = args.test
    if getattr(sys, 'frozen', False) and test_mode_arg:
        print(f"{COLORAMA_COLORS['RED']}Test Mode is not available in exe. Please use terminal to access Test Mode.{COLORAMA_COLORS['RESET']}")
        test_mode_arg = False

    assistant = ACMGAssistant(
        use_2023_guidelines=args.acmg_2023,
        test_mode=test_mode_arg
    )

    if args.batch:
        import csv
        import os
        batch_file = args.batch
        # Dosya yolu kontrolü: önce çalışma dizininde, sonra tam yol olarak dene
        if not os.path.isfile(batch_file):
            # Çalışma dizininde yoksa, tam yol olarak dene
            batch_file_path = os.path.join(os.getcwd(), batch_file)
            if os.path.isfile(batch_file_path):
                batch_file = batch_file_path
            else:
                print(f"{COLORAMA_COLORS['RED']}Batch file not found: {args.batch}{COLORAMA_COLORS['RESET']}")
                sys.exit(1)
        print(f"\n{COLORAMA_COLORS['CYAN']}Batch modunda analiz başlatılıyor: {batch_file}{COLORAMA_COLORS['RESET']}")
        try:
            with open(batch_file, newline='', encoding='utf-8') as csvfile:
                reader = csv.DictReader(csvfile)
                results = []
                for i, row in enumerate(reader, 1):
                    print(f"\n{COLORAMA_COLORS['YELLOW']}[{i}] Varyant analiz ediliyor: {row.get('gene','?')} {row.get('chromosome','?')}:{row.get('position','?')}{COLORAMA_COLORS['RESET']}")
                    # CSV'den gelen verilerle otomatik analiz
                    basic_info = {
                        'chromosome': row.get('chromosome',''),
                        'position': row.get('position',''),
                        'ref_allele': row.get('ref_allele',''),
                        'alt_allele': row.get('alt_allele',''),
                        'variant_type': row.get('variant_type',''),
                        'consequence': row.get('consequence',''),
                        'gene': row.get('gene','')
                    }
                    # Ek alanları otomatik olarak ilgili dict'lere aktar
                    population_data = {}
                    insilico_data = {}
                    genetic_data = {}
                    functional_data = {}
                    patient_phenotypes = None
                    clinvar_data = None

                    # Alan isimleri ve eşleşen dict'ler
                    population_keys = ['gnomad_af', 'exac_af', 'topmed_af']
                    insilico_keys = ['cadd_phred', 'revel', 'sift', 'polyphen2', 'mutation_taster', 'dann', 'fathmm', 'spliceai_ag_score', 'spliceai_al_score', 'spliceai_dg_score', 'spliceai_dl_score']
                    genetic_keys = ['segregation', 'de_novo', 'denovo', 'maternity_confirmed', 'paternity_confirmed', 'inheritance']
                    functional_keys = ['case_control', 'functional_study', 'de_novo', 'denovo']
                    clinvar_keys = ['clinical_significance']

                    for key in row:
                        value = row[key]
                        if key in population_keys:
                            try:
                                population_data[key] = float(value) if value not in ('', None) else None
                            except Exception:
                                population_data[key] = value
                        elif key in insilico_keys:
                            try:
                                insilico_data[key] = float(value) if value not in ('', None) else None
                            except Exception:
                                insilico_data[key] = value
                        elif key in genetic_keys:
                            genetic_data[key] = value
                        elif key in functional_keys:
                            functional_data[key] = value
                        elif key in clinvar_keys:
                            if clinvar_data is None:
                                clinvar_data = {}
                            clinvar_data[key] = value
                        elif key == 'patient_phenotypes':
                            patient_phenotypes = value

                    variant_data = VariantData(
                        basic_info=basic_info,
                        population_data=population_data,
                        insilico_data=insilico_data,
                        genetic_data=genetic_data,
                        functional_data=functional_data,
                        patient_phenotypes=patient_phenotypes,
                        clinvar_data=clinvar_data
                    )
                    evidence_results = assistant._evaluate_evidence(variant_data)
                    classification_result = assistant._classify_variant(evidence_results)
                    report_path = assistant._generate_report(variant_data, evidence_results, classification_result)
                    assistant._display_results(classification_result, report_path, variant_data)
                    results.append({
                        'gene': basic_info['gene'],
                        'chromosome': basic_info['chromosome'],
                        'position': basic_info['position'],
                        'classification': classification_result['classification'],
                        'confidence': classification_result.get('confidence',''),
                        'report': report_path
                    })
                print(f"\n{COLORAMA_COLORS['GREEN']}Batch analiz tamamlandı. Toplam {len(results)} varyant işlendi.{COLORAMA_COLORS['RESET']}")
        except Exception as e:
            print(f"{COLORAMA_COLORS['RED']}Batch dosyası okunamadı veya hata oluştu: {str(e)}{COLORAMA_COLORS['RESET']}")
    else:
        assistant.run()

if __name__ == "__main__":
    main()

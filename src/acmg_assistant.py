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

# Fix working directory when running by double-click
# This ensures imports work correctly regardless of how the script is launched
_script_dir = os.path.dirname(os.path.abspath(__file__))
if _script_dir not in sys.path:
    sys.path.insert(0, _script_dir)
os.chdir(_script_dir)

import argparse
from typing import Dict, Any
from datetime import datetime


from core.variant_data import VariantData
from core.acmg_classifier import ACMGClassifier
from core.evidence_evaluator import EvidenceEvaluator
from utils.input_handler import InputHandler
from utils.api_client import APIClient
from utils.validators import Validators
from utils.report_generator import ReportGenerator
from config.constants import (
    VERSION_INFO,
    COLORAMA_COLORS,
    SECTION_HEADERS,
    CLASSIFICATION_COLORS,
    TEST_SCENARIOS,
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
        self.input_handler = InputHandler(
            test_mode=test_mode, api_client=self.api_client
        )
        self.validators = Validators()
        self.evidence_evaluator = EvidenceEvaluator(
            use_2023_guidelines, test_mode=test_mode
        )
        self.classifier = ACMGClassifier(use_2023_guidelines)

        # Get the correct output directory for reports
        if getattr(sys, "frozen", False):
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
        """Print welcome message with evidence summary and data flow information."""
        print(f"\n{COLORAMA_COLORS['CYAN']}{'='*80}")
        print(
            f"{COLORAMA_COLORS['CYAN']}{COLORAMA_COLORS['BOLD']}🧬 ACMG Variant Classification Assistant v{VERSION_INFO['version']}"
        )
        print(f"{COLORAMA_COLORS['CYAN']}{'='*80}")
        print(
            f"{COLORAMA_COLORS['WHITE']}An educational ACMG/AMP variant interpretation helper with automatic"
        )
        print(
            f"{COLORAMA_COLORS['WHITE']}data retrieval and interactive evidence collection.{COLORAMA_COLORS['RESET']}"
        )
        print(f"{COLORAMA_COLORS['CYAN']}{'-'*80}{COLORAMA_COLORS['RESET']}")
        print(
            f"{COLORAMA_COLORS['YELLOW']}Active Guidelines: {COLORAMA_COLORS['WHITE']}{COLORAMA_COLORS['BOLD']}ACMG/AMP {'2023' if self.use_2023_guidelines else '2015'}{COLORAMA_COLORS['RESET']}"
        )
        print(
            f"{COLORAMA_COLORS['YELLOW']}Mode: {COLORAMA_COLORS['WHITE']}{COLORAMA_COLORS['BOLD']}{'Test Mode' if self.test_mode else 'Interactive Analysis'}{COLORAMA_COLORS['RESET']}"
        )
        print(
            f"{COLORAMA_COLORS['YELLOW']}Date: {COLORAMA_COLORS['WHITE']}{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}{COLORAMA_COLORS['RESET']}"
        )
        print(f"{COLORAMA_COLORS['CYAN']}{'='*80}")
        print(
            f"{COLORAMA_COLORS['GREEN']}Created by: {VERSION_INFO['author']}{COLORAMA_COLORS['RESET']}"
        )

        # Show disclaimer and evidence summary only on initial welcome page
        if show_disclaimer:
            # Data flow note (brief)
            print(f"\n{COLORAMA_COLORS['CYAN']}{COLORAMA_COLORS['BOLD']}📋 DATA FLOW:{COLORAMA_COLORS['RESET']}")
            print(f"{COLORAMA_COLORS['WHITE']}  Predictor and population data are fetched from external APIs when available")
            print(f"  and cached locally. Interactive evidence collection is optional.{COLORAMA_COLORS['RESET']}")
            
            # Disclaimer
            print(
                f"\n{COLORAMA_COLORS['RED']}{COLORAMA_COLORS['BOLD']}⚠️  DISCLAIMER:{COLORAMA_COLORS['RESET']}"
            )
            print(
                f"{COLORAMA_COLORS['RED']}This tool is for educational and research support only and must not be{COLORAMA_COLORS['RESET']}"
            )
            print(
                f"{COLORAMA_COLORS['RED']}used as a standalone clinical decision system. Clinical variant{COLORAMA_COLORS['RESET']}"
            )
            print(
                f"{COLORAMA_COLORS['RED']}interpretations must be validated by certified professionals.{COLORAMA_COLORS['RESET']}"
            )
            print(f"{COLORAMA_COLORS['CYAN']}{'='*80}{COLORAMA_COLORS['RESET']}")

        # Interactive guideline selection (only if not test mode)
        if not self.test_mode:
            print(
                f"\n{COLORAMA_COLORS['CYAN']}Press Enter to start with default guidelines (ACMG/AMP 2015){COLORAMA_COLORS['RESET']}"
            )
            print(
                f"{COLORAMA_COLORS['CYAN']}or type '--acmg-2023' to use ACMG/AMP 2023 guidelines.{COLORAMA_COLORS['RESET']}"
            )

            while True:
                user_input = input(
                    f"{COLORAMA_COLORS['YELLOW']}Guideline selection (Enter for default): {COLORAMA_COLORS['RESET']}"
                ).strip()

                if not user_input:  # Default: ACMG 2015
                    mode_message = f"{COLORAMA_COLORS['CYAN']}{COLORAMA_COLORS['BOLD']}Using ACMG/AMP 2015 guidelines (default).{COLORAMA_COLORS['RESET']}"
                    print(mode_message)
                    break
                elif user_input == "--acmg-2023":
                    self.use_2023_guidelines = True
                    # Reinitialize with 2023 guidelines
                    self.evidence_evaluator = EvidenceEvaluator(
                        self.use_2023_guidelines, test_mode=self.test_mode
                    )
                    self.classifier = ACMGClassifier(self.use_2023_guidelines)
                    mode_message = f"{COLORAMA_COLORS['CYAN']}{COLORAMA_COLORS['BOLD']}Using ACMG/AMP 2023 guidelines.{COLORAMA_COLORS['RESET']}"
                    print(mode_message)
                    break
                else:
                    print(
                        f"{COLORAMA_COLORS['RED']}Invalid input. Please press Enter or type '--acmg-2023'.{COLORAMA_COLORS['RESET']}"
                    )
                    continue

        # Test mode message
        elif self.test_mode:
            print(
                f"\n{COLORAMA_COLORS['MAGENTA']}{COLORAMA_COLORS['BOLD']}Running in Test Mode with {'ACMG/AMP 2023' if self.use_2023_guidelines else 'ACMG/AMP 2015'} guidelines.{COLORAMA_COLORS['RESET']}"
            )
        else:
            print(
                f"\n{COLORAMA_COLORS['CYAN']}{COLORAMA_COLORS['BOLD']}You are now in ACMG/AMP 2015 (default) mode.{COLORAMA_COLORS['RESET']}"
            )

    def run(self):
        """Main execution method."""
        while True:
            try:
                # Collect variant data
                variant_data = self._collect_variant_data()

                # Evaluate evidence
                print(
                    f"\n{COLORAMA_COLORS['CYAN']}{COLORAMA_COLORS['BOLD']}🔍 EVIDENCE EVALUATION{COLORAMA_COLORS['RESET']}"
                )
                print(f"{COLORAMA_COLORS['BLUE']}{'-'*40}{COLORAMA_COLORS['RESET']}")
                evidence_results = self._evaluate_evidence(variant_data)

                # Classify variant
                print(
                    f"\n{COLORAMA_COLORS['MAGENTA']}{COLORAMA_COLORS['BOLD']}⚖️  VARIANT CLASSIFICATION{COLORAMA_COLORS['RESET']}"
                )
                print(f"{COLORAMA_COLORS['BLUE']}{'-'*40}{COLORAMA_COLORS['RESET']}")
                classification_result = self._classify_variant(evidence_results)

                # Generate report
                print(f"\n{SECTION_HEADERS['report_generation']}")
                print(f"{COLORAMA_COLORS['BLUE']}{'-'*40}{COLORAMA_COLORS['RESET']}")
                report_path = self._generate_report(
                    variant_data, evidence_results, classification_result
                )

                # Display results
                self._display_results(classification_result, report_path, variant_data)

                # Ask if user wants to analyze another variant
                if not self.test_mode:
                    self._ask_for_restart()
                else:
                    break

            except KeyboardInterrupt:
                print(
                    f"\n\n{COLORAMA_COLORS['RED']}❌ Operation cancelled by user.{COLORAMA_COLORS['RESET']}"
                )
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
        print(
            f"{COLORAMA_COLORS['GREEN']}{COLORAMA_COLORS['BOLD']}🔄 ANALYSIS COMPLETE{COLORAMA_COLORS['RESET']}"
        )
        print(f"{COLORAMA_COLORS['CYAN']}{'='*80}{COLORAMA_COLORS['RESET']}")

        while True:
            try:
                choice = (
                    input(
                        f"\n{COLORAMA_COLORS['YELLOW']}Would you like to analyze another variant? (y/n) or press Enter to exit: {COLORAMA_COLORS['RESET']}"
                    )
                    .strip()
                    .lower()
                )

                if choice == "" or choice == "n" or choice == "no":
                    print(
                        f"\n{COLORAMA_COLORS['GREEN']}👋 Thank you for using ACMG Variant Classification Assistant!{COLORAMA_COLORS['RESET']}"
                    )
                    self._pause_before_exit()
                    sys.exit(0)
                elif choice == "y" or choice == "yes":
                    print(
                        f"\n{COLORAMA_COLORS['GREEN']}🔄 Starting new analysis...{COLORAMA_COLORS['RESET']}"
                    )
                    return
                else:
                    print(
                        f"{COLORAMA_COLORS['RED']}❌ Please enter 'y' for yes, 'n' for no, or press Enter to exit.{COLORAMA_COLORS['RESET']}"
                    )

            except KeyboardInterrupt:
                print(
                    f"\n\n{COLORAMA_COLORS['RED']}❌ Operation cancelled by user.{COLORAMA_COLORS['RESET']}"
                )
                self._pause_before_exit()
                sys.exit(1)

    def _pause_before_exit(self):
        """Pause before exiting to allow user to see any error messages."""
        try:
            print(
                f"{COLORAMA_COLORS['CYAN']}Press Enter to exit...{COLORAMA_COLORS['RESET']}"
            )
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
        insilico_data = self.input_handler.collect_insilico_data(
            basic_info["variant_type"]
        )

        # Genetic data
        genetic_data = self.input_handler.collect_genetic_data()

        # Functional data
        functional_data = self.input_handler.collect_functional_data()

        # Collect patient phenotypes for PP4/BP5 (interactive mode only)
        # In test mode, phenotypes come from test scenarios
        patient_phenotypes = None
        clinvar_data = {}
        
        if self.test_mode and self.input_handler.selected_scenario:
            scenario = TEST_SCENARIOS[self.input_handler.selected_scenario]
            patient_phenotypes = scenario.get("patient_phenotypes", None)
            clinvar_data = scenario.get("clinvar_data", {})
        else:
            # Interactive mode: collect phenotypes from user for PhenotypeMatcher
            patient_phenotypes = self.input_handler.collect_patient_phenotypes()

        # Create variant data object
        variant_data = VariantData(
            basic_info=basic_info,
            population_data=population_data,
            insilico_data=insilico_data,
            genetic_data=genetic_data,
            functional_data=functional_data,
            patient_phenotypes=patient_phenotypes,
            clinvar_data=clinvar_data,
        )

        # Enrich with API data (only in non-test mode)
        if not self.test_mode:
            self._enrich_with_api_data(variant_data)

        return variant_data

    def _enrich_with_api_data(self, variant_data: VariantData):
        """Enrich variant data with API information."""
        print("\n🌐 Fetching additional data from APIs...")

        # Get ClinVar status
        print("- Checking ClinVar status...")
        clinvar_data = self.api_client.get_clinvar_status(
            variant_data.basic_info["chromosome"],
            variant_data.basic_info["position"],
            variant_data.basic_info["ref_allele"],
            variant_data.basic_info["alt_allele"],
        )
        variant_data.clinvar_data = clinvar_data
        
        # Enhanced auto-enrichment with new APIs
        enriched_data = self.api_client.auto_enrich_variant_data(variant_data)
        
        # Integrate enriched data into variant_data
        if enriched_data.get('alphamissense') and not enriched_data['alphamissense'].get('error'):
            # Add AlphaMissense score to insilico_data if available
            am_score = enriched_data['alphamissense'].get('alphamissense_score')
            if am_score is not None:
                variant_data.insilico_data['alphamissense_score'] = am_score
        
        if enriched_data.get('conservation') and not enriched_data['conservation'].get('error'):
            # Add conservation scores to insilico_data
            cons_scores = enriched_data['conservation'].get('conservation_scores', {})
            for score_type, score_value in cons_scores.items():
                variant_data.insilico_data[f'{score_type}_score'] = score_value
        
        if enriched_data.get('domains') and not enriched_data['domains'].get('error'):
            # Store domain information for PM1 evaluation
            variant_data.domain_data = enriched_data['domains']

    def _evaluate_evidence(self, variant_data: VariantData) -> Dict[str, Any]:
        """Evaluate ACMG evidence criteria with optional interactive literature-based evidence."""
        # Print pre-evaluation status summary
        self._print_evaluation_status()
        
        # First, evaluate automated criteria
        automated_results = self.evidence_evaluator.evaluate_all_criteria(variant_data)
        
        # In interactive mode, offer to collect literature-based evidence
        if not self.test_mode:
            manual_evidence = self._offer_interactive_evidence_collection(variant_data)
            if manual_evidence and manual_evidence.has_evidence():
                # Merge manual evidence with automated results
                return self.evidence_evaluator.merge_manual_with_automated(automated_results)
        
        return automated_results
    
    def _print_evaluation_status(self):
        """Print a brief status summary before evaluation begins."""
        print(f"\n{COLORAMA_COLORS['CYAN']}{'-'*60}{COLORAMA_COLORS['RESET']}")
        print(f"{COLORAMA_COLORS['CYAN']}{COLORAMA_COLORS['BOLD']}📋 PRE-EVALUATION STATUS{COLORAMA_COLORS['RESET']}")
        print(f"{COLORAMA_COLORS['CYAN']}{'-'*60}{COLORAMA_COLORS['RESET']}")
        print(f"{COLORAMA_COLORS['WHITE']}  Predictor data  : multi-source APIs + validated cache (if available){COLORAMA_COLORS['RESET']}")
        print(f"{COLORAMA_COLORS['WHITE']}  Population data : gnomAD (and configured sources) + validated cache{COLORAMA_COLORS['RESET']}")
        interactive_status = "disabled" if self.test_mode else "enabled"
        print(f"{COLORAMA_COLORS['WHITE']}  Interactive mode: {interactive_status}{COLORAMA_COLORS['RESET']}")
        print(f"{COLORAMA_COLORS['CYAN']}{'-'*60}{COLORAMA_COLORS['RESET']}")
    
    def _offer_interactive_evidence_collection(self, variant_data: VariantData):
        """
        Offer the user the option to provide literature-based evidence.
        
        This covers ACMG criteria that require human interpretation:
        - PS3/BS3: Functional studies
        - PS4: Case-control data
        - PP1/BS4: Segregation
        - PS1/PM5: Same codon/amino acid
        - PP5/BP6: Reputable source
        """
        from core.interactive_evidence import ManualEvidence
        
        print(f"\n{COLORAMA_COLORS['CYAN']}{'='*60}{COLORAMA_COLORS['RESET']}")
        print(f"{COLORAMA_COLORS['CYAN']}{COLORAMA_COLORS['BOLD']}📚 INTERACTIVE EVIDENCE COLLECTION (OPTIONAL){COLORAMA_COLORS['RESET']}")
        print(f"{COLORAMA_COLORS['CYAN']}{'='*60}{COLORAMA_COLORS['RESET']}")
        
        # Introductory explanation
        print(f"\n{COLORAMA_COLORS['WHITE']}Now starting interactive evidence collection (optional).{COLORAMA_COLORS['RESET']}")
        print(f"{COLORAMA_COLORS['WHITE']}These questions cover literature-based and judgment-heavy criteria{COLORAMA_COLORS['RESET']}")
        print(f"{COLORAMA_COLORS['WHITE']}(PS3/BS3, PS4, PP1/BS4, PS1/PM5, PP5/BP6).{COLORAMA_COLORS['RESET']}")
        print(f"{COLORAMA_COLORS['YELLOW']}The tool does NOT infer these automatically; your answers determine{COLORAMA_COLORS['RESET']}")
        print(f"{COLORAMA_COLORS['YELLOW']}which of these criteria are applied.{COLORAMA_COLORS['RESET']}")
        
        print(f"\n{COLORAMA_COLORS['CYAN']}Criteria covered:{COLORAMA_COLORS['RESET']}")
        print(f"  • PS3/BS3: Functional studies (in vitro/in vivo)")
        print(f"  • PS4: Case-control / enrichment data")
        print(f"  • PP1/BS4: Family segregation data")
        print(f"  • PS1/PM5: Known pathogenic variants at same position")
        print(f"  • PP5/BP6: Reputable source classifications")
        print(f"{COLORAMA_COLORS['YELLOW']}           (Note: PP5/BP6 are deprecated in recent ClinGen guidance){COLORAMA_COLORS['RESET']}")
        
        print(f"\n{COLORAMA_COLORS['WHITE']}💡 If you have reviewed literature for this variant,{COLORAMA_COLORS['RESET']}")
        print(f"{COLORAMA_COLORS['WHITE']}   you can provide that evidence now. Otherwise, skip.{COLORAMA_COLORS['RESET']}")
        
        while True:
            choice = input(f"\n{COLORAMA_COLORS['GREEN']}Collect interactive evidence? (y/n) [default: n]: {COLORAMA_COLORS['RESET']}").strip().lower()
            
            if choice == '' or choice == 'n' or choice == 'no':
                print(f"{COLORAMA_COLORS['CYAN']}Skipping interactive evidence collection.{COLORAMA_COLORS['RESET']}")
                return None
            elif choice == 'y' or choice == 'yes':
                # Ask which criteria to evaluate
                return self._collect_interactive_evidence(variant_data)
            else:
                print(f"{COLORAMA_COLORS['RED']}Please enter 'y' or 'n'.{COLORAMA_COLORS['RESET']}")
    
    def _collect_interactive_evidence(self, variant_data: VariantData):
        """Collect interactive evidence for selected criteria."""
        from core.interactive_evidence import ManualEvidence
        
        print(f"\n{COLORAMA_COLORS['CYAN']}Select which criteria to evaluate:{COLORAMA_COLORS['RESET']}")
        print(f"  1. All criteria (full guided interview)")
        print(f"  2. PS3/BS3 only (functional studies)")
        print(f"  3. PS4 only (case-control data)")
        print(f"  4. PP1/BS4 only (segregation)")
        print(f"  5. PS1/PM5 only (same position variants)")
        print(f"  6. PP5/BP6 only (reputable source) {COLORAMA_COLORS['YELLOW']}⚠ deprecated{COLORAMA_COLORS['RESET']}")
        print(f"  7. Custom selection")
        
        while True:
            choice = input(f"\n{COLORAMA_COLORS['GREEN']}Enter choice [1-7, default: 1]: {COLORAMA_COLORS['RESET']}").strip()
            
            if choice == '' or choice == '1':
                # Collect all
                return self.evidence_evaluator.collect_manual_evidence(variant_data)
            elif choice == '2':
                return self.evidence_evaluator.collect_manual_evidence(variant_data, criteria=['PS3_BS3'])
            elif choice == '3':
                return self.evidence_evaluator.collect_manual_evidence(variant_data, criteria=['PS4'])
            elif choice == '4':
                return self.evidence_evaluator.collect_manual_evidence(variant_data, criteria=['PP1_BS4'])
            elif choice == '5':
                return self.evidence_evaluator.collect_manual_evidence(variant_data, criteria=['PS1_PM5'])
            elif choice == '6':
                return self.evidence_evaluator.collect_manual_evidence(variant_data, criteria=['PP5_BP6'])
            elif choice == '7':
                # Custom selection
                return self._collect_custom_criteria(variant_data)
            else:
                print(f"{COLORAMA_COLORS['RED']}Please enter a number from 1-7.{COLORAMA_COLORS['RESET']}")
    
    def _collect_custom_criteria(self, variant_data: VariantData):
        """Allow user to select specific criteria to evaluate."""
        print(f"\n{COLORAMA_COLORS['CYAN']}Select criteria (comma-separated):{COLORAMA_COLORS['RESET']}")
        print(f"  Available: PS3_BS3, PS4, PP1_BS4, PS1_PM5, PP5_BP6")
        
        while True:
            selection = input(f"\n{COLORAMA_COLORS['GREEN']}Criteria to evaluate: {COLORAMA_COLORS['RESET']}").strip().upper()
            
            if not selection:
                print(f"{COLORAMA_COLORS['RED']}Please enter at least one criterion.{COLORAMA_COLORS['RESET']}")
                continue
            
            # Parse selection
            criteria = [c.strip() for c in selection.split(',')]
            valid_criteria = ['PS3_BS3', 'PS4', 'PP1_BS4', 'PS1_PM5', 'PP5_BP6']
            
            invalid = [c for c in criteria if c not in valid_criteria]
            if invalid:
                print(f"{COLORAMA_COLORS['RED']}Invalid criteria: {', '.join(invalid)}{COLORAMA_COLORS['RESET']}")
                print(f"Valid options: {', '.join(valid_criteria)}")
                continue
            
            return self.evidence_evaluator.collect_manual_evidence(variant_data, criteria=criteria)

    def _classify_variant(self, evidence_results: Dict[str, Any]) -> Dict[str, Any]:
        """Classify variant based on evidence."""
        return self.classifier.classify(evidence_results)

    def _generate_report(
        self,
        variant_data: VariantData,
        evidence_results: Dict[str, Any],
        classification_result: Dict[str, Any],
    ) -> str:
        """Generate detailed classification report."""
        return self.report_generator.generate_report(
            variant_data, evidence_results, classification_result
        )

    def _display_results(
        self, classification_result: Dict[str, Any], report_path: str, variant_data=None
    ):
        """Display classification results to user."""
        print("\n" + "=" * 80)
        print("🎯 CLASSIFICATION RESULTS")
        print("=" * 80)

        # Main classification
        classification = classification_result["classification"]
        confidence = classification_result.get("confidence", "N/A")

        # Display classification result with color
        classification = classification_result["classification"]
        confidence = classification_result["confidence"]

        # Get color for classification
        if classification in CLASSIFICATION_COLORS:
            color_code = CLASSIFICATION_COLORS[classification]
            print(
                f"\n{color_code}📊 FINAL CLASSIFICATION: {classification}{COLORAMA_COLORS['RESET']}"
            )
        else:
            print(
                f"\n{COLORAMA_COLORS['MAGENTA']}📊 FINAL CLASSIFICATION: {classification}{COLORAMA_COLORS['RESET']}"
            )

        print(
            f"{COLORAMA_COLORS['CYAN']}🔍 CONFIDENCE LEVEL: {confidence}{COLORAMA_COLORS['RESET']}"
        )

        # Applied criteria
        print(
            f"\n{COLORAMA_COLORS['BLUE']}{COLORAMA_COLORS['BOLD']}📋 APPLIED CRITERIA:{COLORAMA_COLORS['RESET']}"
        )
        for category, criteria in classification_result["applied_criteria"].items():
            if criteria:
                print(
                    f"  {COLORAMA_COLORS['WHITE']}{category}: {', '.join(criteria)}{COLORAMA_COLORS['RESET']}"
                )

        # Suggestions
        if classification_result.get("suggestions"):
            print(
                f"\n{COLORAMA_COLORS['YELLOW']}{COLORAMA_COLORS['BOLD']}💡 SUGGESTIONS:{COLORAMA_COLORS['RESET']}"
            )
            for suggestion in classification_result["suggestions"]:
                print(
                    f"  {COLORAMA_COLORS['YELLOW']}• {suggestion}{COLORAMA_COLORS['RESET']}"
                )

        # Report location
        print(
            f"\n{COLORAMA_COLORS['GREEN']}📄 DETAILED REPORT: {report_path}{COLORAMA_COLORS['RESET']}"
        )

        # Varsome URL
        if variant_data:
            varsome_url = self.api_client.generate_varsome_url(
                variant_data.basic_info["chromosome"],
                variant_data.basic_info["position"],
                variant_data.basic_info["ref_allele"],
                variant_data.basic_info["alt_allele"],
            )
            if varsome_url:
                print(
                    f"{COLORAMA_COLORS['CYAN']}🔗 VARSOME URL: {varsome_url}{COLORAMA_COLORS['RESET']}"
                )

        print(f"{COLORAMA_COLORS['CYAN']}{'='*80}{COLORAMA_COLORS['RESET']}")


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="ACMG Variant Classification Assistant - Interactive Single-Variant Analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python acmg_assistant.py                    # Interactive mode with ACMG 2015 guidelines
  python acmg_assistant.py --acmg-2023        # Interactive mode with ACMG 2023 guidelines
  python acmg_assistant.py --test             # Run in test mode with sample data
  python acmg_assistant.py --test --acmg-2023 # Test mode with 2023 guidelines

Note: This tool is designed for careful, interactive analysis of individual variants.
Each variant requires clinical judgment and literature review for accurate classification.
""",
    )

    parser.add_argument(
        "--acmg-2023",
        action="store_true",
        help="Use ACMG/AMP 2023 guidelines instead of 2015",
    )

    parser.add_argument(
        "--test", action="store_true", help="Run in test mode with sample data"
    )

    parser.add_argument(
        "--version",
        action="version",
        version=f'ACMG Assistant v{VERSION_INFO["version"]}',
    )

    args = parser.parse_args()

    # Prevent test mode in exe
    test_mode_arg = args.test
    if getattr(sys, "frozen", False) and test_mode_arg:
        print(
            f"{COLORAMA_COLORS['RED']}Test Mode is not available in exe. Please use terminal to access Test Mode.{COLORAMA_COLORS['RESET']}"
        )
        test_mode_arg = False

    assistant = ACMGAssistant(
        use_2023_guidelines=args.acmg_2023, test_mode=test_mode_arg
    )

    # Run interactive single-variant analysis
    assistant.run()


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"\n{'='*60}")
        print(f"ERROR: {e}")
        print(f"{'='*60}")
        import traceback
        traceback.print_exc()
    finally:
        # Keep console window open when run by double-clicking
        if getattr(sys, "frozen", False) or sys.stdin.isatty():
            input("\nPress Enter to exit...")

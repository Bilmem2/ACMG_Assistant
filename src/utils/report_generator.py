"""
Report Generator Module
======================

This module generates comprehensive reports for ACMG variant classification.
"""

import os
import sys
from datetime import datetime
from typing import Dict, List, Optional, Any
from colorama import Fore, Style, init
from config.constants import (
    OUTPUT_SETTINGS, VERSION_INFO, 
    CLASSIFICATION_COLORS, SECTION_HEADERS,
    COLORAMA_COLORS, SUCCESS_MESSAGES
)

# Initialize colorama
init()

class ReportGenerator:
    """
    Generates detailed classification reports in various formats.
    
    This class creates comprehensive reports including all variant data,
    evidence evaluation, and classification results.
    """
    
    def __init__(self, output_dir: str = "."):
        """
        Initialize the report generator.
        
        Args:
            output_dir (str): Directory for output files
        """
        # For executable compatibility, ensure we always use absolute path
        if output_dir == ".":
            # For executable, use the directory where the executable is located
            if getattr(sys, 'frozen', False):
                # Running as executable
                self.output_dir = os.path.dirname(sys.executable)
            else:
                # Running as script
                script_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
                self.output_dir = script_dir
        else:
            self.output_dir = os.path.abspath(output_dir)
        
        # Ensure output directory exists
        os.makedirs(self.output_dir, exist_ok=True)
        
        self.report_filename = OUTPUT_SETTINGS['report_filename']
        self.log_filename = OUTPUT_SETTINGS['log_filename']
    
    def print_colored_classification(self, classification: str):
        """Print classification result with appropriate colors."""
        if classification in CLASSIFICATION_COLORS:
            color = CLASSIFICATION_COLORS[classification]
            print(f"\n{color}🎯 CLASSIFICATION RESULT: {classification}{COLORAMA_COLORS['RESET']}")
        else:
            print(f"\n{COLORAMA_COLORS['MAGENTA']}🎯 CLASSIFICATION RESULT: {classification}{COLORAMA_COLORS['RESET']}")
    
    def generate_report(self, variant_data, evidence_results: Dict[str, Any], 
                       classification_result: Dict[str, Any]) -> str:
        """
        Generate a comprehensive classification report.
        
        Args:
            variant_data: VariantData object
            evidence_results: Evidence evaluation results
            classification_result: Classification results
            
        Returns:
            str: Path to generated report file
        """
        report_path = os.path.join(self.output_dir, self.report_filename)
        
        with open(report_path, 'w', encoding='utf-8') as f:
            self._write_header(f)
            self._write_variant_summary(f, variant_data)
            self._write_classification_summary(f, classification_result)
            self._write_evidence_details(f, evidence_results)
            self._write_population_analysis(f, variant_data)
            self._write_insilico_analysis(f, variant_data, evidence_results)
            self._write_functional_analysis(f, variant_data)
            self._write_statistical_analysis(f, evidence_results)
            self._write_recommendations(f, classification_result)
            self._write_footer(f)
        
        return report_path
    
    def _write_header(self, f):
        """Write report header."""
        f.write("="*80 + "\n")
        f.write("🧬 ACMG VARIANT CLASSIFICATION REPORT\n")
        f.write("="*80 + "\n\n")
        f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write("\n" + "="*80 + "\n\n")
    
    def _write_variant_summary(self, f, variant_data):
        """Write variant summary section."""
        f.write("📋 VARIANT SUMMARY\n")
        f.write("-"*50 + "\n\n")
        
        basic_info = variant_data.basic_info
        
        f.write(f"Gene: {basic_info.get('gene', 'N/A')}\n")
        f.write(f"Chromosome: {basic_info.get('chromosome', 'N/A')}\n")
        f.write(f"Position: {basic_info.get('position', 'N/A')}\n")
        f.write(f"Reference Allele: {basic_info.get('ref_allele', 'N/A')}\n")
        f.write(f"Alternate Allele: {basic_info.get('alt_allele', 'N/A')}\n")
        f.write(f"Variant Type: {basic_info.get('variant_type', 'N/A')}\n")
        f.write(f"Consequence: {basic_info.get('consequence', 'N/A')}\n")
        
        if basic_info.get('cdna_change'):
            f.write(f"cDNA Change: {basic_info['cdna_change']}\n")
        
        if basic_info.get('protein_change'):
            f.write(f"Protein Change: {basic_info['protein_change']}\n")
        
        # Variant ID
        variant_id = variant_data.get_variant_id()
        f.write(f"Variant ID: {variant_id}\n")
        
        # HGVS nomenclature
        hgvs = variant_data.get_hgvs_nomenclature()
        if hgvs:
            f.write(f"HGVS Nomenclature:\n")
            for hgvs_type, hgvs_string in hgvs.items():
                f.write(f"  {hgvs_type.capitalize()}: {hgvs_string}\n")
        
        f.write("\n")
    
    def _write_classification_summary(self, f, classification_result: Dict[str, Any]):
        """Write classification summary section."""
        f.write("🎯 CLASSIFICATION SUMMARY\n")
        f.write("-"*50 + "\n\n")
        
        classification = classification_result['classification']
        confidence = classification_result['confidence']
        
        f.write(f"FINAL CLASSIFICATION: {classification}\n")
        f.write(f"CONFIDENCE LEVEL: {confidence}\n")
        f.write(f"GUIDELINES VERSION: ACMG/AMP {classification_result.get('guidelines_version', '2015')}\n\n")
        
        # Applied criteria
        applied_criteria = classification_result.get('applied_criteria', {})        
        if applied_criteria.get('pathogenic'):
            f.write(f"PATHOGENIC CRITERIA: {', '.join(applied_criteria['pathogenic'])}\n")
        
        if applied_criteria.get('benign'):
            f.write(f"BENIGN CRITERIA: {', '.join(applied_criteria['benign'])}\n")
        
        # Evidence counts
        pathogenic_counts = classification_result.get('pathogenic_counts', {})
        benign_counts = classification_result.get('benign_counts', {})
        
        f.write(f"\nEVIDENCE COUNTS:\n")
        f.write(f"  Pathogenic Evidence:\n")
        f.write(f"    Very Strong: {pathogenic_counts.get('very_strong', 0)}\n")
        f.write(f"    Strong: {pathogenic_counts.get('strong', 0)}\n")
        f.write(f"    Moderate: {pathogenic_counts.get('moderate', 0)}\n")
        f.write(f"    Supporting: {pathogenic_counts.get('supporting', 0)}\n")
        
        f.write(f"  Benign Evidence:\n")
        f.write(f"    Stand-alone: {benign_counts.get('stand_alone', 0)}\n")
        f.write(f"    Strong: {benign_counts.get('strong', 0)}\n")
        f.write(f"    Supporting: {benign_counts.get('supporting', 0)}\n")
        
        # VAMPP score
        vampp_score = classification_result.get('vampp_score')
        if vampp_score is not None:
            f.write(f"\nCOMPUTATIONAL METASCORE: {vampp_score:.3f}\n")
        
        # Conflicts
        conflicts = classification_result.get('conflicts', [])
        if conflicts:
            f.write(f"\nCONFLICTS:\n")
            for conflict in conflicts:
                f.write(f"  ⚠️  {conflict}\n")
        
        f.write("\n")
    
    def _write_evidence_details(self, f, evidence_results: Dict[str, Any]):
        """Write detailed evidence evaluation."""
        f.write("🔍 EVIDENCE EVALUATION DETAILS\n")
        f.write("-"*50 + "\n\n")
        
        # Pathogenic criteria
        pathogenic_criteria = evidence_results.get('pathogenic_criteria', {})
        if pathogenic_criteria:
            f.write("PATHOGENIC CRITERIA:\n")
            for criterion, result in pathogenic_criteria.items():
                if result.get('applies', False):
                    f.write(f"  ✓ {criterion} ({result.get('strength', 'Unknown')}): {result.get('details', '')}\n")
                else:
                    f.write(f"  ✗ {criterion}: {result.get('details', 'Not applicable')}\n")
            f.write("\n")
        
        # Benign criteria
        benign_criteria = evidence_results.get('benign_criteria', {})
        if benign_criteria:
            f.write("BENIGN CRITERIA:\n")
            for criterion, result in benign_criteria.items():
                if result.get('applies', False):
                    f.write(f"  ✓ {criterion} ({result.get('strength', 'Unknown')}): {result.get('details', '')}\n")
                else:
                    f.write(f"  ✗ {criterion}: {result.get('details', 'Not applicable')}\n")
            f.write("\n")
    
    def _write_population_analysis(self, f, variant_data):
        """Write population frequency analysis."""
        f.write("🌍 POPULATION FREQUENCY ANALYSIS\n")
        f.write("-"*50 + "\n\n")
        
        population_data = variant_data.population_data
        
        if not population_data:
            f.write("No population frequency data available.\n\n")
            return
        
        f.write("FREQUENCY DATA:\n")
        
        if population_data.get('gnomad_af') is not None:
            f.write(f"  gnomAD Overall AF: {population_data['gnomad_af']:.2e}\n")
        
        if population_data.get('gnomad_af_popmax') is not None:
            f.write(f"  gnomAD PopMax AF: {population_data['gnomad_af_popmax']:.2e}\n")
        
        if population_data.get('gnomad_hom_count') is not None:
            f.write(f"  gnomAD Homozygous Count: {population_data['gnomad_hom_count']}\n")
        
        if population_data.get('gnomad_het_count') is not None:
            f.write(f"  gnomAD Heterozygous Count: {population_data['gnomad_het_count']}\n")
        
        if population_data.get('exac_af') is not None:
            f.write(f"  ExAC AF: {population_data['exac_af']:.2e}\n")
        
        if population_data.get('disease_prevalence'):
            f.write(f"  Disease Prevalence: {population_data['disease_prevalence']}\n")
        
        f.write("\n")
    
    def _write_insilico_analysis(self, f, variant_data, evidence_results: Dict[str, Any]):
        """Write in silico prediction analysis."""
        f.write("🤖 IN SILICO PREDICTION ANALYSIS\n")
        f.write("-"*50 + "\n\n")
        
        insilico_data = variant_data.insilico_data
        
        if not insilico_data:
            f.write("No in silico prediction data available.\n\n")
            return
        
        f.write("PREDICTION SCORES:\n")
        
        for predictor, score in insilico_data.items():
            if score is not None:
                f.write(f"  {predictor.upper()}: {score}\n")
        
        # VAMPP score details
        vampp_score = evidence_results.get('vampp_score')
        if vampp_score is not None:
            f.write(f"\nCOMPUTATIONAL METASCORE: {vampp_score:.3f}\n")
            f.write("  This score integrates multiple in silico predictors using\n")
            f.write("  a sophisticated weighted approach with variant-type specific optimization.\n")
        
        f.write("\n")
    
    def _write_functional_analysis(self, f, variant_data):
        """Write functional studies analysis."""
        f.write("🔬 FUNCTIONAL STUDIES ANALYSIS\n")
        f.write("-"*50 + "\n\n")
        
        functional_data = variant_data.functional_data
        
        if not functional_data:
            f.write("No functional studies data available.\n\n")
            return
        
        f.write("FUNCTIONAL DATA:\n")
        
        if functional_data.get('segregation'):
            f.write(f"  Segregation: {functional_data['segregation']}\n")
        
        if functional_data.get('denovo'):
            f.write(f"  De novo Status: {functional_data['denovo']}\n")
        
        if functional_data.get('functional_studies'):
            f.write(f"  Functional Studies: {functional_data['functional_studies']}\n")
        
        if functional_data.get('phenotype_match'):
            f.write(f"  Phenotype Match: {functional_data['phenotype_match']}\n")
        
        # Case-control data
        if functional_data.get('case_control') == 'yes':
            f.write(f"\nCASE-CONTROL DATA:\n")
            f.write(f"  Cases with variant: {functional_data.get('cases_with_variant', 'N/A')}\n")
            f.write(f"  Total cases: {functional_data.get('total_cases', 'N/A')}\n")
            f.write(f"  Controls with variant: {functional_data.get('controls_with_variant', 'N/A')}\n")
            f.write(f"  Total controls: {functional_data.get('total_controls', 'N/A')}\n")
        
        f.write("\n")
    
    def _write_statistical_analysis(self, f, evidence_results: Dict[str, Any]):
        """Write statistical analysis section."""
        f.write("📊 STATISTICAL ANALYSIS\n")
        f.write("-"*50 + "\n\n")
        
        statistical_tests = evidence_results.get('statistical_tests', {})
        
        if not statistical_tests:
            f.write("No statistical tests performed.\n\n")
            return
        
        # Fisher's exact test
        fisher_test = statistical_tests.get('fisher_exact')
        if fisher_test:
            f.write(f"FISHER'S EXACT TEST:\n")
            f.write(f"  P-value: {fisher_test['p_value']:.2e}\n")
            f.write(f"  Significant: {'Yes' if fisher_test['significant'] else 'No'}\n")
            f.write(f"  Interpretation: {'Supports pathogenic classification' if fisher_test['significant'] else 'No significant association'}\n")
        
        f.write("\n")
    
    def _write_recommendations(self, f, classification_result: Dict[str, Any]):
        """Write recommendations section."""
        f.write("💡 RECOMMENDATIONS\n")
        f.write("-"*50 + "\n\n")
        
        suggestions = classification_result.get('suggestions', [])
        
        if suggestions:
            f.write("SUGGESTIONS FOR ADDITIONAL EVIDENCE:\n")
            for i, suggestion in enumerate(suggestions, 1):
                f.write(f"  {i}. {suggestion}\n")
        else:
            f.write("No additional recommendations at this time.\n")
        
        f.write("\nRESOURCES FOR FURTHER INVESTIGATION:\n")
        f.write("  • ClinVar: https://www.ncbi.nlm.nih.gov/clinvar/\n")
        f.write("  • Varsome: https://varsome.com/\n")
        f.write("  • gnomAD: https://gnomad.broadinstitute.org/\n")
        f.write("  • PubMed: https://pubmed.ncbi.nlm.nih.gov/\n")
        f.write("  • Orphanet: https://www.orpha.net/\n")
        
        f.write("\n")
    
    def _write_footer(self, f):
        """Write report footer."""
        f.write("="*80 + "\n")
        f.write("⚠️  DISCLAIMER\n")
        f.write("="*80 + "\n\n")
        
        f.write("This report is generated for research and educational purposes only.\n")
        f.write("Clinical variant interpretations must be validated by certified professionals\n")
        f.write("following institutional and national guidelines. Always cross-check results\n")
        f.write("with primary sources (ClinVar, Varsome, etc.). The author is not responsible\n")
        f.write("for clinical decisions made using this tool.\n\n")
        
        f.write("For questions or feedback, please contact:\n")
        f.write(f"Author: {VERSION_INFO['author']}\n")
        f.write("Email: cansevilmiss@gmail.com\n\n")
        
        f.write("Generated by ACMG Variant Classification Assistant\n")
        f.write(f"Version: {VERSION_INFO['version']}\n")
        f.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write("="*80 + "\n")
    
    def generate_log_entry(self, variant_data, classification_result: Dict[str, Any]) -> str:
        """Generate a log entry for the classification."""
        log_path = os.path.join(self.output_dir, self.log_filename)
        
        # Create log entry
        timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        variant_id = variant_data.get_variant_id()
        classification = classification_result['classification']
        confidence = classification_result['confidence']
        
        log_entry = f"{timestamp} | {variant_id} | {classification} | {confidence}\n"
        
        # Append to log file
        with open(log_path, 'a', encoding='utf-8') as f:
            f.write(log_entry)
        
        return log_path
    
    def generate_json_report(self, variant_data, evidence_results: Dict[str, Any], 
                           classification_result: Dict[str, Any]) -> str:
        """Generate JSON format report."""
        import json
        
        json_data = {
            'metadata': {
                'timestamp': datetime.now().isoformat(),
                'version': VERSION_INFO['version'],
                'guidelines_version': classification_result.get('guidelines_version', '2015')
            },
            'variant_data': variant_data.to_dict(),
            'evidence_results': evidence_results,
            'classification_result': classification_result
        }
        
        json_filename = self.report_filename.replace('.txt', '.json')
        json_path = os.path.join(self.output_dir, json_filename)
        
        with open(json_path, 'w', encoding='utf-8') as f:
            json.dump(json_data, f, indent=2, default=str)
        
        return json_path

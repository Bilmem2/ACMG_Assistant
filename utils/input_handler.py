"""
Input Handler Module
==================

This module handles all user input collection for the ACMG Assistant.
It provides interactive prompts, validation, and test mode functionality.
"""

import sys
import re
from typing import Dict, List, Optional, Any, Union, Callable
from colorama import Fore, Style, init
from config.constants import (
    COLORAMA_COLORS, SECTION_HEADERS, ERROR_MESSAGES, SUCCESS_MESSAGES,
    WARNING_MESSAGES, INFO_MESSAGES, get_colored_message, TEST_MODE_DATA,
    VALIDATION_PATTERNS, ALIASES, ALL_VARIANT_CONSEQUENCES
)

# Initialize colorama
init()

class InputHandler:
    """
    Handles all user input collection and validation.
    
    This class provides methods to collect different types of variant data
    from user input, with validation and error handling.
    """
    
    def __init__(self, test_mode: bool = False, api_client=None):
        """
        Initialize the input handler.
        
        Args:
            test_mode (bool): Whether to run in test mode with sample data
            api_client: API client instance for fetching data
        """
        self.test_mode = test_mode
        self.api_client = api_client
        self.collected_data = {}
    
    def print_header(self, title: str):
        """Print a formatted header."""
        print(f"\n{COLORAMA_COLORS['CYAN']}{'='*60}{COLORAMA_COLORS['RESET']}")
        print(f"{COLORAMA_COLORS['CYAN']}{COLORAMA_COLORS['BOLD']}{title.center(60)}{COLORAMA_COLORS['RESET']}")
        print(f"{COLORAMA_COLORS['CYAN']}{'='*60}{COLORAMA_COLORS['RESET']}\n")
        
    def print_section(self, section_name: str):
        """Print a formatted section header."""
        if section_name in SECTION_HEADERS:
            header = SECTION_HEADERS[section_name]
            print(f"\n{COLORAMA_COLORS['BLUE']}{COLORAMA_COLORS['BOLD']}{header}{COLORAMA_COLORS['RESET']}")
        else:
            print(f"\n{COLORAMA_COLORS['BLUE']}{COLORAMA_COLORS['BOLD']}📋 {section_name}{COLORAMA_COLORS['RESET']}")
        print(f"{COLORAMA_COLORS['BLUE']}{'-'*40}{COLORAMA_COLORS['RESET']}")
        
    def print_success(self, message: str):
        """Print a success message."""
        print(f"{COLORAMA_COLORS['GREEN']}✅ {message}{COLORAMA_COLORS['RESET']}")
        
    def print_warning(self, message: str):
        """Print a warning message."""
        print(f"{COLORAMA_COLORS['YELLOW']}⚠️ {message}{COLORAMA_COLORS['RESET']}")
        
    def print_error(self, message: str):
        """Print an error message."""
        print(f"{COLORAMA_COLORS['RED']}❌ {message}{COLORAMA_COLORS['RESET']}")
        
    def print_info(self, message: str):
        """Print an info message."""
        print(f"{COLORAMA_COLORS['CYAN']}ℹ️ {message}{COLORAMA_COLORS['RESET']}")
    
    def collect_basic_info(self) -> Dict[str, Any]:
        """
        Collect basic variant information.
        
        Returns:
            Dict[str, Any]: Basic variant information
        """
        self.print_section("basic_info")
        
        if self.test_mode:
            self.print_info("Using test mode sample data...")
            return TEST_MODE_DATA['basic_info'].copy()
        
        basic_info = {}
        
        # Gene symbol
        basic_info['gene'] = self._prompt_input(
            "Gene symbol (e.g., BRCA1): ",
            validator=self._validate_gene_symbol,
            required=True
        )
        
        # Chromosome - automatically fetch from Ensembl
        chromosome = None
        if self.api_client:
            try:
                chromosome = self.api_client.get_chromosome_from_ensembl(basic_info['gene'])
            except Exception as e:
                self.print_error(f"Error during chromosome lookup: {e}")
                chromosome = None
        
        if chromosome:
            basic_info['chromosome'] = chromosome
            # Ask user if they want to override the automatic result
            override = self._prompt_input(
                f"Chromosome {chromosome} found. Press Enter to accept or type a different value: ",
                required=False
            )
            if override and override.strip():
                if self._validate_chromosome(override.strip()):
                    basic_info['chromosome'] = override.strip()
                    self.print_success(f"Using manual chromosome: {override.strip()}")
                else:
                    self.print_error("Invalid chromosome format. Using automatic result.")
        else:
            self.print_warning("Could not fetch chromosome automatically.")
            self.print_info("Valid chromosomes: 1-22, X, Y, MT (or chr1, chr2, etc.)")
            self.print_info("Leave empty to skip chromosome-based features")
            chromosome_input = self._prompt_input(
                "Chromosome (1-22, X, Y, MT) or press Enter to skip: ",
                validator=self._validate_chromosome,
                required=False
            )
            if chromosome_input and chromosome_input.strip():
                basic_info['chromosome'] = chromosome_input.strip()
                self.print_success(f"Using chromosome: {chromosome_input.strip()}")
            else:
                basic_info['chromosome'] = None
                self.print_info("Chromosome skipped - some features may be limited")
        
        # Position - Optional but recommended for API data fetching
        self.print_info("Genomic position enables automatic ClinVar status fetching")
        self.print_info("Find this on Varsome.com, UCSC Genome Browser, or Ensembl")
        self.print_info("Can be skipped if not available - will rely on manual data entry")
        position_input = self._prompt_input(
            "Genomic position (numbers only, optional): ",
            validator=self._validate_position,
            required=False,
            convert_type=int
        )
        basic_info['position'] = position_input if position_input else None
        
        # Reference allele - Optional for basic classification
        self.print_info("Reference/alternate alleles enable ClinVar lookup and population analysis")
        self.print_info("Can be skipped if focus is on functional/segregation evidence")
        ref_input = self._prompt_input(
            "Reference allele (e.g., C) or press Enter to skip: ",
            validator=self._validate_allele,
            required=False
        )
        basic_info['ref_allele'] = ref_input.upper() if ref_input else None
        
        # Alternate allele - Optional for basic classification
        alt_input = self._prompt_input(
            "Alternate allele (e.g., T) or press Enter to skip: ",
            validator=self._validate_allele,
            required=False
        )
        basic_info['alt_allele'] = alt_input.upper() if alt_input else None
        
        # Inform user about skipped features
        if not basic_info['ref_allele'] or not basic_info['alt_allele']:
            self.print_info("Allele information skipped - ClinVar lookup and some population features may be limited")
        
        # cDNA change
        basic_info['cdna_change'] = self._prompt_input(
            "cDNA change (HGVS format, e.g., c.5266dupC): ",
            validator=self._validate_hgvs_cdna,
            required=False
        )
        
        # Protein change
        basic_info['protein_change'] = self._prompt_input(
            "Protein change (HGVS format, e.g., p.Gln1756Profs*74): ",
            validator=self._validate_hgvs_protein,
            required=False
        )
        
        # Variant type
        print(f"\n{COLORAMA_COLORS['CYAN']}📊 Variant Type Selection:{COLORAMA_COLORS['RESET']}")
        print(f"{COLORAMA_COLORS['GREEN']}Common Types:{COLORAMA_COLORS['RESET']}")
        print("   • missense (miss)    - Amino acid change")
        print("   • nonsense (non)     - Premature stop codon")
        print("   • frameshift (fs)    - Reading frame shift")
        print("   • splice (spl)       - Affects splicing")
        print("   • synonymous (syn)   - Silent mutation")
        print("   • inframe_indel (indel) - In-frame insertion/deletion")
        print("   • other              - Other variant types")
        print(f"\n{COLORAMA_COLORS['BLUE']}💡 You can use short forms in parentheses{COLORAMA_COLORS['RESET']}")
        
        variant_types = ['missense', 'nonsense', 'frameshift', 'splice', 'synonymous', 'inframe_indel', 'other']
        variant_aliases = {
            'missense': 'miss', 'nonsense': 'non', 'frameshift': 'fs',
            'splice': 'spl', 'synonymous': 'syn', 'inframe_indel': 'indel'
        }
        basic_info['variant_type'] = self._prompt_choice(
            "Variant type: ",
            choices=variant_types,
            required=True,
            aliases=variant_aliases
        )
        
        # VEP consequence
        print(f"\n{COLORAMA_COLORS['CYAN']}📝 Variant Consequence Selection:{COLORAMA_COLORS['RESET']}")
        print(f"{COLORAMA_COLORS['RED']}🔴 High Impact:{COLORAMA_COLORS['RESET']}")
        print("   • frameshift_variant, stop_gained, stop_lost, start_lost")
        print("   • splice_acceptor_variant, splice_donor_variant")
        print(f"{COLORAMA_COLORS['YELLOW']}🟡 Moderate Impact:{COLORAMA_COLORS['RESET']}")
        print("   • missense_variant, inframe_deletion, inframe_insertion")
        print("   • splice_region_variant")
        print(f"{COLORAMA_COLORS['GREEN']}🟢 Low Impact:{COLORAMA_COLORS['RESET']}")
        print("   • synonymous_variant, start_retained_variant, stop_retained_variant")
        print(f"{COLORAMA_COLORS['BLUE']}🔵 Modifier:{COLORAMA_COLORS['RESET']}")
        print("   • intron_variant, 5_prime_UTR_variant, 3_prime_UTR_variant")
        print("   • regulatory_region_variant")
        print(f"{COLORAMA_COLORS['WHITE']}⚪ Other:{COLORAMA_COLORS['RESET']}")
        print("   • other, unknown")
        print(f"\n{COLORAMA_COLORS['BLUE']}💡 Type the exact consequence name or 'other' if not listed{COLORAMA_COLORS['RESET']}")
        
        basic_info['consequence'] = self._prompt_choice(
            "Select variant consequence: ",
            choices=ALL_VARIANT_CONSEQUENCES,
            required=False
        )
        
        # ClinVar status - automatically fetch if all required data available
        if (self.api_client and basic_info.get('chromosome') and basic_info.get('position') and 
            basic_info.get('ref_allele') and basic_info.get('alt_allele')):
            print(f"\n{COLORAMA_COLORS['CYAN']}🌐 Fetching ClinVar status...{COLORAMA_COLORS['RESET']}")
            try:
                clinvar_data = self.api_client.get_clinvar_status(
                    basic_info['chromosome'],
                    basic_info['position'],
                    basic_info['ref_allele'],
                    basic_info['alt_allele']
                )
                
                if clinvar_data and clinvar_data.get('significance'):
                    significance = clinvar_data['significance']
                    print(f"{COLORAMA_COLORS['GREEN']}✓ ClinVar status: {significance}{COLORAMA_COLORS['RESET']}")
                    
                    # Ask user if they want to override
                    override = self._prompt_choice(
                        f"ClinVar shows '{significance}'. Use this or enter manually? ",
                        choices=['use_clinvar', 'manual'],
                        aliases={'use_clinvar': 'use', 'manual': 'man'},
                        required=True
                    )
                    
                    if override == 'use_clinvar':
                        basic_info['clinvar_status'] = significance
                    else:
                        basic_info['clinvar_status'] = self._get_manual_clinvar_status()
                else:
                    print(f"{COLORAMA_COLORS['YELLOW']}⚠️  ClinVar status not found. Please enter manually.{COLORAMA_COLORS['RESET']}")
                    basic_info['clinvar_status'] = self._get_manual_clinvar_status()
                    
            except Exception as e:
                print(f"{COLORAMA_COLORS['YELLOW']}⚠️  Error fetching ClinVar data: {e}{COLORAMA_COLORS['RESET']}")
                basic_info['clinvar_status'] = self._get_manual_clinvar_status()
        else:
            # Missing required data for automatic lookup
            missing_fields = []
            if not basic_info.get('chromosome'):
                missing_fields.append('chromosome')
            if not basic_info.get('position'):
                missing_fields.append('position')
            if not basic_info.get('ref_allele'):
                missing_fields.append('reference allele')
            if not basic_info.get('alt_allele'):
                missing_fields.append('alternate allele')
            
            if missing_fields:
                print(f"\n{COLORAMA_COLORS['YELLOW']}ℹ️  Missing {', '.join(missing_fields)} - ClinVar lookup skipped{COLORAMA_COLORS['RESET']}")
                print(f"{COLORAMA_COLORS['BLUE']}💡 You can still proceed with manual evidence entry{COLORAMA_COLORS['RESET']}")
            
            basic_info['clinvar_status'] = self._get_manual_clinvar_status()
        
        return basic_info
    
    def _get_manual_clinvar_status(self) -> Optional[str]:
        """Get ClinVar status manually from user."""
        print(f"\n{COLORAMA_COLORS['CYAN']}📊 ClinVar Status Selection:{COLORAMA_COLORS['RESET']}")
        print(f"{COLORAMA_COLORS['RED']}• pathogenic (path)        - Definitely causes disease{COLORAMA_COLORS['RESET']}")
        print(f"{COLORAMA_COLORS['YELLOW']}• likely_pathogenic (lpath) - Probably causes disease{COLORAMA_COLORS['RESET']}")
        print(f"{COLORAMA_COLORS['BLUE']}• vus (vus)               - Variant of uncertain significance{COLORAMA_COLORS['RESET']}")
        print(f"{COLORAMA_COLORS['GREEN']}• likely_benign (lben)     - Probably doesn't cause disease{COLORAMA_COLORS['RESET']}")
        print(f"{COLORAMA_COLORS['GREEN']}• benign (ben)            - Definitely doesn't cause disease{COLORAMA_COLORS['RESET']}")
        print(f"{COLORAMA_COLORS['WHITE']}• not_in_clinvar (none)   - Not found in ClinVar database{COLORAMA_COLORS['RESET']}")
        print(f"\n{COLORAMA_COLORS['BLUE']}💡 You can use short forms in parentheses or press Enter to skip{COLORAMA_COLORS['RESET']}")
        
        clinvar_choices = ['pathogenic', 'likely_pathogenic', 'vus', 'likely_benign', 'benign', 'not_in_clinvar']
        clinvar_aliases = {
            'pathogenic': 'path', 'likely_pathogenic': 'lpath', 'vus': 'vus',
            'likely_benign': 'lben', 'benign': 'ben', 'not_in_clinvar': 'none'
        }
        
        return self._prompt_choice(
            "ClinVar status or press Enter to skip: ",
            choices=clinvar_choices,
            required=False,
            aliases=clinvar_aliases
        )
    
    def collect_population_data(self) -> Dict[str, Any]:
        """
        Collect population frequency data.
        
        Returns:
            Dict[str, Any]: Population frequency data
        """
        print("\n🌍 POPULATION FREQUENCY DATA")
        print("-" * 30)
        print("Enter allele frequencies from population databases.")
        print("💡 Use gnomAD.broadinstitute.org or Varsome.com for quick lookup.")
        print("💡 Lower frequencies support pathogenicity.")
        print("💡 Use 'NA' or leave empty if data is not available.")
        
        if self.test_mode:
            print("🧪 Using test mode sample data...")
            return TEST_MODE_DATA['population_data'].copy()
        
        population_data = {}
        
        # gnomAD overall frequency
        population_data['gnomad_af'] = self._prompt_input(
            "gnomAD overall allele frequency (0-1): ",
            validator=self._validate_allele_frequency,
            required=False,
            convert_type=float
        )
        
        # gnomAD popmax frequency
        population_data['gnomad_af_popmax'] = self._prompt_input(
            "gnomAD popmax allele frequency (0-1): ",
            validator=self._validate_allele_frequency,
            required=False,
            convert_type=float
        )
        
        # Homozygous count
        population_data['gnomad_hom_count'] = self._prompt_input(
            "gnomAD homozygous count: ",
            validator=self._validate_count,
            required=False,
            convert_type=int
        )
        
        # Heterozygous count
        population_data['gnomad_het_count'] = self._prompt_input(
            "gnomAD heterozygous count: ",
            validator=self._validate_count,
            required=False,
            convert_type=int
        )
        
        # ExAC frequency
        population_data['exac_af'] = self._prompt_input(
            "ExAC allele frequency (0-1): ",
            validator=self._validate_allele_frequency,
            required=False,
            convert_type=float
        )
        
        # Disease prevalence
        population_data['disease_prevalence'] = self._prompt_input(
            "Disease prevalence (e.g., 1/10000): ",
            required=False
        )
        
        return population_data
    
    def collect_insilico_data(self, variant_type: str) -> Dict[str, Any]:
        """
        Collect in silico prediction scores.
        
        Args:
            variant_type (str): Type of variant
            
        Returns:
            Dict[str, Any]: In silico prediction data
        """
        print("\n🤖 IN SILICO PREDICTION SCORES")
        print("-" * 30)
        print("Enter prediction scores from various algorithms.")
        print("💡 Use Varsome.com or similar tools to get these scores quickly.")
        print("💡 Leave empty if score is not available.")
        print("💡 Higher scores usually = more pathogenic (except SIFT)")
        print("💡 You can skip rarely used predictors if time is limited.")
        print("💡 Some predictors have both raw and ranked scores - choose accordingly.")
        
        if self.test_mode:
            print("🧪 Using test mode sample data...")
            return TEST_MODE_DATA['insilico_data'].copy()
        
        insilico_data = {}
        
        # REVEL score (0-1)
        insilico_data['revel'] = self._prompt_input(
            "REVEL score (0-1): ",
            validator=self._validate_score_0_1,
            required=False,
            convert_type=float
        )
        
        # CADD phred score only (more commonly used)
        insilico_data['cadd_phred'] = self._prompt_input(
            "CADD phred score (0-99): ",
            validator=self._validate_cadd_phred_score,
            required=False,
            convert_type=float
        )
        
        # Conservation scores
        insilico_data['gerp_pp'] = self._prompt_input(
            "GERP++ score (-12.3 to 6.17): ",
            validator=self._validate_gerp_score,
            required=False,
            convert_type=float
        )
        
        # PhyloP scores - multiple types
        print("\n📊 PhyloP Conservation Scores:")
        print("   These scores measure evolutionary conservation")
        print("   Ranked scores (0-1) are preferred for analysis")
        print("   Choose ONE scoring type that will apply to all PhyloP scores:")
        print()
        print("   1. raw_score (or 'raw', 'r') - Raw conservation scores (-20 to 10)")
        print("   2. ranked_score (or 'ranked', 'rank') - Ranked scores (0-1)")
        print("   3. skip (or 's') - Skip all PhyloP scores")
        
        # Ask once for all PhyloP scores
        phylop_choice = self._prompt_choice(
            "Choose PhyloP scoring type for ALL scores: ",
            choices=['raw_score', 'ranked_score', 'skip'],
            required=False,
            aliases={
                '1': 'raw_score', 'raw': 'raw_score', 'r': 'raw_score',
                '2': 'ranked_score', 'ranked': 'ranked_score', 'rank': 'ranked_score',
                '3': 'skip', 's': 'skip'
            }
        )
        
        if phylop_choice == 'raw_score':
            print(f"{COLORAMA_COLORS['GREEN']}✅ Using raw PhyloP scores for all{COLORAMA_COLORS['RESET']}")
            
            # PhyloP vertebrates raw
            insilico_data['phylop_vert'] = self._prompt_input(
                "PhyloP vertebrates raw score (-20 to 10): ",
                validator=self._validate_phylop_score,
                required=False,
                convert_type=float
            )
            
            # PhyloP mammals raw
            insilico_data['phylop_mamm'] = self._prompt_input(
                "PhyloP mammals raw score (-20 to 10): ",
                validator=self._validate_phylop_score,
                required=False,
                convert_type=float
            )
            
            # PhyloP primates raw
            insilico_data['phylop_prim'] = self._prompt_input(
                "PhyloP primates raw score (-20 to 10): ",
                validator=self._validate_phylop_score,
                required=False,
                convert_type=float
            )
            
        elif phylop_choice == 'ranked_score':
            print(f"{COLORAMA_COLORS['GREEN']}✅ Using ranked PhyloP scores for all{COLORAMA_COLORS['RESET']}")
            
            # PhyloP vertebrates ranked
            insilico_data['phylop_vert_ranked'] = self._prompt_input(
                "PhyloP vertebrates ranked score (0-1): ",
                validator=self._validate_score_0_1,
                required=False,
                convert_type=float
            )
            
            # PhyloP mammals ranked
            insilico_data['phylop_mamm_ranked'] = self._prompt_input(
                "PhyloP mammals ranked score (0-1): ",
                validator=self._validate_score_0_1,
                required=False,
                convert_type=float
            )
            
            # PhyloP primates ranked
            insilico_data['phylop_prim_ranked'] = self._prompt_input(
                "PhyloP primates ranked score (0-1): ",
                validator=self._validate_score_0_1,
                required=False,
                convert_type=float
            )
        else:
            print(f"{COLORAMA_COLORS['YELLOW']}⏭️ Skipping all PhyloP scores{COLORAMA_COLORS['RESET']}")
        
        # AlphaMissense score (missense variants only)
        if variant_type == 'missense':
            insilico_data['alphamissense'] = self._prompt_input(
                "AlphaMissense score (0-1): ",
                validator=self._validate_score_0_1,
                required=False,
                convert_type=float
            )
        
        # MetaRNN - choice between score and ranked score
        print("\n📊 MetaRNN options:")
        print("   1. Score (0-1)")
        print("   2. Ranked score (0-1)")
        metarnn_choice = self._prompt_choice(
            "Choose MetaRNN type (or leave empty to skip): ",
            choices=['score', 'ranked_score', 'skip'],
            required=False,
            aliases={'1': 'score', '2': 'ranked_score', 's': 'skip'}
        )
        
        if metarnn_choice == 'score':
            insilico_data['metarnn'] = self._prompt_input(
                "MetaRNN score (0-1): ",
                validator=self._validate_score_0_1,
                required=False,
                convert_type=float
            )
        elif metarnn_choice == 'ranked_score':
            insilico_data['metarnn_ranked'] = self._prompt_input(
                "MetaRNN ranked score (0-1): ",
                validator=self._validate_score_0_1,
                required=False,
                convert_type=float
            )
        
        # ClinPred score
        insilico_data['clinpred'] = self._prompt_input(
            "ClinPred score (0-1): ",
            validator=self._validate_score_0_1,
            required=False,
            convert_type=float
        )
        
        # BayesDel - choice between max AF and no max AF
        print("\n📊 BayesDel options:")
        print("   1. BayesDel addAF (with max AF)")
        print("   2. BayesDel noAF (without max AF)")
        bayesdel_choice = self._prompt_choice(
            "Choose BayesDel type (or leave empty to skip): ",
            choices=['addAF', 'noAF', 'skip'],
            required=False,
            aliases={'1': 'addAF', '2': 'noAF', 's': 'skip'}
        )
        
        if bayesdel_choice == 'addAF':
            insilico_data['bayesdel_addaf'] = self._prompt_input(
                "BayesDel addAF score (-1 to 1): ",
                validator=self._validate_bayesdel_score,
                required=False,
                convert_type=float
            )
        elif bayesdel_choice == 'noAF':
            insilico_data['bayesdel_noaf'] = self._prompt_input(
                "BayesDel noAF score (-1 to 1): ",
                validator=self._validate_bayesdel_score,
                required=False,
                convert_type=float
            )
        
        # SIFT score (missense variants only)
        if variant_type == 'missense':
            insilico_data['sift'] = self._prompt_input(
                "SIFT score (0-1, lower = more pathogenic): ",
                validator=self._validate_score_0_1,
                required=False,
                convert_type=float
            )
        
        # PolyPhen-2 score (missense variants only)
        if variant_type == 'missense':
            insilico_data['polyphen2'] = self._prompt_input(
                "PolyPhen-2 score (0-1): ",
                validator=self._validate_score_0_1,
                required=False,
                convert_type=float
            )
        
        # MutationTaster - choice between score and ranked score
        print("\n📊 MutationTaster options:")
        print("   1. Score (0-1)")
        print("   2. Ranked score (0-1)")
        mutationtaster_choice = self._prompt_choice(
            "Choose MutationTaster type (or leave empty to skip): ",
            choices=['score', 'ranked_score', 'skip'],
            required=False,
            aliases={'1': 'score', '2': 'ranked_score', 's': 'skip'}
        )
        
        if mutationtaster_choice == 'score':
            insilico_data['mutationtaster'] = self._prompt_input(
                "MutationTaster score (0-1): ",
                validator=self._validate_score_0_1,
                required=False,
                convert_type=float
            )
        elif mutationtaster_choice == 'ranked_score':
            insilico_data['mutationtaster_ranked'] = self._prompt_input(
                "MutationTaster ranked score (0-1): ",
                validator=self._validate_score_0_1,
                required=False,
                convert_type=float
            )
        
        # FATHMM score (optional)
        print("\n📊 FATHMM options:")
        print("   1. Score (-16.13 to 10.64)")
        print("   2. Ranked score (0-1)")
        fathmm_choice = self._prompt_choice(
            "Choose FATHMM type (or leave empty to skip): ",
            choices=['score', 'ranked_score', 'skip'],
            required=False,
            aliases={'1': 'score', '2': 'ranked_score', 's': 'skip'}
        )
        
        if fathmm_choice == 'score':
            insilico_data['fathmm'] = self._prompt_input(
                "FATHMM score (-16.13 to 10.64): ",
                validator=self._validate_fathmm_score,
                required=False,
                convert_type=float
            )
        elif fathmm_choice == 'ranked_score':
            insilico_data['fathmm_ranked'] = self._prompt_input(
                "FATHMM ranked score (0-1): ",
                validator=self._validate_score_0_1,
                required=False,
                convert_type=float
            )
        
        # HIGH PRIORITY PREDICTORS - Enhanced coverage
        print(f"\n{COLORAMA_COLORS['GREEN']}🔬 HIGH PRIORITY PREDICTORS{COLORAMA_COLORS['RESET']}")
        print("These predictors significantly improve classification accuracy:")
        
        # VEST4 (Cancer-specific)
        insilico_data['vest4'] = self._prompt_input(
            "VEST4 score (0-1, cancer-specific): ",
            validator=self._validate_score_0_1,
            required=False,
            convert_type=float
        )
        
        # PrimateAI (State-of-the-art)
        insilico_data['primateai'] = self._prompt_input(
            "PrimateAI score (0-1, deep learning): ",
            validator=self._validate_score_0_1,
            required=False,
            convert_type=float
        )
        
        # ESM1b (Protein language model)
        insilico_data['esm1b'] = self._prompt_input(
            "ESM1b score (protein language model): ",
            validator=self._validate_score_0_1,
            required=False,
            convert_type=float
        )
        
        # PROVEAN (Protein variation effect)
        insilico_data['provean'] = self._prompt_input(
            "PROVEAN score (typically -14 to 14, <-2.5 deleterious): ",
            validator=self._validate_provean_score,
            required=False,
            convert_type=float
        )
        
        # Enhanced splice prediction for relevant variants
        if variant_type in ['splice', 'intronic', 'synonymous']:
            print(f"\n{COLORAMA_COLORS['YELLOW']}✂️ ENHANCED SPLICE PREDICTION{COLORAMA_COLORS['RESET']}")
            print("Advanced splice impact prediction:")
            
            # MMSplice (Comprehensive splice modeling)
            insilico_data['mmsplice'] = self._prompt_input(
                "MMSplice score (delta logit PSI): ",
                validator=self._validate_mmsplice_score,
                required=False,
                convert_type=float
            )
        
        return insilico_data
    
    def collect_genetic_data(self) -> Dict[str, Any]:
        """
        Collect genetic and inheritance data.
        
        Returns:
            Dict[str, Any]: Genetic data
        """
        print("\n🧬 GENETIC DATA")
        print("-" * 30)
        
        if self.test_mode:
            print("🧪 Using test mode sample data...")
            return TEST_MODE_DATA['genetic_data'].copy()
        
        genetic_data = {}
        
        # Inheritance pattern
        print("📊 Inheritance Pattern Information:")
        print("   AD = Autosomal Dominant")
        print("   AR = Autosomal Recessive")
        print("   XLD = X-linked Dominant")
        print("   XLR = X-linked Recessive")
        print("   MT = Mitochondrial")
        
        inheritance_choices = ['AD', 'AR', 'XLD', 'XLR', 'MT', 'unknown']
        inheritance_aliases = {
            'AD': 'ad', 'AR': 'ar', 'XLD': 'xld', 'XLR': 'xlr', 'MT': 'mt'
        }
        genetic_data['inheritance'] = self._prompt_choice(
            "Inheritance pattern: ",
            choices=inheritance_choices,
            required=True,
            aliases=inheritance_aliases
        )
        
        # Zygosity
        print("\n📊 Zygosity Information:")
        print("   Homozygous = Both alleles affected (hom)")
        print("   Heterozygous = One allele affected (het)")
        print("   Hemizygous = X-linked in males (hemi)")
        
        zygosity_choices = ['homozygous', 'heterozygous', 'hemizygous', 'unknown']
        zygosity_aliases = {
            'homozygous': 'hom', 'heterozygous': 'het', 'hemizygous': 'hemi'
        }
        genetic_data['zygosity'] = self._prompt_choice(
            "Zygosity: ",
            choices=zygosity_choices,
            required=True,
            aliases=zygosity_aliases
        )
        
        # Phase information
        if genetic_data['inheritance'] == 'AR' and genetic_data['zygosity'] == 'heterozygous':
            genetic_data['compound_het'] = self._prompt_choice(
                "Is this part of compound heterozygous variants? ",
                choices=['yes', 'no'],
                required=True
            )
        
        # Allelic data for recessive diseases
        if genetic_data['inheritance'] == 'AR':
            genetic_data['allelic_data'] = self._prompt_input(
                "Allelic data (describe second variant if compound het): ",
                required=False
            )
        
        return genetic_data
    
    def collect_functional_data(self) -> Dict[str, Any]:
        """
        Collect functional studies and segregation data.
        
        Returns:
            Dict[str, Any]: Functional data
        """
        print("\n🔬 FUNCTIONAL DATA")
        print("-" * 30)
        
        if self.test_mode:
            print("🧪 Using test mode sample data...")
            return TEST_MODE_DATA['functional_data'].copy()
        
        functional_data = {}
        
        # Segregation analysis
        print("📊 Segregation Analysis Information:")
        print("   Cosegregates = Variant tracks with disease in family")
        print("   Does not segregate = Variant does not track with disease")
        
        segregation_choices = ['cosegregates', 'does_not_segregate', 'insufficient_data', 'not_performed']
        segregation_aliases = {
            'cosegregates': 'coseg', 'does_not_segregate': 'noseg', 
            'insufficient_data': 'insuf', 'not_performed': 'none'
        }
        functional_data['segregation'] = self._prompt_choice(
            "Segregation analysis: ",
            choices=segregation_choices,
            required=False,
            aliases=segregation_aliases
        )
        
        # De novo status
        print("\n📊 De Novo Status Information:")
        print("   Confirmed = Parental testing confirms de novo")
        print("   Assumed = Likely de novo but not confirmed")
        print("   Not de novo = Inherited from parent")
        
        denovo_choices = ['confirmed', 'assumed', 'not_denovo', 'unknown']
        denovo_aliases = {
            'confirmed': 'conf', 'assumed': 'assum', 'not_denovo': 'inherit'
        }
        functional_data['denovo'] = self._prompt_choice(
            "De novo status: ",
            choices=denovo_choices,
            required=False,
            aliases=denovo_aliases
        )
        
        # Functional studies
        print("\n📊 Functional Studies Information:")
        print("   Damaging = Studies show harmful effect")
        print("   Benign = Studies show no harmful effect")
        print("   Inconclusive = Studies are unclear")
        
        functional_choices = ['damaging', 'benign', 'inconclusive', 'not_performed']
        functional_aliases = {
            'damaging': 'dam', 'benign': 'ben', 'inconclusive': 'inc', 'not_performed': 'none'
        }
        functional_data['functional_studies'] = self._prompt_choice(
            "Functional studies result: ",
            choices=functional_choices,
            required=False,
            aliases=functional_aliases
        )
        
        # Phenotype match
        print("\n📊 Phenotype Match Information:")
        print("   Specific match = Phenotype highly specific for this gene")
        print("   Partial match = Some phenotype overlap")
        print("   No match = Phenotype doesn't match gene")
        
        phenotype_choices = ['specific_match', 'partial_match', 'no_match', 'unknown']
        phenotype_aliases = {
            'specific_match': 'spec', 'partial_match': 'part', 'no_match': 'none'
        }
        functional_data['phenotype_match'] = self._prompt_choice(
            "Phenotype match: ",
            choices=phenotype_choices,
            required=False,
            aliases=phenotype_aliases
        )
        
        # Case-control data
        print("\n📊 Case-Control Data Information:")
        print("   This is for research studies comparing patients vs controls")
        print("   Helps establish statistical significance (PS4 criterion)")
        
        functional_data['case_control'] = self._prompt_choice(
            "Case-control data available? ",
            choices=['yes', 'no'],
            aliases={'yes': 'y', 'no': 'n'},
            required=False
        )
        
        if functional_data['case_control'] == 'yes':
            functional_data['cases_with_variant'] = self._prompt_input(
                "Number of cases with variant: ",
                validator=self._validate_count,
                required=False,
                convert_type=int
            )
            
            functional_data['total_cases'] = self._prompt_input(
                "Total number of cases: ",
                validator=self._validate_count,
                required=False,
                convert_type=int
            )
            
            functional_data['controls_with_variant'] = self._prompt_input(
                "Number of controls with variant: ",
                validator=self._validate_count,
                required=False,
                convert_type=int
            )
            
            functional_data['total_controls'] = self._prompt_input(
                "Total number of controls: ",
                validator=self._validate_count,
                required=False,
                convert_type=int
            )
        
        return functional_data
    
    def _prompt_input(self, prompt: str, validator: Optional[Callable] = None, 
                     required: bool = False, convert_type: Optional[type] = None) -> Any:
        """
        Prompt user for input with validation.
        
        Args:
            prompt (str): Prompt message
            validator (callable, optional): Validation function
            required (bool): Whether input is required
            convert_type (type, optional): Type to convert input to
            
        Returns:
            Any: User input (validated and converted if specified)
        """
        while True:
            try:
                value = input(prompt).strip()
                
                # Handle empty input
                if not value:
                    if required:
                        print("❌ This field is required. Please enter a value.")
                        continue
                    else:
                        return None
                
                # Handle 'NA' input
                if value.upper() in ['NA', 'N/A', 'NOT AVAILABLE']:
                    return None
                
                # Validate input
                if validator:
                    if not validator(value):
                        # Give specific error messages for common validators
                        if validator.__name__ == '_validate_position':
                            print("❌ Please enter a valid genomic position (positive integer, e.g., 41276045)")
                        elif validator.__name__ == '_validate_gene_symbol':
                            print("❌ Please enter a valid gene symbol (e.g., BRCA1, TP53)")
                        elif validator.__name__ == '_validate_chromosome':
                            print("❌ Please enter a valid chromosome (1-22, X, Y, or MT)")
                        elif validator.__name__ == '_validate_allele':
                            print("❌ Please enter valid DNA bases (A, T, C, G only)")
                        elif 'score' in validator.__name__:
                            print("❌ Please enter a valid score within the specified range")
                        else:
                            print("❌ Invalid input format. Please check and try again.")
                        continue
                
                # Convert type if specified
                if convert_type:
                    try:
                        value = convert_type(value)
                    except ValueError:
                        print(f"❌ Invalid {convert_type.__name__} format. Please try again.")
                        continue
                
                return value
                
            except KeyboardInterrupt:
                print("\n❌ Operation cancelled.")
                sys.exit(1)
            except Exception as e:
                print(f"❌ Error: {e}")
                continue
    
    def _prompt_choice(self, prompt: str, choices: List[str], required: bool = False, 
                      aliases: Optional[Dict[str, str]] = None) -> Optional[str]:
        """
        Prompt user to choose from a list of options with alias support.
        
        Args:
            prompt (str): Prompt message
            choices (List[str]): List of valid choices
            required (bool): Whether choice is required
            aliases (dict, optional): Dictionary of aliases for choices
            
        Returns:
            Optional[str]: Selected choice
        """
        # For long choice lists, don't display all options repeatedly
        if len(choices) > 10:
            print(f"\n{COLORAMA_COLORS['BLUE']}💡 Available options shown above. Type your choice below.{COLORAMA_COLORS['RESET']}")
        else:
            # Create display text with aliases for short lists
            display_choices = []
            for choice in choices:
                if aliases and choice in aliases:
                    display_choices.append(f"{choice} ({aliases[choice]})")
                else:
                    display_choices.append(choice)
            
            print(f"\nOptions: {', '.join(display_choices)}")
            if aliases:
                print("💡 You can use short forms shown in parentheses")
        
        while True:
            try:
                value = input(f"{prompt}").strip().lower()
                
                if not value:
                    if required:
                        print("❌ This field is required. Please select an option.")
                        continue
                    else:
                        return None
                
                # Check direct matches first
                if value in [choice.lower() for choice in choices]:
                    return next(choice for choice in choices if choice.lower() == value)
                
                # Check aliases
                if aliases:
                    for choice, alias in aliases.items():
                        if value == alias.lower():
                            return choice
                
                # For long lists, show only first few options in error message
                if len(choices) > 10:
                    common_choices = choices[:8]  # Show first 8 options
                    print(f"❌ Invalid choice. Common options: {', '.join(common_choices)}, ...")
                    print(f"   Or type 'other' if not in the list above.")
                else:
                    # Show all options for short lists
                    display_choices = []
                    for choice in choices:
                        if aliases and choice in aliases:
                            display_choices.append(f"{choice} ({aliases[choice]})")
                        else:
                            display_choices.append(choice)
                    print(f"❌ Invalid choice. Please select from: {', '.join(display_choices)}")
                continue
                    
            except KeyboardInterrupt:
                print("\n❌ Operation cancelled.")
                sys.exit(1)
    
    # Validation methods
    def _validate_gene_symbol(self, value: str) -> bool:
        """Validate gene symbol format."""
        return bool(re.match(VALIDATION_PATTERNS['gene_symbol'], value.upper()))
    
    def _validate_chromosome(self, value: str) -> bool:
        """Validate chromosome format."""
        if not value:
            return False
        
        # Remove 'chr' prefix if present
        clean_value = value.upper().replace('CHR', '')
        
        # Check if it matches the pattern
        valid_chromosomes = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', 
                           '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', 
                           '21', '22', 'X', 'Y', 'MT', 'M']
        
        return clean_value in valid_chromosomes
    
    def _validate_position(self, value: str) -> bool:
        """Validate genomic position."""
        try:
            # Remove any non-digit characters except for scientific notation
            cleaned_value = ''.join(c for c in value if c.isdigit())
            if not cleaned_value:
                return False
            pos = int(cleaned_value)
            return pos > 0 and pos < 1000000000  # Reasonable range for genomic positions
        except ValueError:
            return False
    
    def _validate_allele(self, value: str) -> bool:
        """Validate allele format."""
        return bool(re.match(r'^[ATCG]+$', value.upper()))
    
    def _validate_hgvs_cdna(self, value: str) -> bool:
        """Validate HGVS cDNA format."""
        if not value:
            return True  # Optional field
        return bool(re.match(VALIDATION_PATTERNS['hgvs_cdna'], value))
    
    def _validate_hgvs_protein(self, value: str) -> bool:
        """Validate HGVS protein format."""
        if not value:
            return True  # Optional field
        
        # Common patterns for protein variants
        patterns = [
            # Three-letter amino acid code: p.Arg123Gln, p.Met107Val
            r'^p\.[A-Z][a-z]{2}\d+[A-Z][a-z]{2}$',
            # Three-letter code with stop: p.Arg123*
            r'^p\.[A-Z][a-z]{2}\d+\*$',
            # Single letter code: p.R123Q  
            r'^p\.[A-Z]\d+[A-Z]$',
            # Single letter with stop: p.R123*
            r'^p\.[A-Z]\d+\*$',
            # Stop to amino acid: p.*123Arg
            r'^p\.\*\d+[A-Z][a-z]{2}$',
            # Three letter synonymous: p.Arg123=
            r'^p\.[A-Z][a-z]{2}\d+=$',
            # Single letter synonymous: p.R123=
            r'^p\.[A-Z]\d+=$',
            # Frameshift: p.Arg123fs or p.Arg123Alafs*10
            r'^p\.[A-Z][a-z]{2}\d+([A-Z][a-z]{2})?fs(\*\d+)?$',
            # Extension: p.*123Argext*10
            r'^p\.\*\d+[A-Z][a-z]{2}ext\*\d+$'
        ]
        
        for pattern in patterns:
            if re.match(pattern, value):
                return True
        
        return False
    
    def _validate_allele_frequency(self, value: str) -> bool:
        """Validate allele frequency."""
        try:
            freq = float(value)
            return 0 <= freq <= 1
        except ValueError:
            return False
    
    def _validate_count(self, value: str) -> bool:
        """Validate count value."""
        try:
            count = int(value)
            return count >= 0
        except ValueError:
            return False
    
    def _validate_score_0_1(self, value: str) -> bool:
        """Validate score between 0 and 1."""
        try:
            score = float(value)
            return 0 <= score <= 1
        except ValueError:
            return False
    
    def _validate_cadd_score(self, value: str) -> bool:
        """Validate CADD score."""
        try:
            score = float(value)
            return 0 <= score <= 50
        except ValueError:
            return False
    
    def _validate_cadd_phred_score(self, value: str) -> bool:
        """Validate CADD phred score."""
        try:
            score = float(value)
            return 0 <= score <= 99
        except ValueError:
            return False
    
    def _validate_gerp_score(self, value: str) -> bool:
        """Validate GERP++ score."""
        try:
            score = float(value)
            return -12.3 <= score <= 6.17
        except ValueError:
            return False
    
    def _validate_phylop_score(self, value: str) -> bool:
        """Validate PhyloP score."""
        try:
            score = float(value)
            return -20 <= score <= 10
        except ValueError:
            return False
    
    def _validate_fathmm_score(self, value: str) -> bool:
        """Validate FATHMM score."""
        try:
            score = float(value)
            return -16.13 <= score <= 10.64
        except ValueError:
            return False
    
    def _validate_bayesdel_score(self, value: str) -> bool:
        """Validate BayesDel score."""
        try:
            score = float(value)
            return -1 <= score <= 1
        except ValueError:
            return False
    
    def _validate_provean_score(self, value: str) -> bool:
        """Validate PROVEAN score."""
        try:
            score = float(value)
            return -14 <= score <= 14
        except ValueError:
            return False
    
    def _validate_mmsplice_score(self, value: str) -> bool:
        """Validate MMSplice score."""
        try:
            score = float(value)
            return -10 <= score <= 10  # Delta logit PSI range
        except ValueError:
            return False

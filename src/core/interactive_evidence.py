"""
Interactive Evidence Collector Module
=====================================

This module handles collection of literature-based, judgment-heavy ACMG evidence
criteria through interactive user input. These criteria require human interpretation
of clinical and research data, and should NOT be inferred automatically.

Covered criteria:
- PS3 / BS3: Functional studies (well-established in vitro/in vivo)
- PS4: Case-control / enrichment data
- PP1 / BS4: Segregation with disease in families
- PS1 / PM5: Same codon / same amino acid change as known pathogenic variant
- PP5 / BP6: Reputable source / conflicting clinical assertions

IMPORTANT PHILOSOPHY:
The program's job is ONLY to:
1. Ask structured questions
2. Transform user answers into ACMG evidence codes and strengths
3. Combine them with automatically derived evidence

It does NOT perform automatic literature mining, PubMed querying, or 
hard-coded database lookups for these criteria.

Author: Can Sevilmiş
License: MIT
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Callable, Any, Tuple
from enum import Enum


# =============================================================================
# Evidence Strength Levels (ACMG)
# =============================================================================
class EvidenceStrength(Enum):
    """ACMG evidence strength levels."""
    VERY_STRONG = "Very Strong"
    STRONG = "Strong"
    MODERATE = "Moderate"
    SUPPORTING = "Supporting"


# =============================================================================
# Manual Evidence Data Structure
# =============================================================================
@dataclass
class ManualEvidence:
    """
    Data structure representing user-provided ACMG evidence.
    
    Attributes:
        codes: List of evidence codes with strength suffixes 
               (e.g., ["PS3_strong", "PP1_supporting"])
        explanations: Dictionary mapping each code to a human-readable explanation
                     describing the basis for the evidence assignment
    """
    codes: List[str] = field(default_factory=list)
    explanations: Dict[str, str] = field(default_factory=dict)
    
    def add_evidence(self, code: str, explanation: str) -> None:
        """
        Add an evidence code with its explanation.
        
        Args:
            code: ACMG evidence code (e.g., "PS3_strong")
            explanation: Human-readable explanation for the evidence
        """
        if code not in self.codes:
            self.codes.append(code)
        self.explanations[code] = explanation
    
    def has_evidence(self) -> bool:
        """Check if any manual evidence has been collected."""
        return len(self.codes) > 0
    
    def get_pathogenic_codes(self) -> List[str]:
        """Get all pathogenic evidence codes (P*, PM*, PP*)."""
        return [c for c in self.codes if c.startswith(('PVS', 'PS', 'PM', 'PP'))]
    
    def get_benign_codes(self) -> List[str]:
        """Get all benign evidence codes (BA*, BS*, BP*)."""
        return [c for c in self.codes if c.startswith(('BA', 'BS', 'BP'))]
    
    def merge_with(self, other: 'ManualEvidence') -> 'ManualEvidence':
        """Merge this evidence with another ManualEvidence object."""
        merged = ManualEvidence()
        merged.codes = list(set(self.codes + other.codes))
        merged.explanations = {**self.explanations, **other.explanations}
        return merged


# =============================================================================
# Thresholds and Constants for Evidence Mapping
# =============================================================================
# These thresholds are EDUCATIONAL APPROXIMATIONS based on ACMG guidelines.
# Adjust as needed for specific clinical contexts.

# PS3/BS3: Functional Studies Thresholds
# Based on number of independent studies and quality
FUNCTIONAL_STUDY_THRESHOLDS = {
    'PS3_strong': {'min_damaging_studies': 2, 'min_quality': 'high'},
    'PS3_moderate': {'min_damaging_studies': 1, 'min_quality': 'high'},
    'PS3_supporting': {'min_damaging_studies': 1, 'min_quality': 'moderate'},
    'BS3_strong': {'min_benign_studies': 2, 'min_quality': 'high'},
    'BS3_moderate': {'min_benign_studies': 1, 'min_quality': 'high'},
    'BS3_supporting': {'min_benign_studies': 1, 'min_quality': 'moderate'},
}

# PP1/BS4: Segregation Thresholds
# Based on number of meioses (informative family members)
# Reference: Jarvik & Browning, 2016 - Cosegregation evidence
SEGREGATION_THRESHOLDS = {
    'PP1_strong': {'min_affected_carriers': 5},
    'PP1_moderate': {'min_affected_carriers': 3},
    'PP1_supporting': {'min_affected_carriers': 2},
    'BS4_supporting': {'min_unaffected_carriers': 2},
    'BS4_moderate': {'min_unaffected_carriers': 3},
}

# PS4: Case-Control Thresholds
# Based on odds ratio and statistical significance
PS4_THRESHOLDS = {
    'PS4_strong': {'min_cases': 5, 'min_or': 5.0},
    'PS4_moderate': {'min_cases': 3, 'min_or': 3.0},
    'PS4_supporting': {'min_cases': 2, 'min_or': 2.0},
}


# =============================================================================
# Input Provider Abstraction (for testability)
# =============================================================================
class InputProvider:
    """
    Abstraction for input collection, enabling testability.
    
    In production, this uses real input() calls.
    In tests, a mock can provide scripted responses.
    """
    
    def prompt(self, message: str) -> str:
        """Prompt user for input and return the response."""
        return input(message).strip()
    
    def prompt_yes_no(self, message: str, default: Optional[bool] = None) -> bool:
        """
        Prompt user for yes/no response.
        
        Args:
            message: Prompt message
            default: Default value if user presses Enter (None = require input)
            
        Returns:
            True for yes, False for no
        """
        suffix = " [y/n]"
        if default is True:
            suffix = " [Y/n]"
        elif default is False:
            suffix = " [y/N]"
        
        while True:
            response = self.prompt(f"{message}{suffix}: ").lower()
            
            if not response and default is not None:
                return default
            
            if response in ('y', 'yes'):
                return True
            elif response in ('n', 'no'):
                return False
            
            print("Please enter 'y' or 'n'.")
    
    def prompt_choice(self, message: str, choices: List[str], 
                      default: Optional[str] = None) -> str:
        """
        Prompt user to select from choices.
        
        Args:
            message: Prompt message
            choices: List of valid choices
            default: Default choice if user presses Enter
            
        Returns:
            Selected choice
        """
        choices_str = "/".join(choices)
        if default:
            choices_str = choices_str.replace(default, f"[{default}]")
        
        while True:
            response = self.prompt(f"{message} ({choices_str}): ").lower()
            
            if not response and default:
                return default
            
            if response in [c.lower() for c in choices]:
                return response
            
            print(f"Please select from: {', '.join(choices)}")
    
    def prompt_integer(self, message: str, min_val: int = 0, 
                       max_val: int = 1000, default: Optional[int] = None) -> Optional[int]:
        """
        Prompt user for an integer value.
        
        Args:
            message: Prompt message
            min_val: Minimum valid value
            max_val: Maximum valid value
            default: Default value if user presses Enter
            
        Returns:
            Integer value or None if skipped
        """
        suffix = f" ({min_val}-{max_val})"
        if default is not None:
            suffix += f" [default: {default}]"
        
        while True:
            response = self.prompt(f"{message}{suffix}: ")
            
            if not response:
                return default
            
            try:
                value = int(response)
                if min_val <= value <= max_val:
                    return value
                print(f"Please enter a value between {min_val} and {max_val}.")
            except ValueError:
                print("Please enter a valid integer.")


class MockInputProvider(InputProvider):
    """
    Mock input provider for testing.
    
    Provides scripted responses to prompts.
    """
    
    def __init__(self, responses: List[str]):
        """
        Initialize with list of responses to provide.
        
        Args:
            responses: List of responses in order they will be requested
        """
        self.responses = responses
        self.index = 0
    
    def prompt(self, message: str) -> str:
        """Return next scripted response."""
        if self.index < len(self.responses):
            response = self.responses[self.index]
            self.index += 1
            return response
        return ""


# =============================================================================
# Interactive Evidence Collector
# =============================================================================
class InteractiveEvidenceCollector:
    """
    Collects literature-based ACMG evidence through interactive user input.
    
    This class handles all judgment-heavy criteria that require human 
    interpretation of clinical/research data:
    - PS3/BS3: Functional studies
    - PS4: Case-control data
    - PP1/BS4: Segregation
    - PS1/PM5: Same codon/amino acid change
    - PP5/BP6: Reputable source
    
    The collector is designed to be:
    1. Conservative: Prefer "no evidence" over over-calling
    2. Transparent: Show clear reasoning for evidence assignment
    3. Testable: Input can be mocked for unit testing
    
    Example:
        >>> collector = InteractiveEvidenceCollector()
        >>> evidence = collector.collect_all()
        >>> print(evidence.codes)
        ['PS3_moderate', 'PP1_supporting']
    """
    
    def __init__(self, input_provider: Optional[InputProvider] = None,
                 show_prompts: bool = True):
        """
        Initialize the collector.
        
        Args:
            input_provider: Custom input provider (for testing)
            show_prompts: Whether to show detailed prompts and help text
        """
        self.input = input_provider or InputProvider()
        self.show_prompts = show_prompts
        self.evidence = ManualEvidence()
    
    def _print_header(self, title: str) -> None:
        """Print a section header."""
        if self.show_prompts:
            print(f"\n{'='*60}")
            print(f"📋 {title}")
            print(f"{'='*60}")
    
    def _print_info(self, message: str) -> None:
        """Print an informational message."""
        if self.show_prompts:
            print(f"ℹ️  {message}")
    
    def _print_success(self, message: str) -> None:
        """Print a success message."""
        if self.show_prompts:
            print(f"✅ {message}")
    
    def _print_note(self, message: str) -> None:
        """Print a note/tip."""
        if self.show_prompts:
            print(f"💡 {message}")
    
    # =========================================================================
    # PS3 / BS3: Functional Studies
    # =========================================================================
    def collect_ps3_bs3(self) -> ManualEvidence:
        """
        Collect PS3/BS3 evidence based on functional studies.
        
        PS3: Well-established in vitro or in vivo functional studies supportive
             of a damaging effect on the gene or gene product.
        BS3: Well-established in vitro or in vivo functional studies show no
             damaging effect on protein function or splicing.
        
        Returns:
            ManualEvidence with PS3 or BS3 codes if applicable
        """
        self._print_header("PS3/BS3: Functional Studies Evidence")
        self._print_info("Functional studies include: enzyme assays, cell-based assays,")
        self._print_info("animal models, splicing assays, protein stability tests, etc.")
        self._print_note("Only include published, peer-reviewed studies.")
        
        has_studies = self.input.prompt_yes_no(
            "Are there published functional studies for this specific variant?",
            default=False
        )
        
        if not has_studies:
            self._print_info("No functional studies → No PS3/BS3 evidence")
            return ManualEvidence()
        
        # Count damaging studies
        damaging_count = self.input.prompt_integer(
            "How many independent studies show DAMAGING effect?",
            min_val=0, max_val=50, default=0
        ) or 0
        
        # Count benign studies
        benign_count = self.input.prompt_integer(
            "How many independent studies show NO EFFECT (benign)?",
            min_val=0, max_val=50, default=0
        ) or 0
        
        # Quality assessment
        self._print_info("Study quality assessment:")
        self._print_info("  HIGH: Well-validated assay, good controls, reproducible")
        self._print_info("  MODERATE: Reasonable assay but some limitations")
        self._print_info("  LOW: Preliminary data, poor controls, or unvalidated assay")
        
        quality = self.input.prompt_choice(
            "Overall quality of the functional evidence?",
            choices=['high', 'moderate', 'low'],
            default='moderate'
        )
        
        # Determine evidence code
        evidence = ManualEvidence()
        
        if damaging_count > 0 and benign_count > 0:
            # Conflicting evidence - be conservative
            self._print_info(f"Conflicting results: {damaging_count} damaging vs {benign_count} benign")
            if damaging_count > benign_count * 2:
                evidence.add_evidence(
                    'PS3_supporting',
                    f"Functional studies show predominantly damaging effect "
                    f"({damaging_count} damaging vs {benign_count} benign studies, {quality} quality), "
                    f"but some conflicting evidence exists."
                )
                self._print_success("Assigned: PS3_supporting (conflicting but predominantly damaging)")
            elif benign_count > damaging_count * 2:
                evidence.add_evidence(
                    'BS3_supporting',
                    f"Functional studies show predominantly no effect "
                    f"({benign_count} benign vs {damaging_count} damaging studies, {quality} quality), "
                    f"but some conflicting evidence exists."
                )
                self._print_success("Assigned: BS3_supporting (conflicting but predominantly benign)")
            else:
                self._print_info("Conflicting functional evidence with no clear direction → No PS3/BS3")
        
        elif damaging_count > 0:
            # Damaging effect - assign PS3
            code = self._map_ps3_strength(damaging_count, quality)
            if code:
                evidence.add_evidence(
                    code,
                    f"Functional studies ({damaging_count} independent studies, {quality} quality) "
                    f"demonstrate damaging effect on protein function."
                )
                self._print_success(f"Assigned: {code}")
        
        elif benign_count > 0:
            # Benign effect - assign BS3
            code = self._map_bs3_strength(benign_count, quality)
            if code:
                evidence.add_evidence(
                    code,
                    f"Functional studies ({benign_count} independent studies, {quality} quality) "
                    f"show no damaging effect on protein function."
                )
                self._print_success(f"Assigned: {code}")
        
        return evidence
    
    def _map_ps3_strength(self, study_count: int, quality: str) -> Optional[str]:
        """Map study count and quality to PS3 strength level."""
        if quality == 'high':
            if study_count >= 2:
                return 'PS3_strong'
            elif study_count >= 1:
                return 'PS3_moderate'
        elif quality == 'moderate':
            if study_count >= 2:
                return 'PS3_moderate'
            elif study_count >= 1:
                return 'PS3_supporting'
        elif quality == 'low':
            if study_count >= 2:
                return 'PS3_supporting'
        return None
    
    def _map_bs3_strength(self, study_count: int, quality: str) -> Optional[str]:
        """Map study count and quality to BS3 strength level."""
        if quality == 'high':
            if study_count >= 2:
                return 'BS3_strong'
            elif study_count >= 1:
                return 'BS3_moderate'
        elif quality == 'moderate':
            if study_count >= 2:
                return 'BS3_moderate'
            elif study_count >= 1:
                return 'BS3_supporting'
        elif quality == 'low':
            if study_count >= 2:
                return 'BS3_supporting'
        return None
    
    # =========================================================================
    # PS4: Case-Control / Enrichment Data
    # =========================================================================
    def collect_ps4(self) -> ManualEvidence:
        """
        Collect PS4 evidence based on case-control data.
        
        PS4: The prevalence of the variant in affected individuals is
             significantly increased compared with the prevalence in controls.
        
        Returns:
            ManualEvidence with PS4 code if applicable
        """
        self._print_header("PS4: Case-Control / Enrichment Evidence")
        self._print_info("PS4 applies when variant is significantly enriched in cases vs controls.")
        self._print_note("This requires published case-control or cohort study data.")
        
        has_data = self.input.prompt_yes_no(
            "Is there published case-control data for this variant?",
            default=False
        )
        
        if not has_data:
            self._print_info("No case-control data → No PS4 evidence")
            return ManualEvidence()
        
        # Collect case-control numbers
        cases_with = self.input.prompt_integer(
            "Number of AFFECTED individuals WITH the variant?",
            min_val=0, max_val=10000, default=0
        ) or 0
        
        total_cases = self.input.prompt_integer(
            "Total number of AFFECTED individuals studied?",
            min_val=1, max_val=100000, default=1
        ) or 1
        
        controls_with = self.input.prompt_integer(
            "Number of CONTROL individuals WITH the variant?",
            min_val=0, max_val=10000, default=0
        ) or 0
        
        total_controls = self.input.prompt_integer(
            "Total number of CONTROL individuals studied?",
            min_val=1, max_val=100000, default=1
        ) or 1
        
        # Calculate odds ratio
        evidence = ManualEvidence()
        
        if cases_with > 0 and total_controls > 0:
            # Calculate odds ratio with continuity correction
            a = cases_with
            b = total_cases - cases_with
            c = max(controls_with, 0.5)  # Continuity correction
            d = total_controls - controls_with
            
            if b > 0 and c > 0:
                odds_ratio = (a * d) / (b * c)
                
                self._print_info(f"Cases: {cases_with}/{total_cases} ({100*cases_with/total_cases:.1f}%)")
                self._print_info(f"Controls: {controls_with}/{total_controls} ({100*controls_with/total_controls:.2f}%)")
                self._print_info(f"Calculated odds ratio: {odds_ratio:.2f}")
                
                # Map to evidence strength
                code = self._map_ps4_strength(cases_with, odds_ratio)
                if code:
                    evidence.add_evidence(
                        code,
                        f"Case-control data shows variant enrichment: "
                        f"{cases_with}/{total_cases} cases vs {controls_with}/{total_controls} controls "
                        f"(OR={odds_ratio:.2f})."
                    )
                    self._print_success(f"Assigned: {code}")
                else:
                    self._print_info("Enrichment not significant enough for PS4")
            else:
                self._print_info("Cannot calculate odds ratio with given data")
        else:
            self._print_info("Insufficient data for PS4 assessment")
        
        return evidence
    
    def _map_ps4_strength(self, case_count: int, odds_ratio: float) -> Optional[str]:
        """Map case count and odds ratio to PS4 strength level."""
        if case_count >= 5 and odds_ratio >= 5.0:
            return 'PS4_strong'
        elif case_count >= 3 and odds_ratio >= 3.0:
            return 'PS4_moderate'
        elif case_count >= 2 and odds_ratio >= 2.0:
            return 'PS4_supporting'
        return None
    
    # =========================================================================
    # PP1 / BS4: Segregation Data
    # =========================================================================
    def collect_pp1_bs4(self) -> ManualEvidence:
        """
        Collect PP1/BS4 evidence based on segregation data.
        
        PP1: Co-segregation with disease in multiple affected family members
             in a gene definitively known to cause the disease.
        BS4: Lack of segregation in affected members of a family.
        
        Returns:
            ManualEvidence with PP1 or BS4 codes if applicable
        """
        self._print_header("PP1/BS4: Segregation Evidence")
        self._print_info("Segregation analysis examines whether the variant tracks")
        self._print_info("with disease status in family members.")
        self._print_note("Each affected carrier who inherited the variant = more evidence")
        
        has_family = self.input.prompt_yes_no(
            "Is there family segregation data available?",
            default=False
        )
        
        if not has_family:
            self._print_info("No family data → No PP1/BS4 evidence")
            return ManualEvidence()
        
        # Segregation direction
        self._print_info("Segregation patterns:")
        self._print_info("  COSEGREGATES: Variant tracks with disease in family")
        self._print_info("  DOES_NOT_SEGREGATE: Affected family members lack variant")
        self._print_info("  MIXED: Some affected have variant, some don't")
        
        pattern = self.input.prompt_choice(
            "What is the segregation pattern?",
            choices=['cosegregates', 'does_not_segregate', 'mixed'],
            default='cosegregates'
        )
        
        evidence = ManualEvidence()
        
        if pattern == 'cosegregates':
            # Count affected carriers (meioses)
            affected_carriers = self.input.prompt_integer(
                "How many AFFECTED family members carry the variant?",
                min_val=1, max_val=50, default=1
            ) or 1
            
            unaffected_noncarriers = self.input.prompt_integer(
                "How many UNAFFECTED family members lack the variant?",
                min_val=0, max_val=50, default=0
            ) or 0
            
            self._print_info(f"Affected carriers: {affected_carriers}")
            self._print_info(f"Unaffected non-carriers: {unaffected_noncarriers}")
            
            # Map to PP1 strength
            code = self._map_pp1_strength(affected_carriers)
            if code:
                evidence.add_evidence(
                    code,
                    f"Variant co-segregates with disease in {affected_carriers} affected "
                    f"family members. {unaffected_noncarriers} unaffected family members "
                    f"do not carry the variant."
                )
                self._print_success(f"Assigned: {code}")
            else:
                self._print_info("Insufficient segregation data for PP1")
        
        elif pattern == 'does_not_segregate':
            # Lack of segregation suggests benign
            affected_noncarriers = self.input.prompt_integer(
                "How many AFFECTED family members LACK the variant?",
                min_val=1, max_val=50, default=1
            ) or 1
            
            self._print_info(f"Affected non-carriers: {affected_noncarriers}")
            
            # Map to BS4 strength
            code = self._map_bs4_strength(affected_noncarriers)
            if code:
                evidence.add_evidence(
                    code,
                    f"Variant does not segregate with disease: {affected_noncarriers} "
                    f"affected family members do not carry the variant."
                )
                self._print_success(f"Assigned: {code}")
            else:
                self._print_info("Insufficient non-segregation data for BS4")
        
        else:  # mixed
            self._print_info("Mixed segregation pattern → No clear PP1/BS4 evidence")
            self._print_note("Consider incomplete penetrance or phenocopies")
        
        return evidence
    
    def _map_pp1_strength(self, affected_carriers: int) -> Optional[str]:
        """Map number of affected carriers to PP1 strength level."""
        if affected_carriers >= 5:
            return 'PP1_strong'
        elif affected_carriers >= 3:
            return 'PP1_moderate'
        elif affected_carriers >= 2:
            return 'PP1_supporting'
        return None
    
    def _map_bs4_strength(self, affected_noncarriers: int) -> Optional[str]:
        """Map number of affected non-carriers to BS4 strength level."""
        if affected_noncarriers >= 3:
            return 'BS4_moderate'
        elif affected_noncarriers >= 2:
            return 'BS4_supporting'
        return None
    
    # =========================================================================
    # PS1 / PM5: Same Codon / Same Amino Acid Change
    # =========================================================================
    def collect_ps1_pm5(self) -> ManualEvidence:
        """
        Collect PS1/PM5 evidence based on known pathogenic variants.
        
        PS1: Same amino acid change as a previously established pathogenic variant
             regardless of nucleotide change.
        PM5: Novel missense change at an amino acid residue where a different
             missense change has been determined to be pathogenic.
        
        Returns:
            ManualEvidence with PS1 or PM5 codes if applicable
        """
        self._print_header("PS1/PM5: Same Codon / Same Amino Acid Evidence")
        self._print_info("PS1: Exact same amino acid change as known pathogenic variant")
        self._print_info("PM5: Different amino acid change at same position as known pathogenic")
        self._print_note("Check ClinVar, HGMD, or literature for known pathogenic variants")
        
        # Check for same amino acid change (PS1)
        has_same_aa = self.input.prompt_yes_no(
            "Is there a PATHOGENIC variant with the SAME amino acid change\n"
            "(same ref AA → same alt AA, but possibly different nucleotide)?",
            default=False
        )
        
        evidence = ManualEvidence()
        
        if has_same_aa:
            self._print_info("Checking applicability...")
            
            # Check for splicing impact (PS1 caveat)
            splice_concern = self.input.prompt_yes_no(
                "Could the nucleotide change affect splicing differently\n"
                "than the known pathogenic variant?",
                default=False
            )
            
            if splice_concern:
                self._print_info("Splicing concern exists → PS1 not applicable")
                self._print_note("Consider PM5 if at same position with different AA change")
            else:
                evidence.add_evidence(
                    'PS1_strong',
                    "Same amino acid change as a known pathogenic variant "
                    "(different nucleotide change). No splicing concern identified."
                )
                self._print_success("Assigned: PS1_strong")
                return evidence  # PS1 takes precedence over PM5
        
        # Check for different AA change at same position (PM5)
        has_same_position = self.input.prompt_yes_no(
            "Is there a PATHOGENIC variant at the SAME amino acid position\n"
            "but with a DIFFERENT amino acid change?",
            default=False
        )
        
        if has_same_position:
            # Check how pathogenic the known variant is
            known_classification = self.input.prompt_choice(
                "Classification of the known pathogenic variant at same position?",
                choices=['pathogenic', 'likely_pathogenic'],
                default='pathogenic'
            )
            
            if known_classification == 'pathogenic':
                evidence.add_evidence(
                    'PM5_moderate',
                    "Different missense change at same amino acid position as a known "
                    "Pathogenic variant. Suggests this position is critical for protein function."
                )
                self._print_success("Assigned: PM5_moderate")
            else:
                evidence.add_evidence(
                    'PM5_supporting',
                    "Different missense change at same amino acid position as a Likely Pathogenic "
                    "variant. Position may be important for protein function."
                )
                self._print_success("Assigned: PM5_supporting")
        else:
            self._print_info("No known pathogenic variants at this position → No PS1/PM5")
        
        return evidence
    
    # =========================================================================
    # PP5 / BP6: Reputable Source
    # =========================================================================
    def collect_pp5_bp6(self) -> ManualEvidence:
        """
        Collect PP5/BP6 evidence based on reputable clinical assertions.
        
        PP5: Reputable source recently reports variant as pathogenic, but the
             evidence is not available to the laboratory.
        BP6: Reputable source recently reports variant as benign, but the
             evidence is not available to the laboratory.
        
        Note: These criteria are deprecated in some frameworks (e.g., ClinGen)
        because they can lead to circular evidence. Use cautiously.
        
        Returns:
            ManualEvidence with PP5 or BP6 codes if applicable
        """
        self._print_header("PP5/BP6: Reputable Source Evidence")
        self._print_info("PP5/BP6 apply when expert labs have classified this variant")
        self._print_info("but their full evidence is not available for review.")
        self._print_note("⚠️  These criteria are DEPRECATED by some guidelines (ClinGen).")
        self._print_note("Use only when independent evidence evaluation is not possible.")
        
        use_criteria = self.input.prompt_yes_no(
            "Do you want to apply PP5/BP6 criteria?\n"
            "(Not recommended if you can evaluate the underlying evidence directly)",
            default=False
        )
        
        if not use_criteria:
            self._print_info("PP5/BP6 not applied (recommended)")
            return ManualEvidence()
        
        evidence = ManualEvidence()
        
        # Check for pathogenic assertion
        has_path_assertion = self.input.prompt_yes_no(
            "Has a reputable source (e.g., expert lab, ClinVar 2+ stars)\n"
            "classified this variant as PATHOGENIC or LIKELY PATHOGENIC?",
            default=False
        )
        
        if has_path_assertion:
            # Check for conflicts
            has_conflicts = self.input.prompt_yes_no(
                "Are there conflicting interpretations (some benign, some pathogenic)?",
                default=False
            )
            
            if has_conflicts:
                self._print_info("Conflicting interpretations → PP5 not applicable")
            else:
                star_level = self.input.prompt_choice(
                    "ClinVar review status?",
                    choices=['4_star', '3_star', '2_star', '1_star', 'unknown'],
                    default='2_star'
                )
                
                if star_level in ['4_star', '3_star', '2_star']:
                    evidence.add_evidence(
                        'PP5_supporting',
                        f"Reputable source classifies variant as pathogenic "
                        f"(ClinVar {star_level.replace('_', ' ')})."
                    )
                    self._print_success("Assigned: PP5_supporting")
                else:
                    self._print_info("Low review status → PP5 not applied")
        
        # Check for benign assertion (only if no pathogenic)
        if not evidence.has_evidence():
            has_benign_assertion = self.input.prompt_yes_no(
                "Has a reputable source (e.g., expert lab, ClinVar 2+ stars)\n"
                "classified this variant as BENIGN or LIKELY BENIGN?",
                default=False
            )
            
            if has_benign_assertion:
                has_conflicts = self.input.prompt_yes_no(
                    "Are there conflicting interpretations?",
                    default=False
                )
                
                if has_conflicts:
                    self._print_info("Conflicting interpretations → BP6 not applicable")
                else:
                    evidence.add_evidence(
                        'BP6_supporting',
                        "Reputable source classifies variant as benign."
                    )
                    self._print_success("Assigned: BP6_supporting")
        
        if not evidence.has_evidence():
            self._print_info("No reputable source assertion → No PP5/BP6")
        
        return evidence
    
    # =========================================================================
    # Main Collection Methods
    # =========================================================================
    def collect_all(self, variant_data: Any = None) -> ManualEvidence:
        """
        Collect all interactive evidence through a guided interview.
        
        Args:
            variant_data: Optional VariantData object for context
            
        Returns:
            ManualEvidence with all collected evidence codes
        """
        self._print_header("INTERACTIVE LITERATURE-BASED EVIDENCE COLLECTION")
        self._print_info("You will be asked about literature-based evidence that requires")
        self._print_info("human interpretation: functional studies, segregation, case-control,")
        self._print_info("known pathogenic variants at same position, and expert assertions.")
        self._print_note("Answer conservatively: if unsure, choose 'no' or skip.")
        print()
        
        # Collect each evidence type
        ps3_bs3 = self.collect_ps3_bs3()
        ps4 = self.collect_ps4()
        pp1_bs4 = self.collect_pp1_bs4()
        ps1_pm5 = self.collect_ps1_pm5()
        pp5_bp6 = self.collect_pp5_bp6()
        
        # Merge all evidence
        self.evidence = ManualEvidence()
        for ev in [ps3_bs3, ps4, pp1_bs4, ps1_pm5, pp5_bp6]:
            self.evidence = self.evidence.merge_with(ev)
        
        # Summary
        self._print_header("MANUAL EVIDENCE SUMMARY")
        if self.evidence.has_evidence():
            self._print_info(f"Collected {len(self.evidence.codes)} evidence code(s):")
            for code in self.evidence.codes:
                print(f"  • {code}")
                if code in self.evidence.explanations:
                    print(f"    → {self.evidence.explanations[code][:80]}...")
        else:
            self._print_info("No manual evidence collected.")
        
        return self.evidence
    
    def collect_selective(self, criteria: List[str], 
                          variant_data: Any = None) -> ManualEvidence:
        """
        Collect only selected criteria.
        
        Args:
            criteria: List of criteria to collect (e.g., ['PS3_BS3', 'PP1_BS4'])
            variant_data: Optional VariantData object for context
            
        Returns:
            ManualEvidence with collected evidence codes
        """
        self.evidence = ManualEvidence()
        
        criteria_map = {
            'PS3_BS3': self.collect_ps3_bs3,
            'PS4': self.collect_ps4,
            'PP1_BS4': self.collect_pp1_bs4,
            'PS1_PM5': self.collect_ps1_pm5,
            'PP5_BP6': self.collect_pp5_bp6,
        }
        
        for criterion in criteria:
            if criterion in criteria_map:
                ev = criteria_map[criterion]()
                self.evidence = self.evidence.merge_with(ev)
        
        return self.evidence


# =============================================================================
# Pure Functions for Evidence Mapping (for unit testing)
# =============================================================================
def map_functional_studies_to_evidence(
    damaging_count: int,
    benign_count: int,
    quality: str
) -> Tuple[Optional[str], str]:
    """
    Pure function to map functional study data to evidence code.
    
    Args:
        damaging_count: Number of studies showing damaging effect
        benign_count: Number of studies showing no effect
        quality: Study quality ('high', 'moderate', 'low')
        
    Returns:
        Tuple of (evidence_code or None, explanation)
    """
    if damaging_count == 0 and benign_count == 0:
        return None, "No functional studies available"
    
    if damaging_count > 0 and benign_count > 0:
        # Conflicting - be conservative
        if damaging_count > benign_count * 2:
            return 'PS3_supporting', f"Conflicting but predominantly damaging ({damaging_count} vs {benign_count})"
        elif benign_count > damaging_count * 2:
            return 'BS3_supporting', f"Conflicting but predominantly benign ({benign_count} vs {damaging_count})"
        else:
            return None, "Conflicting functional evidence with no clear direction"
    
    if damaging_count > 0:
        if quality == 'high':
            if damaging_count >= 2:
                return 'PS3_strong', f"{damaging_count} high-quality damaging studies"
            else:
                return 'PS3_moderate', f"{damaging_count} high-quality damaging study"
        elif quality == 'moderate':
            if damaging_count >= 2:
                return 'PS3_moderate', f"{damaging_count} moderate-quality damaging studies"
            else:
                return 'PS3_supporting', f"{damaging_count} moderate-quality damaging study"
        else:
            if damaging_count >= 2:
                return 'PS3_supporting', f"{damaging_count} low-quality damaging studies"
            else:
                return None, "Single low-quality study insufficient"
    
    if benign_count > 0:
        if quality == 'high':
            if benign_count >= 2:
                return 'BS3_strong', f"{benign_count} high-quality benign studies"
            else:
                return 'BS3_moderate', f"{benign_count} high-quality benign study"
        elif quality == 'moderate':
            if benign_count >= 2:
                return 'BS3_moderate', f"{benign_count} moderate-quality benign studies"
            else:
                return 'BS3_supporting', f"{benign_count} moderate-quality benign study"
        else:
            if benign_count >= 2:
                return 'BS3_supporting', f"{benign_count} low-quality benign studies"
            else:
                return None, "Single low-quality study insufficient"
    
    return None, "No evidence criteria met"


def map_segregation_to_evidence(
    pattern: str,
    affected_carriers: int = 0,
    affected_noncarriers: int = 0
) -> Tuple[Optional[str], str]:
    """
    Pure function to map segregation data to evidence code.
    
    Args:
        pattern: 'cosegregates', 'does_not_segregate', or 'mixed'
        affected_carriers: Number of affected individuals with variant
        affected_noncarriers: Number of affected individuals without variant
        
    Returns:
        Tuple of (evidence_code or None, explanation)
    """
    if pattern == 'cosegregates':
        if affected_carriers >= 5:
            return 'PP1_strong', f"Cosegregation in {affected_carriers} affected carriers"
        elif affected_carriers >= 3:
            return 'PP1_moderate', f"Cosegregation in {affected_carriers} affected carriers"
        elif affected_carriers >= 2:
            return 'PP1_supporting', f"Cosegregation in {affected_carriers} affected carriers"
        else:
            return None, "Insufficient segregation data (need ≥2 affected carriers)"
    
    elif pattern == 'does_not_segregate':
        if affected_noncarriers >= 3:
            return 'BS4_moderate', f"Non-segregation: {affected_noncarriers} affected lack variant"
        elif affected_noncarriers >= 2:
            return 'BS4_supporting', f"Non-segregation: {affected_noncarriers} affected lack variant"
        else:
            return None, "Insufficient non-segregation data"
    
    else:  # mixed
        return None, "Mixed segregation pattern - no clear evidence"


def map_case_control_to_evidence(
    cases_with: int,
    total_cases: int,
    controls_with: int,
    total_controls: int
) -> Tuple[Optional[str], str]:
    """
    Pure function to map case-control data to evidence code.
    
    Args:
        cases_with: Number of cases with variant
        total_cases: Total cases studied
        controls_with: Number of controls with variant
        total_controls: Total controls studied
        
    Returns:
        Tuple of (evidence_code or None, explanation)
    """
    if cases_with == 0 or total_cases == 0 or total_controls == 0:
        return None, "Insufficient case-control data"
    
    # Calculate odds ratio with continuity correction
    a = cases_with
    b = total_cases - cases_with
    c = max(controls_with, 0.5)
    d = total_controls - controls_with
    
    if b <= 0 or c <= 0:
        return None, "Cannot calculate odds ratio"
    
    odds_ratio = (a * d) / (b * c)
    
    if cases_with >= 5 and odds_ratio >= 5.0:
        return 'PS4_strong', f"Strong enrichment in cases (OR={odds_ratio:.2f}, n={cases_with})"
    elif cases_with >= 3 and odds_ratio >= 3.0:
        return 'PS4_moderate', f"Moderate enrichment in cases (OR={odds_ratio:.2f}, n={cases_with})"
    elif cases_with >= 2 and odds_ratio >= 2.0:
        return 'PS4_supporting', f"Supporting enrichment in cases (OR={odds_ratio:.2f}, n={cases_with})"
    else:
        return None, f"Enrichment not significant (OR={odds_ratio:.2f}, n={cases_with})"

"""
Phenotype Matcher Module
========================

This module provides local, offline phenotype-matching capabilities for PP4/BP5-style
ACMG evidence evaluation. It uses local JSON data files for gene-phenotype associations
and HPO term normalization.

IMPORTANT: This is an EDUCATIONAL APPROXIMATION of phenotype matching, not a
production-grade HPO semantic engine. For clinical use, consider integrating
with proper HPO ontology tools or semantic similarity APIs.

Classes:
    - HPOClient: Local HPO term normalization layer
    - GenePhenotypeDatabase: Local gene-phenotype association database
    - PhenotypeSimilarityCalculator: Weighted Jaccard similarity scoring
    - PhenotypeMatcher: Main interface for phenotype-based evidence evaluation
"""

import json
import os
import re
from pathlib import Path
from typing import Dict, List, Optional, Set, Any, Union

# Import thresholds and low-information HPOs from config
try:
    from config.constants import (
        PHENOTYPE_SIMILARITY_THRESHOLDS,
        LOW_INFORMATION_HPO,
        LOW_INFORMATION_HPO_WEIGHT,
        MIN_TERMS_FOR_PHENOTYPE_EVIDENCE
    )
except ImportError:
    # Fallback defaults if config import fails
    PHENOTYPE_SIMILARITY_THRESHOLDS = {
        'PP4': 0.8,
        'PP4_SUPPORTING': 0.5,
        'BP5': 0.2,
    }
    LOW_INFORMATION_HPO = {
        "HP:0002664",  # Neoplasm
        "HP:0000118",  # Phenotypic abnormality
    }
    LOW_INFORMATION_HPO_WEIGHT = 0.3
    MIN_TERMS_FOR_PHENOTYPE_EVIDENCE = 3


class HPOClient:
    """
    Local HPO term normalization layer.
    
    This class provides offline normalization of phenotype terms to HPO IDs.
    It uses a local synonyms dictionary to map free-text phenotype descriptions
    to standardized HPO identifiers.
    
    NOTE: This is a simple text-matching approach. For production use, consider
    using proper HPO ontology tools with semantic similarity.
    
    Example:
        >>> client = HPOClient()
        >>> client.get_phenotype_terms("breast cancer")
        {'HP:0003002'}
        >>> client.get_phenotype_terms("HP:0001250")
        {'HP:0001250'}
    """
    
    # Regex pattern for valid HPO IDs (e.g., HP:0000001)
    HPO_PATTERN = re.compile(r'^HP:\d{7}$')
    
    def __init__(self, synonyms_path: Optional[str] = None):
        """
        Initialize the HPO client with local synonyms data.
        
        Args:
            synonyms_path: Optional path to HPO synonyms JSON file.
                          Defaults to src/data/hpo_synonyms.json
        """
        self._synonyms_cache: Optional[Dict] = None
        self._synonyms_path = synonyms_path
        
    def _get_default_path(self) -> Path:
        """Get the default path to the HPO synonyms file."""
        # Try multiple locations to find the data file
        possible_paths = [
            Path(__file__).parent.parent / 'data' / 'hpo_synonyms.json',
            Path(__file__).parent / 'data' / 'hpo_synonyms.json',
            Path('src/data/hpo_synonyms.json'),
            Path('data/hpo_synonyms.json'),
        ]
        
        for path in possible_paths:
            if path.exists():
                return path
        
        # Return the first path even if it doesn't exist (will create empty cache)
        return possible_paths[0]
    
    def _load_synonyms(self) -> Dict:
        """Load the synonyms dictionary from JSON file."""
        if self._synonyms_cache is not None:
            return self._synonyms_cache
        
        path = Path(self._synonyms_path) if self._synonyms_path else self._get_default_path()
        
        try:
            with open(path, 'r', encoding='utf-8') as f:
                data = json.load(f)
                self._synonyms_cache = data.get('text_to_hpo', {})
        except (FileNotFoundError, json.JSONDecodeError) as e:
            # Silently fall back to empty cache in production
            self._synonyms_cache = {}
        
        return self._synonyms_cache
    
    def is_hpo_id(self, term: str) -> bool:
        """
        Check if a term is a valid HPO ID.
        
        Args:
            term: The term to check
            
        Returns:
            True if the term matches HPO ID pattern (HP:NNNNNNN)
        """
        if not term:
            return False
        return bool(self.HPO_PATTERN.match(term.strip().upper()))
    
    def normalize_term(self, term: str) -> str:
        """
        Normalize a phenotype term for matching.
        
        Applies the following transformations:
        1. Convert to lowercase
        2. Strip leading/trailing whitespace
        3. Collapse multiple spaces to single space
        
        Args:
            term: The phenotype term to normalize
            
        Returns:
            Normalized term string
        """
        if not term:
            return ""
        
        # Lowercase and strip
        normalized = term.strip().lower()
        
        # Remove extra whitespace
        normalized = ' '.join(normalized.split())
        
        return normalized
    
    def _try_singular_form(self, term: str) -> Optional[str]:
        """
        Attempt to convert plural term to singular for synonym lookup.
        
        Simple heuristic: if term ends with 's' or 'es', try without suffix.
        This helps match terms like "seizures" -> "seizure".
        
        Args:
            term: Normalized term to try converting
            
        Returns:
            Singular form if different from input, None otherwise
        """
        if not term or len(term) < 3:
            return None
        
        # Try removing 'es' suffix first (e.g., "seizures" -> "seizur" won't work, 
        # but "masses" -> "mass" might)
        if term.endswith('es') and len(term) > 3:
            singular = term[:-2]
            if singular != term:
                return singular
        
        # Try removing 's' suffix (e.g., "seizures" -> "seizure")
        if term.endswith('s') and not term.endswith('ss'):
            singular = term[:-1]
            if singular != term:
                return singular
        
        return None
    
    def get_phenotype_terms(self, phenotype: Union[str, List[str]]) -> Set[str]:
        """
        Convert phenotype input to a set of normalized HPO IDs/tokens.
        
        This method handles:
        - Already valid HPO IDs: returned as-is
        - Free text descriptions: mapped via local synonyms dictionary
        - Plural forms: attempts singular fallback if plural lookup fails
        - Unknown terms: wrapped as normalized tokens for set comparison
        
        Args:
            phenotype: A single phenotype term or list of terms.
                      Can be HPO IDs (HP:NNNNNNN) or text descriptions.
                      
        Returns:
            Set of normalized HPO IDs or tokens.
            
        Example:
            >>> client = HPOClient()
            >>> client.get_phenotype_terms("HP:0001250")
            {'HP:0001250'}
            >>> client.get_phenotype_terms(["breast cancer", "HP:0000137"])
            {'HP:0003002', 'HP:0000137'}
            >>> client.get_phenotype_terms("Seizures")  # uppercase + plural
            {'HP:0001250'}
        """
        if not phenotype:
            return set()
        
        # Handle single term vs list
        if isinstance(phenotype, str):
            terms = [phenotype]
        else:
            terms = list(phenotype)
        
        result = set()
        synonyms = self._load_synonyms()
        
        for term in terms:
            if not term:
                continue
                
            term = str(term).strip()
            
            # Case 1: Already a valid HPO ID
            if self.is_hpo_id(term):
                result.add(term.upper())
                continue
            
            # Case 2: Try to map via synonyms (exact match after normalization)
            normalized = self.normalize_term(term)
            if normalized in synonyms:
                result.add(synonyms[normalized])
                continue
            
            # Case 3: Try singular form if plural lookup failed
            # (e.g., "seizures" -> "seizure", "Carcinomas" -> "carcinoma")
            singular = self._try_singular_form(normalized)
            if singular and singular in synonyms:
                result.add(synonyms[singular])
                continue
            
            # Case 4: Try partial matching (for compound terms)
            matched = False
            for syn_text, hpo_id in synonyms.items():
                if syn_text in normalized or normalized in syn_text:
                    result.add(hpo_id)
                    matched = True
                    break
            
            if matched:
                continue
            
            # Case 5: Unknown term - wrap as token for set comparison
            # This allows the similarity calculator to still work with unknown terms
            result.add(f"TEXT:{normalized}")
        
        return result


class GenePhenotypeDatabase:
    """
    Local gene-phenotype association database.
    
    This class loads gene-phenotype mappings from a local JSON file and provides
    lookup functionality. Data is cached in memory after first load.
    
    The database maps gene symbols to sets of HPO IDs representing the
    phenotypes associated with that gene in disease.
    
    Example:
        >>> db = GenePhenotypeDatabase()
        >>> db.get_gene_phenotypes("BRCA1")
        {'HP:0003002', 'HP:0100013', 'HP:0010619', 'HP:0000137', 'HP:0030075'}
        >>> db.get_gene_phenotypes("UNKNOWN_GENE")
        set()
    """
    
    def __init__(self, data_path: Optional[str] = None):
        """
        Initialize the database with local gene-phenotype data.
        
        Args:
            data_path: Optional path to gene phenotypes JSON file.
                      Defaults to src/data/gene_phenotypes.json
        """
        self._cache: Optional[Dict] = None
        self._data_path = data_path
    
    def _get_default_path(self) -> Path:
        """Get the default path to the gene phenotypes file."""
        possible_paths = [
            Path(__file__).parent.parent / 'data' / 'gene_phenotypes.json',
            Path(__file__).parent / 'data' / 'gene_phenotypes.json',
            Path('src/data/gene_phenotypes.json'),
            Path('data/gene_phenotypes.json'),
        ]
        
        for path in possible_paths:
            if path.exists():
                return path
        
        return possible_paths[0]
    
    def _load_data(self) -> Dict:
        """Load the gene-phenotype data from JSON file."""
        if self._cache is not None:
            return self._cache
        
        path = Path(self._data_path) if self._data_path else self._get_default_path()
        
        try:
            with open(path, 'r', encoding='utf-8') as f:
                data = json.load(f)
                # Remove metadata keys (those starting with _)
                self._cache = {
                    k: v for k, v in data.items() 
                    if not k.startswith('_')
                }
        except (FileNotFoundError, json.JSONDecodeError) as e:
            # Silently fall back to empty cache
            self._cache = {}
        
        return self._cache
    
    def get_gene_phenotypes(self, gene: str) -> Set[str]:
        """
        Get the set of HPO IDs associated with a gene.
        
        Args:
            gene: Gene symbol (case-insensitive, e.g., 'BRCA1', 'brca1')
            
        Returns:
            Set of HPO IDs associated with the gene.
            Returns empty set if gene is unknown (no exception raised).
            
        Example:
            >>> db = GenePhenotypeDatabase()
            >>> db.get_gene_phenotypes("TP53")
            {'HP:0002664', 'HP:0100526', 'HP:0003002', 'HP:0030416', ...}
        """
        if not gene:
            return set()
        
        data = self._load_data()
        
        # Normalize gene symbol (uppercase)
        gene_upper = gene.strip().upper()
        
        # Lookup in database
        gene_data = data.get(gene_upper)
        
        if gene_data is None:
            return set()
        
        # Extract HPO terms from gene data
        if isinstance(gene_data, dict):
            hpo_terms = gene_data.get('hpo_terms', [])
        elif isinstance(gene_data, list):
            hpo_terms = gene_data
        else:
            hpo_terms = []
        
        return set(hpo_terms)
    
    def get_gene_info(self, gene: str) -> Optional[Dict[str, Any]]:
        """
        Get full information about a gene's phenotype associations.
        
        Args:
            gene: Gene symbol
            
        Returns:
            Dictionary with hpo_terms, phenotype_names, inheritance, disease,
            or None if gene is unknown.
        """
        if not gene:
            return None
        
        data = self._load_data()
        gene_upper = gene.strip().upper()
        
        return data.get(gene_upper)
    
    def has_gene(self, gene: str) -> bool:
        """Check if a gene exists in the database."""
        if not gene:
            return False
        
        data = self._load_data()
        return gene.strip().upper() in data
    
    def get_all_genes(self) -> List[str]:
        """Get list of all genes in the database."""
        data = self._load_data()
        return list(data.keys())


class PhenotypeSimilarityCalculator:
    """
    Weighted Jaccard phenotype similarity calculator.
    
    This class computes weighted similarity between two sets of phenotype terms.
    Low-information HPO terms (e.g., "Neoplasm", "Phenotypic abnormality") are
    down-weighted to prevent false PP4 evidence from overly generic terms.
    
    WHY WEIGHTED: Terms like HP:0002664 (Neoplasm) match many genes non-specifically.
    If a patient presents only with "cancer", standard Jaccard would give high
    similarity to any cancer-associated gene, causing false PP4. Down-weighting
    generic terms ensures that specific phenotypes drive the similarity score.
    
    Weighted Jaccard formula:
        - Intersection weight = sum of weights for common HPOs
        - Union weight = sum of weights for all unique HPOs
        - Similarity = intersection_weight / union_weight
    
    Example:
        >>> calc = PhenotypeSimilarityCalculator()
        >>> # Generic "Neoplasm" alone should not trigger high similarity
        >>> calc.calculate_similarity({'HP:0002664'}, {'HP:0002664', 'HP:0003002'})
        0.2...  # Low because HP:0002664 is down-weighted
    """
    
    def __init__(
        self,
        low_info_hpos: Optional[Set[str]] = None,
        low_info_weight: float = LOW_INFORMATION_HPO_WEIGHT
    ):
        """
        Initialize the calculator with low-information HPO configuration.
        
        Args:
            low_info_hpos: Set of HPO IDs to down-weight. Defaults to LOW_INFORMATION_HPO.
            low_info_weight: Weight for low-information terms (0.0-1.0). Default 0.3.
        """
        self.low_info_hpos = low_info_hpos or LOW_INFORMATION_HPO
        self.low_info_weight = low_info_weight
    
    def _get_term_weight(self, term: str) -> float:
        """
        Get the weight for a phenotype term.
        
        Normal terms get weight 1.0, low-information terms get reduced weight.
        TEXT: prefixed terms (unknown phenotypes) also get reduced weight.
        
        Args:
            term: HPO ID or TEXT: prefixed token
            
        Returns:
            Weight between 0.0 and 1.0
        """
        # Unknown text terms get reduced weight
        if term.startswith("TEXT:"):
            return self.low_info_weight
        
        # Low-information HPO terms get reduced weight
        if term in self.low_info_hpos:
            return self.low_info_weight
        
        # Normal terms get full weight
        return 1.0
    
    def calculate_similarity(
        self, 
        patient_terms: Union[Set[str], List[str]], 
        gene_terms: Union[Set[str], List[str]]
    ) -> float:
        """
        Calculate weighted Jaccard similarity between two sets of phenotype terms.
        
        Weighted Jaccard = sum(weights of intersection) / sum(weights of union)
        
        Low-information HPO terms (e.g., HP:0002664 "Neoplasm") receive reduced
        weight to prevent false PP4 from generic cancer-related terms alone.
        
        Args:
            patient_terms: Set or list of patient HPO terms/tokens
            gene_terms: Set or list of gene-associated HPO terms/tokens
            
        Returns:
            Float between 0.0 (no overlap) and 1.0 (identical sets with full weights).
            Returns 0.0 if either set is empty.
            
        Example:
            >>> calc = PhenotypeSimilarityCalculator()
            >>> # Specific terms give higher similarity
            >>> calc.calculate_similarity(
            ...     ['HP:0003002', 'HP:0000137'],  # Breast carcinoma, Ovarian neoplasm
            ...     ['HP:0003002', 'HP:0100013', 'HP:0000137']
            ... )
            0.666...
            >>> # Generic "Neoplasm" alone gives low similarity
            >>> calc.calculate_similarity(
            ...     ['HP:0002664'],  # Just "Neoplasm" (low-info)
            ...     ['HP:0002664', 'HP:0003002']
            ... )
            0.23...  # Low due to down-weighting
        """
        # Convert to sets if needed
        set_a = set(patient_terms) if not isinstance(patient_terms, set) else patient_terms
        set_b = set(gene_terms) if not isinstance(gene_terms, set) else gene_terms
        
        # Handle empty sets
        if not set_a or not set_b:
            return 0.0
        
        # Calculate weighted intersection and union
        intersection = set_a & set_b
        union = set_a | set_b
        
        if not union:
            return 0.0
        
        # Sum weights for intersection
        intersection_weight = sum(self._get_term_weight(term) for term in intersection)
        
        # Sum weights for union
        union_weight = sum(self._get_term_weight(term) for term in union)
        
        if union_weight == 0:
            return 0.0
        
        # Return weighted Jaccard, clamped to [0, 1]
        similarity = intersection_weight / union_weight
        return max(0.0, min(1.0, similarity))
    
    def calculate_overlap_ratio(
        self,
        patient_terms: Union[Set[str], List[str]],
        gene_terms: Union[Set[str], List[str]]
    ) -> float:
        """
        Calculate what fraction of patient terms match gene terms (weighted).
        
        This is an asymmetric measure using weights.
        Useful when patient phenotypes are the reference set.
        
        Args:
            patient_terms: Set of patient HPO terms
            gene_terms: Set of gene-associated HPO terms
            
        Returns:
            Float between 0.0 and 1.0 representing weighted fraction of patient
            terms that match gene phenotypes.
        """
        set_a = set(patient_terms) if not isinstance(patient_terms, set) else patient_terms
        set_b = set(gene_terms) if not isinstance(gene_terms, set) else gene_terms
        
        if not set_a:
            return 0.0
        
        intersection = set_a & set_b
        
        intersection_weight = sum(self._get_term_weight(term) for term in intersection)
        patient_weight = sum(self._get_term_weight(term) for term in set_a)
        
        if patient_weight == 0:
            return 0.0
        
        return intersection_weight / patient_weight


class PhenotypeMatcher:
    """
    Main interface for phenotype-based ACMG evidence evaluation.
    
    This class combines HPO normalization, gene-phenotype lookup, and weighted
    similarity calculation to provide PP4/BP5-style evidence based on phenotype matching.
    
    Thresholds are loaded from config.constants.PHENOTYPE_SIMILARITY_THRESHOLDS.
    Low-information HPO terms are down-weighted to prevent false PP4 from generic terms.
    
    Usage:
        >>> matcher = PhenotypeMatcher()
        >>> result = matcher.evaluate_phenotype_match(
        ...     variant_data,  # VariantData with gene info
        ...     ['HP:0003002', 'HP:0000137']  # Patient phenotypes
        ... )
        >>> print(result)
        {
            'similarity': 0.6,
            'evidence_code': 'PP4_supporting',
            'gene': 'BRCA1',
            'disease': 'Hereditary breast and ovarian cancer syndrome',
            'inheritance': 'autosomal_dominant',
            'explanation': 'Phenotype similarity 0.60 with BRCA1-associated phenotype set...'
        }
    """
    
    def __init__(
        self,
        gene_db_path: Optional[str] = None,
        synonyms_path: Optional[str] = None,
        thresholds: Optional[Dict[str, float]] = None
    ):
        """
        Initialize the phenotype matcher with its components.
        
        Args:
            gene_db_path: Optional path to gene phenotypes JSON file
            synonyms_path: Optional path to HPO synonyms JSON file
            thresholds: Optional custom thresholds for PP4/BP5 assignment.
                       Uses PHENOTYPE_SIMILARITY_THRESHOLDS from config if not provided.
        """
        self.hpo_client = HPOClient(synonyms_path)
        self.gene_phenotype_db = GenePhenotypeDatabase(gene_db_path)
        self.similarity_calculator = PhenotypeSimilarityCalculator()
        
        # Use thresholds from config, allowing override
        self.thresholds = thresholds or {
            'PP4': PHENOTYPE_SIMILARITY_THRESHOLDS.get('PP4', 0.8),
            'PP4_SUPPORTING': PHENOTYPE_SIMILARITY_THRESHOLDS.get('PP4_SUPPORTING', 0.5),
            'BP5': PHENOTYPE_SIMILARITY_THRESHOLDS.get('BP5', 0.2),
        }
        
        # Minimum terms for reliable evidence
        self.min_terms = MIN_TERMS_FOR_PHENOTYPE_EVIDENCE
    
    def evaluate_phenotype_match(
        self,
        variant_data,
        patient_phenotypes: Optional[Union[List[str], Set[str]]] = None
    ) -> Dict[str, Any]:
        """
        Evaluate phenotype-genotype correlation for PP4/BP5 evidence.
        
        This method:
        1. Gets the gene from variant_data
        2. Looks up gene-associated phenotypes from local database
        3. Normalizes patient phenotypes via HPO client
        4. Calculates weighted similarity score
        5. Returns evidence code based on configurable thresholds
        
        Args:
            variant_data: VariantData object containing at least basic_info.gene
            patient_phenotypes: List of patient HPO terms or text descriptions.
                              If None, tries to get from variant_data.patient_phenotypes
                              
        Returns:
            Dictionary with:
                - similarity: float (0.0 to 1.0)
                - evidence_code: str or None ('PP4', 'PP4_supporting', 'BP5', or None)
                - strength: str ('Supporting' or None)
                - gene: str (gene symbol)
                - disease: str (associated disease name, if available)
                - inheritance: str (inheritance pattern, if available)
                - explanation: str (rich human-readable explanation)
                - patient_terms: set of normalized patient HPO terms
                - gene_terms: set of gene-associated HPO terms
                
        Example:
            >>> matcher = PhenotypeMatcher()
            >>> result = matcher.evaluate_phenotype_match(
            ...     variant_data,  # BRCA1 variant
            ...     ['breast cancer', 'ovarian neoplasm']
            ... )
            >>> result['evidence_code']
            'PP4'  # High match for BRCA1
        """
        # Get gene from variant data
        gene = None
        if hasattr(variant_data, 'gene'):
            gene = variant_data.gene
        elif hasattr(variant_data, 'basic_info') and variant_data.basic_info:
            gene = variant_data.basic_info.get('gene')
        
        if not gene:
            return {
                'similarity': 0.0,
                'evidence_code': None,
                'strength': None,
                'gene': None,
                'disease': None,
                'inheritance': None,
                'explanation': 'No gene information available for phenotype matching.',
                'patient_terms': set(),
                'gene_terms': set()
            }
        
        # Get patient phenotypes from argument or variant_data
        if patient_phenotypes is None:
            patient_phenotypes = getattr(variant_data, 'patient_phenotypes', None)
        
        # Get gene info for disease/inheritance context
        gene_info = self.gene_phenotype_db.get_gene_info(gene)
        disease = gene_info.get('disease') if gene_info else None
        inheritance = gene_info.get('inheritance') if gene_info else None
        
        # Get gene-associated phenotypes
        gene_terms = self.gene_phenotype_db.get_gene_phenotypes(gene)
        
        if not gene_terms:
            return {
                'similarity': 0.0,
                'evidence_code': None,
                'strength': None,
                'gene': gene,
                'disease': None,
                'inheritance': None,
                'explanation': f'No curated phenotype data available for {gene}; phenotype-based evidence not applied.',
                'patient_terms': set(),
                'gene_terms': set()
            }
        
        # Handle missing patient phenotypes
        if not patient_phenotypes:
            return {
                'similarity': 0.0,
                'evidence_code': None,
                'strength': None,
                'gene': gene,
                'disease': disease,
                'inheritance': inheritance,
                'explanation': 'No patient phenotype data provided.',
                'patient_terms': set(),
                'gene_terms': gene_terms
            }
        
        # Normalize patient phenotypes
        patient_terms = self.hpo_client.get_phenotype_terms(patient_phenotypes)
        
        if not patient_terms:
            return {
                'similarity': 0.0,
                'evidence_code': None,
                'strength': None,
                'gene': gene,
                'disease': disease,
                'inheritance': inheritance,
                'explanation': 'Could not normalize patient phenotype terms.',
                'patient_terms': set(),
                'gene_terms': gene_terms
            }
        
        # Calculate weighted similarity
        similarity = self.similarity_calculator.calculate_similarity(
            patient_terms, gene_terms
        )
        
        # Determine evidence code based on thresholds
        evidence_code, strength, explanation = self._determine_evidence(
            similarity, gene, disease, inheritance, patient_terms, gene_terms
        )
        
        return {
            'similarity': similarity,
            'evidence_code': evidence_code,
            'strength': strength,
            'gene': gene,
            'disease': disease,
            'inheritance': inheritance,
            'explanation': explanation,
            'patient_terms': patient_terms,
            'gene_terms': gene_terms
        }
    
    def _determine_evidence(
        self,
        similarity: float,
        gene: str,
        disease: Optional[str],
        inheritance: Optional[str],
        patient_terms: Set[str],
        gene_terms: Set[str]
    ) -> tuple:
        """
        Determine evidence code based on similarity score.
        
        Returns tuple of (evidence_code, strength, explanation).
        
        Safety guards:
        - No evidence if total terms < min_terms (sparse data protection)
        - Explanation includes threshold info for transparency
        """
        # Calculate overlap and union for explanation
        overlap = patient_terms & gene_terms
        overlap_count = len(overlap)
        union = patient_terms | gene_terms
        total_terms = len(union)
        
        # Build disease/inheritance context string
        disease_str = disease if disease else f"{gene}-associated disease"
        inheritance_str = f", inheritance: {inheritance}" if inheritance else ""
        
        # Safety guard: Check for sparse phenotype data
        # If total terms are too few, similarity calculations may be unreliable
        if total_terms < self.min_terms:
            return (
                None,
                None,
                f"Phenotype information too sparse for reliable similarity-based evidence "
                f"(only {total_terms} total terms, minimum {self.min_terms} required)."
            )
        
        # PP4: Strong phenotype match (>= 0.8)
        pp4_threshold = self.thresholds['PP4']
        if similarity >= pp4_threshold:
            return (
                'PP4',
                'Supporting',
                f"Phenotype similarity {similarity:.2f} (>={pp4_threshold} PP4 threshold) "
                f"with {gene}-associated phenotype set ({disease_str}{inheritance_str}). "
                f"{overlap_count} of {len(patient_terms)} patient phenotypes match. "
                f"Evidence: PP4."
            )
        
        # PP4_supporting: Moderate match (0.5 - 0.8)
        pp4_supp_threshold = self.thresholds['PP4_SUPPORTING']
        if similarity >= pp4_supp_threshold:
            return (
                'PP4_supporting',
                'Supporting',
                f"Phenotype similarity {similarity:.2f} (>={pp4_supp_threshold} PP4_supporting threshold) "
                f"with {gene}-associated phenotype set ({disease_str}{inheritance_str}). "
                f"{overlap_count} of {len(patient_terms)} patient phenotypes match. "
                f"Evidence: PP4_supporting."
            )
        
        # BP5: Very low match (<= 0.2)
        bp5_threshold = self.thresholds['BP5']
        if similarity <= bp5_threshold:
            return (
                'BP5',
                'Supporting',
                f"Phenotype similarity {similarity:.2f} (<={bp5_threshold} BP5 threshold); "
                f"phenotype is inconsistent with {gene}-associated phenotype set "
                f"({disease_str}{inheritance_str}). Evidence: BP5."
            )
        
        # Neutral zone: no evidence
        return (
            None,
            None,
            f"Phenotype similarity {similarity:.2f} with {gene}; no phenotype-based evidence. "
            f"Similarity is between BP5 (<={bp5_threshold}) and PP4_supporting "
            f"(>={pp4_supp_threshold}) thresholds."
        )
    
    def get_gene_disease_info(self, gene: str) -> Optional[Dict[str, Any]]:
        """
        Get disease information for a gene.
        
        Args:
            gene: Gene symbol
            
        Returns:
            Dictionary with disease name, inheritance pattern, and phenotypes,
            or None if gene is unknown.
        """
        return self.gene_phenotype_db.get_gene_info(gene)
    
    def normalize_phenotypes(self, phenotypes: List[str]) -> Set[str]:
        """
        Normalize a list of phenotype terms to HPO IDs.
        
        Args:
            phenotypes: List of phenotype terms (HPO IDs or text)
            
        Returns:
            Set of normalized HPO IDs/tokens
        """
        return self.hpo_client.get_phenotype_terms(phenotypes)

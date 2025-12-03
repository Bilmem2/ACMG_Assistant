"""
Variant Data Structure
=====================

This module defines the VariantData class that holds all information
about a genetic variant being classified.

CANONICAL IMPLEMENTATION: This is the primary VariantData implementation.

NOTE: A legacy VariantData class also exists in core/evidence_evaluator.py.
TODO: Migrate all usages to this dataclass implementation and remove the
      legacy class from evidence_evaluator.py in a future refactor.

Changes in this documentation pass (Dec 2024):
- Added module-level documentation clarifying canonical status
- Added type hints to all public methods
- Improved docstrings with Args/Returns sections

Changes (Jan 2025):
- Added typed predictor_scores and population_stats fields for multi-source
  API data flow (PredictorScore and PopulationStats dataclasses)
"""

from typing import Dict, List, Optional, Any, TYPE_CHECKING
from dataclasses import dataclass, field
import json

# Import typed data structures (avoiding circular imports)
if TYPE_CHECKING:
    from config.predictors import PredictorScore, PopulationStats


@dataclass
class VariantData:
    @property
    def gene(self) -> Optional[str]:
        """Get gene symbol from basic_info."""
        return self.basic_info.get('gene', None)

    @property
    def hgvs_c(self) -> Optional[str]:
        """Get HGVS cDNA notation from basic_info."""
        return self.basic_info.get('hgvs_c', None)
    """
    Comprehensive data structure for genetic variant information.
    
    This class holds all the information needed for ACMG classification
    including basic variant info, population data, in silico predictions,
    genetic context, and functional studies.
    """
    
    # Basic variant information
    basic_info: Dict[str, Any] = field(default_factory=dict)
    
    # Population frequency data
    population_data: Dict[str, Any] = field(default_factory=dict)
    
    # In silico prediction scores
    insilico_data: Dict[str, Any] = field(default_factory=dict)
    
    # Genetic and inheritance data
    genetic_data: Dict[str, Any] = field(default_factory=dict)
    
    # Functional studies and segregation data
    functional_data: Dict[str, Any] = field(default_factory=dict)
    
    # Patient phenotype information for PP4/BP5 evaluation
    # Can be: List of HPO terms (HP:0001250), free-text phenotypes, or None
    patient_phenotypes: Optional[List[str]] = None
    
    # External database information
    clinvar_data: Dict[str, Any] = field(default_factory=dict)
    
    # Additional annotations
    annotations: Dict[str, Any] = field(default_factory=dict)
    
    # Metadata
    metadata: Dict[str, Any] = field(default_factory=dict)
    
    # ==========================================================================
    # Typed Multi-Source API Data (New Jan 2025)
    # ==========================================================================
    # These fields hold pre-fetched data from external APIs using typed
    # dataclasses. They enable the "pure interpreter" pattern where evaluators
    # only interpret data without making API calls themselves.
    
    # Typed predictor scores from multi-source APIs (PredictorScore objects)
    # Keys are predictor names: 'revel', 'cadd_phred', 'alphamissense', etc.
    predictor_scores: Optional[Dict[str, Any]] = None
    
    # Typed population statistics from multi-source APIs (PopulationStats objects)
    # Keys are population names: 'gnomad_v4', 'exac', 'topmed', etc.
    population_stats: Optional[Dict[str, Any]] = None
    
    def __post_init__(self):
        """Initialize metadata after object creation."""
        if not self.metadata:
            from datetime import datetime
            self.metadata = {
                'created_at': datetime.now().isoformat(),
                'version': '1.0.0',
                'source': 'ACMG Assistant'
            }
    
    def validate(self) -> List[str]:
        """
        Validate the variant data and return list of validation errors.
        
        Returns:
            List[str]: List of validation error messages
        """
        errors = []
        
        # Validate basic info
        if not self.basic_info:
            errors.append("Basic variant information is required")
        else:
            required_fields = ['gene', 'chromosome', 'position', 'ref_allele', 'alt_allele']
            for field in required_fields:
                if field not in self.basic_info or not self.basic_info[field]:
                    errors.append(f"Required field '{field}' is missing from basic info")
        
        # Validate chromosome format
        if 'chromosome' in self.basic_info:
            chr_val = str(self.basic_info['chromosome']).upper()
            valid_chromosomes = [str(i) for i in range(1, 23)] + ['X', 'Y', 'MT']
            if chr_val not in valid_chromosomes:
                errors.append(f"Invalid chromosome: {chr_val}")
        
        # Validate position
        if 'position' in self.basic_info:
            try:
                pos = int(self.basic_info['position'])
                if pos < 1:
                    errors.append("Position must be a positive integer")
            except (ValueError, TypeError):
                errors.append("Position must be a valid integer")
        
        # Validate alleles
        if 'ref_allele' in self.basic_info and 'alt_allele' in self.basic_info:
            ref = self.basic_info['ref_allele']
            alt = self.basic_info['alt_allele']
            valid_bases = set('ATCG')
            if not all(base in valid_bases for base in ref.upper()):
                errors.append(f"Invalid reference allele: {ref}")
            if not all(base in valid_bases for base in alt.upper()):
                errors.append(f"Invalid alternate allele: {alt}")
        
        # Validate population frequencies
        if self.population_data:
            for key, value in self.population_data.items():
                if key.endswith('_af') and value is not None:
                    try:
                        freq = float(value)
                        if not (0 <= freq <= 1):
                            errors.append(f"Allele frequency {key} must be between 0 and 1")
                    except (ValueError, TypeError):
                        errors.append(f"Allele frequency {key} must be a number")
        
        # Validate in silico scores
        if self.insilico_data:
            for predictor, score in self.insilico_data.items():
                if score is not None:
                    try:
                        score_val = float(score)
                        # Check predictor-specific ranges
                        if predictor in ['revel', 'alphamissense', 'metarnn', 'clinpred']:
                            if not (0 <= score_val <= 1):
                                errors.append(f"{predictor} score must be between 0 and 1")
                        elif predictor == 'cadd':
                            if not (0 <= score_val <= 50):
                                errors.append(f"CADD score must be between 0 and 50")
                    except (ValueError, TypeError):
                        errors.append(f"{predictor} score must be a number")
        
        return errors
    
    def get_variant_id(self) -> str:
        """
        Generate a unique identifier for this variant.
        
        Returns:
            str: Variant identifier in format chr:pos:ref:alt
        """
        if not all(key in self.basic_info for key in ['chromosome', 'position', 'ref_allele', 'alt_allele']):
            return "unknown_variant"
        
        return f"{self.basic_info['chromosome']}:{self.basic_info['position']}:{self.basic_info['ref_allele']}:{self.basic_info['alt_allele']}"
    
    def get_hgvs_nomenclature(self) -> Dict[str, str]:
        """
        Get HGVS nomenclature for this variant.
        
        Returns:
            Dict[str, str]: Dictionary with genomic, cDNA, and protein HGVS
        """
        hgvs = {}
        
        # Genomic HGVS
        if all(key in self.basic_info for key in ['chromosome', 'position', 'ref_allele', 'alt_allele']):
            chr_val = self.basic_info['chromosome']
            pos = self.basic_info['position']
            ref = self.basic_info['ref_allele']
            alt = self.basic_info['alt_allele']
            hgvs['genomic'] = f"chr{chr_val}:g.{pos}{ref}>{alt}"
        
        # cDNA HGVS
        if 'cdna_change' in self.basic_info:
            hgvs['cdna'] = self.basic_info['cdna_change']
        
        # Protein HGVS
        if 'protein_change' in self.basic_info:
            hgvs['protein'] = self.basic_info['protein_change']
        
        return hgvs
    
    def get_summary_info(self) -> Dict[str, Any]:
        """
        Get summary information about this variant.
        
        Returns:
            Dict[str, Any]: Summary information
        """
        summary = {
            'variant_id': self.get_variant_id(),
            'gene': self.basic_info.get('gene', 'Unknown'),
            'variant_type': self.basic_info.get('variant_type', 'Unknown'),
            'consequence': self.basic_info.get('consequence', 'Unknown'),
            'hgvs': self.get_hgvs_nomenclature(),
            'has_population_data': bool(self.population_data),
            'has_insilico_data': bool(self.insilico_data),
            'has_functional_data': bool(self.functional_data),
            'has_clinvar_data': bool(self.clinvar_data)
        }
        
        return summary
    
    def to_dict(self) -> Dict[str, Any]:
        """
        Convert variant data to dictionary format.
        
        Returns:
            Dict[str, Any]: Dictionary representation
        """
        result = {
            'basic_info': self.basic_info,
            'population_data': self.population_data,
            'insilico_data': self.insilico_data,
            'genetic_data': self.genetic_data,
            'functional_data': self.functional_data,
            'clinvar_data': self.clinvar_data,
            'annotations': self.annotations,
            'metadata': self.metadata
        }
        
        # Include typed API data if present (serialize to dicts)
        if self.predictor_scores:
            result['predictor_scores'] = {
                name: (
                    {'predictor': score.predictor, 'value': score.value, 
                     'source': score.source, 'version': score.version,
                     'is_inverted': score.is_inverted}
                    if hasattr(score, 'predictor') else score
                )
                for name, score in self.predictor_scores.items()
            }
        
        if self.population_stats:
            result['population_stats'] = {
                name: (
                    {'population': stats.population, 'af': stats.af,
                     'an': stats.an, 'ac': stats.ac, 'source': stats.source,
                     'version': stats.version}
                    if hasattr(stats, 'population') else stats
                )
                for name, stats in self.population_stats.items()
            }
        
        return result
    
    def to_json(self, indent: int = 2) -> str:
        """
        Convert variant data to JSON string.
        
        Args:
            indent (int): JSON indentation level
            
        Returns:
            str: JSON representation
        """
        return json.dumps(self.to_dict(), indent=indent, default=str)
    
    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'VariantData':
        """
        Create VariantData instance from dictionary.
        
        Args:
            data (Dict[str, Any]): Dictionary containing variant data
            
        Returns:
            VariantData: New instance
        """
        instance = cls(
            basic_info=data.get('basic_info', {}),
            population_data=data.get('population_data', {}),
            insilico_data=data.get('insilico_data', {}),
            genetic_data=data.get('genetic_data', {}),
            functional_data=data.get('functional_data', {}),
            clinvar_data=data.get('clinvar_data', {}),
            annotations=data.get('annotations', {}),
            metadata=data.get('metadata', {})
        )
        
        # Restore typed API data if present
        if 'predictor_scores' in data and data['predictor_scores']:
            instance.predictor_scores = data['predictor_scores']
        if 'population_stats' in data and data['population_stats']:
            instance.population_stats = data['population_stats']
        
        return instance
    
    @classmethod
    def from_json(cls, json_str: str) -> 'VariantData':
        """
        Create VariantData instance from JSON string.
        
        Args:
            json_str (str): JSON string containing variant data
            
        Returns:
            VariantData: New instance
        """
        data = json.loads(json_str)
        return cls.from_dict(data)
    
    def update_field(self, category: str, field: str, value: Any) -> None:
        """
        Update a specific field in the variant data.
        
        Args:
            category (str): Data category (e.g., 'basic_info', 'population_data')
            field (str): Field name
            value (Any): New value
        """
        if category == 'basic_info':
            self.basic_info[field] = value
        elif category == 'population_data':
            self.population_data[field] = value
        elif category == 'insilico_data':
            self.insilico_data[field] = value
        elif category == 'genetic_data':
            self.genetic_data[field] = value
        elif category == 'functional_data':
            self.functional_data[field] = value
        elif category == 'clinvar_data':
            self.clinvar_data[field] = value
        elif category == 'annotations':
            self.annotations[field] = value
        elif category == 'metadata':
            self.metadata[field] = value
        else:
            raise ValueError(f"Unknown category: {category}")
    
    def get_field(self, category: str, field: str, default: Any = None) -> Any:
        """
        Get a specific field from the variant data.
        
        Args:
            category (str): Data category
            field (str): Field name
            default (Any): Default value if field not found
            
        Returns:
            Any: Field value or default
        """
        if category == 'basic_info':
            return self.basic_info.get(field, default)
        elif category == 'population_data':
            return self.population_data.get(field, default)
        elif category == 'insilico_data':
            return self.insilico_data.get(field, default)
        elif category == 'genetic_data':
            return self.genetic_data.get(field, default)
        elif category == 'functional_data':
            return self.functional_data.get(field, default)
        elif category == 'clinvar_data':
            return self.clinvar_data.get(field, default)
        elif category == 'annotations':
            return self.annotations.get(field, default)
        elif category == 'metadata':
            return self.metadata.get(field, default)
        else:
            return default
    
    # ==========================================================================
    # Multi-Source API Data Helpers (New Jan 2025)
    # ==========================================================================
    
    def has_predictor_scores(self) -> bool:
        """Check if typed predictor scores are available."""
        return bool(self.predictor_scores)
    
    def get_predictor_score(self, predictor: str) -> Optional[Any]:
        """
        Get a typed predictor score by name.
        
        Args:
            predictor: Predictor name (e.g., 'revel', 'cadd_phred')
            
        Returns:
            PredictorScore object or None if not available
        """
        if not self.predictor_scores:
            return None
        return self.predictor_scores.get(predictor)
    
    def get_available_predictors(self) -> List[str]:
        """Get list of predictors with valid scores."""
        if not self.predictor_scores:
            return []
        return [
            name for name, score in self.predictor_scores.items()
            if score is not None and (
                hasattr(score, 'value') and score.value is not None
                or isinstance(score, (int, float))
            )
        ]
    
    def has_population_stats(self) -> bool:
        """Check if typed population statistics are available."""
        return bool(self.population_stats)
    
    def get_population_stat(self, population: str) -> Optional[Any]:
        """
        Get population statistics by name.
        
        Args:
            population: Population name (e.g., 'gnomad_v4', 'exac')
            
        Returns:
            PopulationStats object or None if not available
        """
        if not self.population_stats:
            return None
        return self.population_stats.get(population)
    
    def get_max_population_af(self) -> Optional[float]:
        """
        Get maximum allele frequency across all population sources.
        
        Returns:
            Maximum AF or None if no population data
        """
        if not self.population_stats:
            return None
        
        afs = []
        for stats in self.population_stats.values():
            if hasattr(stats, 'af') and stats.af is not None:
                afs.append(stats.af)
            elif isinstance(stats, dict) and stats.get('af') is not None:
                afs.append(stats['af'])
        
        return max(afs) if afs else None
    
    def is_absent_from_populations(self) -> bool:
        """
        Check if variant is absent from all population databases.
        
        Returns:
            True if absent from all populations or no data available
        """
        if not self.population_stats:
            return True
        
        for stats in self.population_stats.values():
            if hasattr(stats, 'is_absent'):
                if not stats.is_absent():
                    return False
            elif isinstance(stats, dict):
                af = stats.get('af')
                if af is not None and af > 0:
                    return False
        
        return True

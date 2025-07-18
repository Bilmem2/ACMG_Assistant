
"""
Variant Data Structure
=====================

This module defines the VariantData class that holds all information
about a genetic variant being classified.
"""

from typing import Dict, List, Optional, Any
from dataclasses import dataclass, field
import json


@dataclass
class VariantData:
    @property
    def gene(self):
        return self.basic_info.get('gene', None)

    @property
    def hgvs_c(self):
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
    
    # External database information
    clinvar_data: Dict[str, Any] = field(default_factory=dict)
    
    # Additional annotations
    annotations: Dict[str, Any] = field(default_factory=dict)
    
    # Metadata
    metadata: Dict[str, Any] = field(default_factory=dict)
    
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
        return {
            'basic_info': self.basic_info,
            'population_data': self.population_data,
            'insilico_data': self.insilico_data,
            'genetic_data': self.genetic_data,
            'functional_data': self.functional_data,
            'clinvar_data': self.clinvar_data,
            'annotations': self.annotations,
            'metadata': self.metadata
        }
    
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
        return cls(
            basic_info=data.get('basic_info', {}),
            population_data=data.get('population_data', {}),
            insilico_data=data.get('insilico_data', {}),
            genetic_data=data.get('genetic_data', {}),
            functional_data=data.get('functional_data', {}),
            clinvar_data=data.get('clinvar_data', {}),
            annotations=data.get('annotations', {}),
            metadata=data.get('metadata', {})
        )
    
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

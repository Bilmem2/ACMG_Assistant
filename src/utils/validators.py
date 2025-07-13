"""
Validators Module
================

This module contains validation functions for various data types
used in the ACMG classification process.
"""

import re
from typing import Any, Dict, List, Optional, Tuple
from config.constants import VALIDATION_PATTERNS, VARIANT_CONSEQUENCES


class Validators:
    """
    Collection of validation methods for variant data.
    
    Provides comprehensive validation for all types of input data
    used in ACMG variant classification.
    """
    
    def __init__(self):
        """Initialize validators."""
        pass
    
    def validate_basic_info(self, basic_info: Dict[str, Any]) -> List[str]:
        """
        Validate basic variant information.
        
        Args:
            basic_info (Dict[str, Any]): Basic variant information
            
        Returns:
            List[str]: List of validation errors
        """
        errors = []
        
        # Required fields
        required_fields = ['gene', 'chromosome', 'position', 'ref_allele', 'alt_allele']
        for field in required_fields:
            if field not in basic_info or not basic_info[field]:
                errors.append(f"Required field '{field}' is missing")
        
        # Validate individual fields
        if 'gene' in basic_info:
            gene_errors = self.validate_gene_symbol(basic_info['gene'])
            errors.extend(gene_errors)
        
        if 'chromosome' in basic_info:
            chr_errors = self.validate_chromosome(basic_info['chromosome'])
            errors.extend(chr_errors)
        
        if 'position' in basic_info:
            pos_errors = self.validate_position(basic_info['position'])
            errors.extend(pos_errors)
        
        if 'ref_allele' in basic_info:
            ref_errors = self.validate_allele(basic_info['ref_allele'], 'reference')
            errors.extend(ref_errors)
        
        if 'alt_allele' in basic_info:
            alt_errors = self.validate_allele(basic_info['alt_allele'], 'alternate')
            errors.extend(alt_errors)
        
        # Validate HGVS nomenclature if provided
        if 'cdna_change' in basic_info and basic_info['cdna_change']:
            hgvs_errors = self.validate_hgvs_cdna(basic_info['cdna_change'])
            errors.extend(hgvs_errors)
        
        if 'protein_change' in basic_info and basic_info['protein_change']:
            hgvs_errors = self.validate_hgvs_protein(basic_info['protein_change'])
            errors.extend(hgvs_errors)
        
        return errors
    
    def validate_gene_symbol(self, gene_symbol: str) -> List[str]:
        """
        Validate gene symbol format.
        
        Args:
            gene_symbol (str): Gene symbol
            
        Returns:
            List[str]: List of validation errors
        """
        errors = []
        
        if not gene_symbol:
            errors.append("Gene symbol is required")
            return errors
        
        # Check format
        if not re.match(VALIDATION_PATTERNS['gene_symbol'], gene_symbol.upper()):
            errors.append(f"Invalid gene symbol format: {gene_symbol}")
        
        # Check length
        if len(gene_symbol) > 20:
            errors.append(f"Gene symbol too long: {gene_symbol}")
        
        return errors
    
    def validate_chromosome(self, chromosome: str) -> List[str]:
        """
        Validate chromosome format.
        
        Args:
            chromosome (str): Chromosome
            
        Returns:
            List[str]: List of validation errors
        """
        errors = []
        
        if not chromosome:
            errors.append("Chromosome is required")
            return errors
        
        # Normalize chromosome format
        chr_val = str(chromosome).upper().replace('CHR', '')
        # Convert M to MT for mitochondrial chromosome
        if chr_val == 'M':
            chr_val = 'MT'
        
        # Check if valid chromosome
        valid_chromosomes = [str(i) for i in range(1, 23)] + ['X', 'Y', 'MT']
        if chr_val not in valid_chromosomes:
            errors.append(f"Invalid chromosome: {chromosome}")
        
        return errors
    
    def validate_position(self, position: Any) -> List[str]:
        """
        Validate genomic position.
        
        Args:
            position (Any): Genomic position
            
        Returns:
            List[str]: List of validation errors
        """
        errors = []
        
        if position is None:
            errors.append("Position is required")
            return errors
        
        try:
            pos_int = int(position)
            if pos_int < 1:
                errors.append("Position must be a positive integer")
            elif pos_int > 250000000:  # Approximate max chromosome length
                errors.append("Position is too large")
        except (ValueError, TypeError):
            errors.append(f"Position must be an integer: {position}")
        
        return errors
    
    def validate_allele(self, allele: str, allele_type: str) -> List[str]:
        """
        Validate allele format.
        
        Args:
            allele (str): Allele sequence
            allele_type (str): Type of allele (reference or alternate)
            
        Returns:
            List[str]: List of validation errors
        """
        errors = []
        
        if not allele:
            errors.append(f"{allele_type.capitalize()} allele is required")
            return errors
        
        # Check for valid DNA bases
        valid_bases = set('ATCG')
        if not set(allele.upper()).issubset(valid_bases):
            errors.append(f"Invalid {allele_type} allele - contains non-DNA bases: {allele}")
        
        # Check length
        if len(allele) > 1000:
            errors.append(f"{allele_type.capitalize()} allele too long: {len(allele)} bases")
        
        return errors
    
    def validate_hgvs_cdna(self, hgvs_cdna: str) -> List[str]:
        """
        Validate HGVS cDNA nomenclature.
        
        Args:
            hgvs_cdna (str): HGVS cDNA change
            
        Returns:
            List[str]: List of validation errors
        """
        errors = []
        
        if not hgvs_cdna:
            return errors  # Optional field
        
        # Basic format check
        if not hgvs_cdna.startswith('c.'):
            errors.append("HGVS cDNA must start with 'c.'")
        
        # More detailed validation could be added here
        # For now, just check basic pattern
        if not re.match(VALIDATION_PATTERNS['hgvs_cdna'], hgvs_cdna):
            errors.append(f"Invalid HGVS cDNA format: {hgvs_cdna}")
        
        return errors
    
    def validate_hgvs_protein(self, hgvs_protein: str) -> List[str]:
        """
        Validate HGVS protein nomenclature.
        
        Args:
            hgvs_protein (str): HGVS protein change
            
        Returns:
            List[str]: List of validation errors
        """
        errors = []
        
        if not hgvs_protein:
            return errors  # Optional field
        
        # Basic format check
        if not hgvs_protein.startswith('p.'):
            errors.append("HGVS protein must start with 'p.'")
        
        # More detailed validation could be added here
        if not re.match(VALIDATION_PATTERNS['hgvs_protein'], hgvs_protein):
            errors.append(f"Invalid HGVS protein format: {hgvs_protein}")
        
        return errors
    
    def validate_population_data(self, population_data: Dict[str, Any]) -> List[str]:
        """
        Validate population frequency data.
        
        Args:
            population_data (Dict[str, Any]): Population frequency data
            
        Returns:
            List[str]: List of validation errors
        """
        errors = []
        
        # Validate allele frequencies
        frequency_fields = ['gnomad_af', 'gnomad_af_popmax', 'exac_af']
        for field in frequency_fields:
            if field in population_data and population_data[field] is not None:
                freq_errors = self.validate_allele_frequency(population_data[field], field)
                errors.extend(freq_errors)
        
        # Validate counts
        count_fields = ['gnomad_hom_count', 'gnomad_het_count']
        for field in count_fields:
            if field in population_data and population_data[field] is not None:
                count_errors = self.validate_count(population_data[field], field)
                errors.extend(count_errors)
        
        return errors
    
    def validate_allele_frequency(self, frequency: Any, field_name: str) -> List[str]:
        """
        Validate allele frequency value.
        
        Args:
            frequency (Any): Allele frequency
            field_name (str): Field name for error messages
            
        Returns:
            List[str]: List of validation errors
        """
        errors = []
        
        try:
            freq_val = float(frequency)
            if not (0 <= freq_val <= 1):
                errors.append(f"{field_name} must be between 0 and 1")
        except (ValueError, TypeError):
            errors.append(f"{field_name} must be a number")
        
        return errors
    
    def validate_count(self, count: Any, field_name: str) -> List[str]:
        """
        Validate count value.
        
        Args:
            count (Any): Count value
            field_name (str): Field name for error messages
            
        Returns:
            List[str]: List of validation errors
        """
        errors = []
        
        try:
            count_val = int(count)
            if count_val < 0:
                errors.append(f"{field_name} must be non-negative")
        except (ValueError, TypeError):
            errors.append(f"{field_name} must be an integer")
        
        return errors
    
    def validate_insilico_data(self, insilico_data: Dict[str, Any]) -> List[str]:
        """
        Validate in silico prediction data.
        
        Args:
            insilico_data (Dict[str, Any]): In silico prediction data
            
        Returns:
            List[str]: List of validation errors
        """
        errors = []
        
        # Validate individual scores
        score_ranges = {
            'revel': (0, 1),
            'alphamissense': (0, 1),
            'metarnn': (0, 1),
            'clinpred': (0, 1),
            'sift': (0, 1),
            'polyphen2': (0, 1),
            'mutationtaster': (0, 1),
            'cadd': (0, 50),
            'bayesdel': (-1, 1)
        }
        
        for predictor, score in insilico_data.items():
            if score is not None and predictor in score_ranges:
                min_val, max_val = score_ranges[predictor]
                score_errors = self.validate_score_range(score, predictor, min_val, max_val)
                errors.extend(score_errors)
        
        return errors
    
    def validate_score_range(self, score: Any, predictor: str, 
                           min_val: float, max_val: float) -> List[str]:
        """
        Validate score within specified range.
        
        Args:
            score (Any): Score value
            predictor (str): Predictor name
            min_val (float): Minimum valid value
            max_val (float): Maximum valid value
            
        Returns:
            List[str]: List of validation errors
        """
        errors = []
        
        try:
            score_val = float(score)
            if not (min_val <= score_val <= max_val):
                errors.append(f"{predictor} score must be between {min_val} and {max_val}")
        except (ValueError, TypeError):
            errors.append(f"{predictor} score must be a number")
        
        return errors
    
    def validate_genetic_data(self, genetic_data: Dict[str, Any]) -> List[str]:
        """
        Validate genetic and inheritance data.
        
        Args:
            genetic_data (Dict[str, Any]): Genetic data
            
        Returns:
            List[str]: List of validation errors
        """
        errors = []
        
        # Validate inheritance pattern
        if 'inheritance' in genetic_data:
            valid_patterns = ['AD', 'AR', 'XLD', 'XLR', 'MT', 'unknown']
            if genetic_data['inheritance'] not in valid_patterns:
                errors.append(f"Invalid inheritance pattern: {genetic_data['inheritance']}")
        
        # Validate zygosity
        if 'zygosity' in genetic_data:
            valid_zygosity = ['homozygous', 'heterozygous', 'hemizygous', 'unknown']
            if genetic_data['zygosity'] not in valid_zygosity:
                errors.append(f"Invalid zygosity: {genetic_data['zygosity']}")
        
        # Check logical consistency
        if (genetic_data.get('inheritance') == 'AR' and 
            genetic_data.get('zygosity') == 'heterozygous'):
            if 'compound_het' not in genetic_data:
                errors.append("Compound heterozygous status required for AR inheritance with heterozygous variant")
        
        return errors
    
    def validate_functional_data(self, functional_data: Dict[str, Any]) -> List[str]:
        """
        Validate functional studies data.
        
        Args:
            functional_data (Dict[str, Any]): Functional data
            
        Returns:
            List[str]: List of validation errors
        """
        errors = []
        
        # Validate segregation data
        if 'segregation' in functional_data:
            valid_segregation = ['cosegregates', 'does_not_segregate', 'insufficient_data', 'not_performed']
            if functional_data['segregation'] not in valid_segregation:
                errors.append(f"Invalid segregation value: {functional_data['segregation']}")
        
        # Validate de novo status
        if 'denovo' in functional_data:
            valid_denovo = ['confirmed', 'assumed', 'not_denovo', 'unknown']
            if functional_data['denovo'] not in valid_denovo:
                errors.append(f"Invalid de novo status: {functional_data['denovo']}")
        
        # Validate functional studies
        if 'functional_studies' in functional_data:
            valid_functional = ['damaging', 'benign', 'inconclusive', 'not_performed']
            if functional_data['functional_studies'] not in valid_functional:
                errors.append(f"Invalid functional studies result: {functional_data['functional_studies']}")
        
        # Validate case-control data
        if functional_data.get('case_control') == 'yes':
            required_fields = ['cases_with_variant', 'total_cases', 'controls_with_variant', 'total_controls']
            for field in required_fields:
                if field not in functional_data or functional_data[field] is None:
                    errors.append(f"Case-control field '{field}' is required")
                else:
                    count_errors = self.validate_count(functional_data[field], field)
                    errors.extend(count_errors)
        
        return errors
    
    def validate_variant_consistency(self, variant_data: Dict[str, Any]) -> List[str]:
        """
        Validate consistency between different data fields.
        
        Args:
            variant_data (Dict[str, Any]): Complete variant data
            
        Returns:
            List[str]: List of validation errors
        """
        errors = []
        
        # Check if reference and alternate alleles are different
        basic_info = variant_data.get('basic_info', {})
        if (basic_info.get('ref_allele') and basic_info.get('alt_allele') and
            basic_info['ref_allele'].upper() == basic_info['alt_allele'].upper()):
            errors.append("Reference and alternate alleles cannot be identical")
        
        # Check variant type consistency
        variant_type = basic_info.get('variant_type')
        ref_allele = basic_info.get('ref_allele', '')
        alt_allele = basic_info.get('alt_allele', '')
        
        if variant_type and ref_allele and alt_allele:
            if variant_type == 'missense' and len(ref_allele) != len(alt_allele):
                errors.append("Missense variants should have equal-length alleles")
            elif variant_type == 'synonymous' and len(ref_allele) != len(alt_allele):
                errors.append("Synonymous variants should have equal-length alleles")
        
        # Check population frequency consistency
        population_data = variant_data.get('population_data', {})
        gnomad_af = population_data.get('gnomad_af')
        gnomad_af_popmax = population_data.get('gnomad_af_popmax')
        
        if gnomad_af is not None and gnomad_af_popmax is not None:
            if gnomad_af > gnomad_af_popmax:
                errors.append("Overall gnomAD frequency cannot be higher than popmax frequency")
        
        return errors
    
    def validate_all_data(self, variant_data: Dict[str, Any]) -> Dict[str, List[str]]:
        """
        Validate all variant data comprehensively.
        
        Args:
            variant_data (Dict[str, Any]): Complete variant data
            
        Returns:
            Dict[str, List[str]]: Validation results by category
        """
        results = {
            'basic_info': [],
            'population_data': [],
            'insilico_data': [],
            'genetic_data': [],
            'functional_data': [],
            'consistency': []
        }
        
        # Validate each category
        if 'basic_info' in variant_data:
            results['basic_info'] = self.validate_basic_info(variant_data['basic_info'])
        
        if 'population_data' in variant_data:
            results['population_data'] = self.validate_population_data(variant_data['population_data'])
        
        if 'insilico_data' in variant_data:
            results['insilico_data'] = self.validate_insilico_data(variant_data['insilico_data'])
        
        if 'genetic_data' in variant_data:
            results['genetic_data'] = self.validate_genetic_data(variant_data['genetic_data'])
        
        if 'functional_data' in variant_data:
            results['functional_data'] = self.validate_functional_data(variant_data['functional_data'])
        
        # Validate consistency
        results['consistency'] = self.validate_variant_consistency(variant_data)
        
        return results
    
    def has_errors(self, validation_results: Dict[str, List[str]]) -> bool:
        """
        Check if validation results contain any errors.
        
        Args:
            validation_results (Dict[str, List[str]]): Validation results
            
        Returns:
            bool: True if errors found
        """
        return any(errors for errors in validation_results.values())
    
    def get_all_errors(self, validation_results: Dict[str, List[str]]) -> List[str]:
        """
        Get all validation errors as a flat list.
        
        Args:
            validation_results (Dict[str, List[str]]): Validation results
            
        Returns:
            List[str]: All errors
        """
        all_errors = []
        for category, errors in validation_results.items():
            for error in errors:
                all_errors.append(f"{category}: {error}")
        return all_errors

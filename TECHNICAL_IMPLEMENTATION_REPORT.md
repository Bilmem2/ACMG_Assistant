# ACMG Variant Classification Assistant v3.0.0 - Technical Implementation Report

**Prepared for**: Academic and Clinical Genetics Community  
**Date**: July 9, 2025  
**Version**: 3.0.0  
**Author**: Can Sevilmiş  
**Repository**: https://github.com/Bilmem2/acmg-assessor

---

## EXECUTIVE SUMMARY

This technical report provides a comprehensive overview of the ACMG Variant Classification Assistant v3.0.0, detailing the complete implementation of all 28 ACMG/AMP evidence criteria, computational methodologies, statistical frameworks, and algorithmic approaches used in automated variant pathogenicity assessment. This document serves as both a technical reference for researchers and clinicians seeking to understand the underlying methodologies and a validation resource for academic review and clinical implementation.

**Key Implementation Highlights:**
- Complete 28/28 ACMG/AMP criteria implementation (100% coverage)
- VAMPP-score inspired computational metascore framework
- Interactive literature-based evidence evaluation system
- Population genetics integration with gene-specific thresholds
- Comprehensive quality control and validation mechanisms

## 1. OVERVIEW AND METHODOLOGY

### 1.1 Algorithm Purpose
The ACMG Variant Classification Assistant implements the complete 28 ACMG/AMP evidence criteria framework as defined by Richards et al. (2015) and updated by ClinGen SVI working groups (2023). The algorithm provides systematic, evidence-based variant pathogenicity assessment using computational predictions, population frequency data, and interactive clinical evaluation.

### 1.2 Core Classification Framework
The algorithm operates on a multi-tiered evidence evaluation system:

**Evidence Strength Hierarchy:**
- **Very Strong (PVS)**: 8 points
- **Strong (PS)**: 4 points  
- **Moderate (PM)**: 2 points
- **Supporting (PP)**: 1 point
- **Stand Alone Benign (BA)**: -8 points
- **Strong Benign (BS)**: -4 points
- **Supporting Benign (BP)**: -1 point

**Final Classification Logic:**
- **Pathogenic**: ≥8 points (e.g., PVS1 + PS1 OR PS1 + PS2 + PM1 + PM2)
- **Likely Pathogenic**: 4-7 points
- **Uncertain Significance**: -3 to 3 points
- **Likely Benign**: -7 to -4 points
- **Benign**: ≤-8 points

---

## 2. DETAILED CRITERIA IMPLEMENTATION

### 2.1 PATHOGENIC EVIDENCE CRITERIA

#### PVS1 - Null Variants in LOF Intolerant Genes
**Implementation:**
```python
def _evaluate_pvs1(variant_data):
    variant_type = variant_data.variant_type.lower()
    consequence = variant_data.consequence.lower()
    
    # LOF variant types
    lof_types = ['nonsense', 'frameshift', 'splice_donor', 'splice_acceptor', 
                 'start_lost', 'stop_gained']
    
    # SpliceAI integration for splice variants
    if 'splice' in consequence:
        spliceai_scores = [ag_score, al_score, dg_score, dl_score]
        max_splice_score = max(spliceai_scores) if spliceai_scores else 0
        if max_splice_score > 0.5:  # High confidence splice disruption
            return {'applies': True, 'strength': 'Very Strong'}
```

**Statistical Thresholds:**
- SpliceAI threshold: >0.5 for high confidence splice disruption
- Consequence annotation: Direct mapping from VEP/SnpEff annotations

#### PS2/PM6 - De Novo Variants
**Implementation:**
```python
def _evaluate_ps2_pm6(variant_data):
    de_novo_status = variant_data.de_novo_status
    parental_confirmation = variant_data.parental_confirmation
    
    if de_novo_status == 'confirmed' and parental_confirmation == 'confirmed':
        return {'applies': True, 'strength': 'Strong', 'criterion': 'PS2'}
    elif de_novo_status == 'assumed' and parental_confirmation == 'not_confirmed':
        return {'applies': True, 'strength': 'Moderate', 'criterion': 'PM6'}
```

**Logic:**
- **PS2**: Confirmed de novo with verified paternity/maternity
- **PM6**: Assumed de novo without parental confirmation

#### PM2 - Population Frequency Analysis
**Implementation:**
```python
def _evaluate_pm2(variant_data):
    gnomad_af = variant_data.gnomad_af
    gene = variant_data.gene
    
    # Gene-specific frequency thresholds
    gene_thresholds = {
        'BRCA1': 0.00001,   # Very rare disease genes
        'BRCA2': 0.00001,
        'TP53': 0.000001,   # Ultra-rare for dominant conditions
        'default': 0.0001   # General threshold
    }
    
    threshold = gene_thresholds.get(gene, gene_thresholds['default'])
    
    if gnomad_af is None or gnomad_af < threshold:
        return {'applies': True, 'strength': 'Moderate'}
```

**Statistical Framework:**
- Default threshold: <0.01% (1 in 10,000)
- Gene-specific thresholds based on disease prevalence
- Uses gnomAD v3.1.2 population frequencies

#### PP3/BP4 - Computational Evidence (Enhanced Metascore)
**Implementation:**
```python
def _calculate_metascore(variant_data):
    predictors = variant_data.insilico_data
    
    # Weighted scoring system inspired by VAMPP-score
    weights = {
        'revel_score': 0.25,      # Primary ensemble predictor
        'cadd_phred': 0.20,       # Conservation + functional impact
        'alphamissense': 0.15,    # Protein structure-based
        'sift_score': 0.10,       # Evolutionary conservation (inverted)
        'polyphen2': 0.10,        # Structural impact
        'metasvm': 0.10,          # SVM ensemble
        'vest4': 0.05,            # Variant effect
        'fathmm': 0.05            # Hidden Markov Model
    }
    
    weighted_sum = 0
    total_weight = 0
    
    for predictor, weight in weights.items():
        if predictor in predictors and predictors[predictor] is not None:
            score = predictors[predictor]
            
            # Normalize scores to 0-1 range
            if predictor == 'sift_score':
                score = 1 - score  # SIFT is inverted (lower = more damaging)
            elif predictor == 'cadd_phred':
                score = min(score / 40.0, 1.0)  # CADD normalization
            
            weighted_sum += score * weight
            total_weight += weight
    
    if total_weight > 0:
        metascore = weighted_sum / total_weight
        
        # Apply evidence criteria
        if metascore >= 0.7:
            return {'criterion': 'PP3', 'applies': True, 'score': metascore}
        elif metascore <= 0.3:
            return {'criterion': 'BP4', 'applies': True, 'score': metascore}
    
    return {'applies': False, 'score': None}
```

**Statistical Methodology:**
- VAMPP-score inspired weighted ensemble approach
- Dynamic predictor weighting based on availability
- Normalization to 0-1 scale for all predictors
- Evidence thresholds: PP3 ≥0.7, BP4 ≤0.3

### 2.2 BENIGN EVIDENCE CRITERIA

#### BA1/BS1 - Population Frequency Analysis
**Implementation:**
```python
def _evaluate_ba1_bs1(variant_data):
    gnomad_af = variant_data.gnomad_af
    disease_prevalence = variant_data.disease_prevalence
    
    # BA1: Allele frequency >5% in population databases
    if gnomad_af and gnomad_af > 0.05:
        return {'applies': True, 'strength': 'Stand Alone', 'criterion': 'BA1'}
    
    # BS1: Frequency greater than expected for disorder
    if disease_prevalence:
        expected_max_af = disease_prevalence * 2  # Assuming full penetrance
        if gnomad_af and gnomad_af > expected_max_af:
            return {'applies': True, 'strength': 'Strong', 'criterion': 'BS1'}
```

**Statistical Framework:**
- **BA1**: Fixed 5% threshold for common variants
- **BS1**: Disease prevalence-based calculation with penetrance adjustment

#### BP7 - Synonymous Variants with Splice Analysis
**Implementation:**
```python
def _evaluate_bp7(variant_data):
    if variant_data.variant_type == 'synonymous':
        spliceai_scores = [
            variant_data.spliceai_ag_score,
            variant_data.spliceai_al_score, 
            variant_data.spliceai_dg_score,
            variant_data.spliceai_dl_score
        ]
        
        max_splice_score = max([s for s in spliceai_scores if s is not None])
        
        if max_splice_score < 0.1:  # Low splice impact
            return {'applies': True, 'strength': 'Supporting'}
```

**Threshold:**
- SpliceAI maximum score <0.1 for no predicted splice impact

---

## 3. INTERACTIVE CRITERIA EVALUATION

### 3.1 Literature-Based Criteria (PS1, PM5, PP4, etc.)
These criteria require manual literature review and are implemented with guided user interaction:

**PS1 Implementation:**
```python
def _evaluate_ps1_interactive(variant_data):
    gene = variant_data.gene
    aa_change = variant_data.amino_acid_change
    
    # Provide search recommendations
    search_terms = [
        f"PubMed: '{gene} {aa_change} pathogenic'",
        f"ClinVar: '{gene} {aa_change}'",
        f"HGMD: Same amino acid change"
    ]
    
    # User evaluation with guided questions
    user_response = prompt_user_with_guidance(
        question="Have you found the exact same amino acid change reported as pathogenic?",
        search_recommendations=search_terms,
        guidance="PS1 applies if the EXACT same amino acid change has been reported as pathogenic with sufficient evidence."
    )
    
    return process_user_response(user_response, 'PS1', 'Strong')
```

### 3.2 Guidance Framework
Each interactive criterion provides:
- Specific search recommendations
- Literature database suggestions
- Evidence evaluation guidance
- Conservative application principles

---

## 4. STATISTICAL METHODS AND THRESHOLDS

### 4.1 Population Frequency Statistics
**Data Sources:**
- gnomAD v3.1.2 (primary)
- 1000 Genomes Project
- ESP6500
- ExAC

**Frequency Calculations:**
```python
# Gene-specific frequency thresholds
frequency_thresholds = {
    'autosomal_dominant': {
        'high_penetrance': disease_prevalence / 2,
        'moderate_penetrance': disease_prevalence / 1.5,
        'low_penetrance': disease_prevalence
    },
    'autosomal_recessive': {
        'carrier_frequency': sqrt(disease_prevalence),
        'pathogenic_threshold': disease_prevalence / 100
    }
}
```

### 4.2 Conservation Analysis
**PhyloP Scoring:**
- 100-way vertebrate conservation: threshold >2.5
- 30-way mammalian conservation: threshold >1.5
- 17-way primate conservation: threshold >1.0

**GERP++ Scoring:**
- Rejected substitutions score: threshold >4.0
- Element score: threshold >2.0

### 4.3 Splice Prediction Integration
**SpliceAI v1.3:**
- Acceptor gain (AG): threshold >0.5 for high confidence
- Acceptor loss (AL): threshold >0.5 for high confidence
- Donor gain (DG): threshold >0.5 for high confidence
- Donor loss (DL): threshold >0.5 for high confidence

**MaxEntScan Integration:**
- Reference score vs. alternate score difference
- Threshold: >10% reduction for splice impact

---

## 5. COMPUTATIONAL PREDICTOR INTEGRATION

### 5.1 Primary Ensemble Predictors

#### REVEL (Rare Exome Variant Ensemble Learner)
**Implementation:**
- Weight: 25% of metascore
- Range: 0-1 (higher = more pathogenic)
- Threshold: >0.75 for pathogenic evidence

#### CADD (Combined Annotation Dependent Depletion)
**Implementation:**
- Uses CADD-Phred scores
- Normalization: score/40 (capped at 1.0)
- Weight: 20% of metascore
- Threshold: >20 for potential pathogenicity

#### AlphaMissense (DeepMind)
**Implementation:**
- Protein structure-based prediction
- Range: 0-1 (higher = more pathogenic)
- Weight: 15% of metascore
- Threshold: >0.564 for likely pathogenic

### 5.2 Traditional Predictors

#### SIFT (Sorting Intolerant From Tolerant)
**Implementation:**
- Score inversion: 1 - SIFT_score (for consistency)
- Weight: 10% of metascore
- Threshold: <0.05 original score (>0.95 inverted)

#### PolyPhen-2 (HumDiv)
**Implementation:**
- Uses HumDiv model for disease variants
- Range: 0-1 (higher = more damaging)
- Weight: 10% of metascore
- Categories: Benign (<0.15), Possibly Damaging (0.15-0.85), Probably Damaging (>0.85)

### 5.3 Metascore Calculation Algorithm
```python
def calculate_enhanced_metascore(predictors, variant_type='missense'):
    """
    Enhanced metascore calculation with dynamic weighting
    """
    
    # Variant-type specific weights
    if variant_type == 'missense':
        weights = MISSENSE_WEIGHTS
    elif variant_type == 'splice':
        weights = SPLICE_WEIGHTS
    else:
        weights = DEFAULT_WEIGHTS
    
    # Calculate weighted average with available predictors
    total_score = 0
    total_weight = 0
    
    for predictor, value in predictors.items():
        if value is not None and predictor in weights:
            normalized_score = normalize_predictor_score(predictor, value)
            weight = weights[predictor]
            
            total_score += normalized_score * weight
            total_weight += weight
    
    if total_weight == 0:
        return None
    
    metascore = total_score / total_weight
    
    # Apply confidence adjustment based on number of predictors
    confidence_factor = min(len([p for p in predictors.values() if p is not None]) / 5, 1.0)
    adjusted_score = metascore * confidence_factor
    
    return adjusted_score
```

---

## 6. EVIDENCE CONSOLIDATION AND CLASSIFICATION

### 6.1 Evidence Strength Mapping
```python
EVIDENCE_WEIGHTS = {
    'PVS1': 8,   # Very Strong Pathogenic
    'PS1': 4, 'PS2': 4, 'PS3': 4, 'PS4': 4,  # Strong Pathogenic
    'PM1': 2, 'PM2': 2, 'PM3': 2, 'PM4': 2, 'PM5': 2, 'PM6': 2,  # Moderate Pathogenic
    'PP1': 1, 'PP2': 1, 'PP3': 1, 'PP4': 1, 'PP5': 1,  # Supporting Pathogenic
    'BA1': -8,  # Stand Alone Benign
    'BS1': -4, 'BS2': -4, 'BS3': -4, 'BS4': -4,  # Strong Benign
    'BP1': -1, 'BP2': -1, 'BP3': -1, 'BP4': -1, 'BP5': -1, 'BP6': -1, 'BP7': -1  # Supporting Benign
}
```

### 6.2 Classification Algorithm
```python
def classify_variant(applied_criteria):
    """
    Final variant classification based on applied criteria
    """
    total_score = sum(EVIDENCE_WEIGHTS[criterion] for criterion in applied_criteria)
    
    # ACMG/AMP classification rules
    if total_score >= 8:
        return 'Pathogenic'
    elif 4 <= total_score <= 7:
        return 'Likely Pathogenic'
    elif -3 <= total_score <= 3:
        return 'Uncertain Significance'
    elif -7 <= total_score <= -4:
        return 'Likely Benign'
    else:  # total_score <= -8
        return 'Benign'
```

### 6.3 Confidence Scoring
```python
def calculate_confidence(applied_criteria, data_completeness):
    """
    Calculate classification confidence based on evidence strength and data quality
    """
    evidence_count = len(applied_criteria)
    strong_evidence = sum(1 for c in applied_criteria if c.startswith(('PVS', 'PS', 'BS', 'BA')))
    
    base_confidence = min(evidence_count / 3, 1.0)  # Normalized by typical evidence count
    strength_bonus = strong_evidence * 0.2  # Bonus for strong evidence
    data_penalty = (1 - data_completeness) * 0.3  # Penalty for missing data
    
    confidence = max(0, min(1, base_confidence + strength_bonus - data_penalty))
    
    if confidence >= 0.8:
        return 'High'
    elif confidence >= 0.6:
        return 'Medium'
    else:
        return 'Low'
```

---

## 7. DATA SOURCES AND INTEGRATION

### 7.1 Population Databases
- **gnomAD v3.1.2**: Primary population frequency source (763,406 exomes + genomes)
- **1000 Genomes**: Additional population frequency validation
- **ESP6500**: European and African American populations
- **ExAC**: Predecessor to gnomAD (retained for comparison)

### 7.2 Clinical Databases
- **ClinVar**: Clinical significance annotations and submissions
- **HGMD**: Human Gene Mutation Database (literature references)
- **OMIM**: Gene-disease associations and inheritance patterns

### 7.3 Functional Prediction Databases
- **dbNSFP v4.3**: Comprehensive functional prediction scores
- **Varsome API**: Aggregated variant annotation
- **Ensembl VEP**: Consequence prediction and annotation

---

## 8. QUALITY CONTROL AND VALIDATION

### 8.1 Input Validation
```python
def validate_variant_input(variant_data):
    """
    Comprehensive input validation with error handling
    """
    validators = {
        'chromosome': validate_chromosome,
        'position': validate_genomic_position,
        'ref_alt': validate_alleles,
        'frequencies': validate_frequency_range,
        'scores': validate_predictor_scores
    }
    
    for field, validator in validators.items():
        try:
            validator(variant_data[field])
        except ValidationError as e:
            raise InvalidInputError(f"Invalid {field}: {e}")
```

### 8.2 Cross-Validation
- **Internal consistency**: Checks between population frequencies and clinical significance
- **Predictor correlation**: Flags discordant computational predictions
- **Literature consistency**: Warns about conflicts between automated and manual evidence

### 8.3 Error Handling
```python
def handle_missing_data(criterion, data_availability):
    """
    Graceful handling of missing data with appropriate fallbacks
    """
    if data_availability < 0.5:  # Less than 50% data available
        return {
            'applies': False,
            'reason': 'insufficient_data',
            'recommendation': f'Manual review required for {criterion}'
        }
    else:
        return evaluate_with_available_data(criterion, data_availability)
```

---

## 9. TECHNICAL ARCHITECTURE

### 9.1 Core Classes
```python
class VariantData:
    """
    Structured data container for variant information
    """
    def __init__(self, variant_dict):
        self.basic_info = variant_dict.get('basic_info', {})
        self.population_data = variant_dict.get('population_data', {})
        self.insilico_data = variant_dict.get('insilico_data', {})
        self.genetic_data = variant_dict.get('genetic_data', {})
        self.functional_data = variant_dict.get('functional_data', {})

class EvidenceEvaluator:
    """
    Core evaluation engine for ACMG/AMP criteria
    """
    def __init__(self, test_mode=False, use_2023_guidelines=False):
        self.test_mode = test_mode
        self.use_2023_guidelines = use_2023_guidelines
        self.applied_criteria = {}
    
    def evaluate_all_criteria(self, variant_data):
        """
        Main evaluation pipeline
        """
        pathogenic_evidence = self._evaluate_pathogenic_criteria(variant_data)
        benign_evidence = self._evaluate_benign_criteria(variant_data)
        
        return self._consolidate_results(pathogenic_evidence, benign_evidence)
```

### 9.2 Modular Architecture
```
src/
├── core/
│   ├── evidence_evaluator.py    # Main evaluation logic
│   ├── acmg_classifier.py       # Classification algorithms
│   └── variant_data.py          # Data structures
├── utils/
│   ├── input_handler.py         # Data validation and processing
│   ├── report_generator.py      # Output formatting
│   └── validators.py            # Input validation functions
└── config/
    └── constants.py             # Thresholds and configuration
```

---

## 10. LIMITATIONS AND FUTURE DIRECTIONS

### 10.1 Current Limitations
1. **Manual Score Entry**: All predictor scores require manual input
2. **Literature Review**: Interactive criteria depend on user expertise
3. **Population Specificity**: Limited ethnicity-specific frequency data
4. **Functional Studies**: No automated assessment of experimental evidence

### 10.2 Statistical Considerations
1. **Multiple Testing**: No correction for multiple criteria evaluation
2. **Predictor Correlation**: Some predictors are not independent
3. **Population Structure**: Limited consideration of population stratification
4. **Penetrance Variability**: Fixed penetrance assumptions for frequency calculations

### 10.3 Future Enhancements
1. **Machine Learning Integration**: Automated literature mining
2. **Real-time Database Queries**: API integration for current data
3. **Ethnic-Specific Analysis**: Population-stratified frequency analysis
4. **Uncertainty Quantification**: Bayesian confidence intervals

---

## 11. REFERENCES AND IMPLEMENTATION STANDARDS

### 11.1 Primary Guidelines
1. Richards, S. et al. (2015). Standards and guidelines for the interpretation of sequence variants. *Genetics in Medicine* 17:405-424.
2. ClinGen Sequence Variant Interpretation Working Group. Updates to evidence criteria (2023).

### 11.2 Computational Methods
1. Ioannidis, N.M. et al. (2016). REVEL: An ensemble method for predicting the pathogenicity of rare missense variants. *American Journal of Human Genetics* 99:877-885.
2. Cheng, J. et al. (2023). Accurate proteome-wide missense variant effect prediction with AlphaMissense. *Science* 381:eadg7492.
3. Jaganathan, K. et al. (2019). Predicting splicing from primary sequence with deep learning. *Cell* 176:535-548.

### 11.3 Population Genetics
1. Karczewski, K.J. et al. (2020). The mutational constraint spectrum quantified from variation in 141,456 humans. *Nature* 581:434-443.
2. Lek, M. et al. (2016). Analysis of protein-coding genetic variation in 60,706 humans. *Nature* 536:285-291.

---

## 12. ACADEMIC IMPACT AND COMMUNITY CONTRIBUTIONS

### 12.1 Open Source Implementation
This implementation represents the first complete, open-source implementation of all 28 ACMG/AMP criteria with full technical transparency. The codebase is available for academic review, validation, and extension by the research community.

### 12.2 Clinical Validation Opportunities
The systematic implementation provides a standardized platform for:
- Multi-institutional validation studies
- Comparative analysis with clinical laboratory interpretations
- Benchmarking against expert panel classifications
- Training dataset development for machine learning approaches

### 12.3 Research Applications
**Computational Biology Research:**
- Ensemble predictor development and optimization
- Population-specific frequency threshold optimization
- Machine learning model training and validation
- Uncertainty quantification methodologies

**Clinical Genetics Research:**
- Large-scale variant reclassification studies
- Inter-laboratory concordance analysis
- Evidence criteria modification proposals
- Phenotype-genotype correlation studies

### 12.4 Educational Value
**For Graduate Students and Residents:**
- Complete reference implementation of ACMG/AMP guidelines
- Hands-on learning tool for variant interpretation principles
- Understanding of computational prediction integration
- Interactive learning for literature-based evidence evaluation

**For Clinical Laboratory Personnel:**
- Standardized workflow reference
- Quality control validation tool
- Training resource for new staff
- Benchmark for internal classification systems

### 12.5 Community Feedback and Validation
**Encouraged Contributions:**
- Clinical laboratory validation reports
- Population-specific threshold optimizations
- Additional computational predictor integrations
- Enhanced statistical methodologies
- Bug reports and feature requests

**Academic Collaboration Opportunities:**
- Multi-center validation studies
- Comparative effectiveness research
- Guideline development committee input
- Educational curriculum integration

### 12.6 Publication and Citation
**Academic Publications:**
This work builds upon and cites foundational research in variant interpretation, computational prediction, and population genetics. Users of this tool in academic research are encouraged to cite both this implementation and the underlying methodological publications referenced throughout this document.

**Reproducible Research:**
All algorithms, thresholds, and statistical methods are fully documented and implemented in open-source code, supporting reproducible research principles and enabling independent validation.

---

## 13. CONTACT AND COLLABORATION

### 13.1 Technical Support
For technical questions, bug reports, or implementation clarifications:
- **Email**: cansevilmiss@gmail.com
- **GitHub Issues**: https://github.com/Bilmem2/acmg-assessor/issues
- **LinkedIn**: https://linkedin.com/in/cansevilmiss

### 13.2 Academic Collaboration
We welcome collaboration with:
- Clinical genetics laboratories seeking validation studies
- Academic researchers developing variant interpretation methods
- Educational institutions incorporating variant interpretation training
- International consortiums standardizing variant classification

### 13.3 Data Sharing and Validation
**Validation Data Sharing:**
Institutions conducting validation studies are encouraged to share anonymized results to contribute to community knowledge and method improvement.

**Best Practices Documentation:**
We seek to collaborate with clinical experts to document best practices for interactive criteria evaluation and interpretation workflows.

---

**Document Version**: 1.0  
**Last Updated**: July 9, 2025  
**Total Implementation**: 28/28 ACMG/AMP Criteria (100%)  
**Repository**: https://github.com/Bilmem2/acmg-assessor  
**License**: MIT License (Open Source)  

**For technical questions, academic collaboration, or clinical validation partnerships, contact**: cansevilmiss@gmail.com

---

*This technical report is provided to support transparency, reproducibility, and academic validation of automated variant classification methodologies. We encourage peer review, independent validation, and collaborative improvement of these methods for the benefit of the clinical genetics community.*

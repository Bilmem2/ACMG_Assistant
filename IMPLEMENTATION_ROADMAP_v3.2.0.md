# ACMG Assistant Implementation Roadmap v3.2.0

## Executive Summary

Bu dokümantasyon, ACMG Variant Classification Assistant'ın test sonuçlarına dayalı sistematik geliştirme planını içermektedir. 30 gerçek varyant üzerinde yapılan kapsamlı testler sonucunda belirlenen eksiklikler ve iyileştirme alanları için detaylı implementasyon stratejileri sunulmuştur.

## Current Performance Analysis

### Test Results Summary
- **Total Variants Tested**: 30
- **Exact Match Accuracy**: 40.0% (12/30)
- **Near Match Accuracy**: 60.0% (18/30)
- **Expected Accuracy**: 43.3% (13/30)

### Performance by Variant Type
```
Stop gained:         100% (6/6)  - EXCELLENT
Frameshift variant:  100% (4/4)  - EXCELLENT
Missense variant:    12.5% (2/16) - NEEDS MAJOR IMPROVEMENT
Inframe deletion:    0% (0/3)     - NEEDS COMPLETE REWRITE
Synonymous variant:  0% (0/1)     - NEEDS IMPLEMENTATION
```

### Performance by ACMG Criteria
```
EXCELLENT (100%):    BA1, BS4, PM2, PP1, PP3, PS2, PS4
GOOD (75%):          PVS1
MODERATE (50%):      PS1
CRITICAL (0%):       PM1, PM3, PM4, PM5, PP2, PP4, PS3
```

## Phase 1: Critical Fixes (Priority 1)

### 1.1 Missense Variant Metascore Engine

**Problem**: Missense variants have only 12.5% accuracy
**Impact**: High - affects 16/30 test cases
**Timeline**: 2-3 weeks

#### Implementation Strategy

```python
# New file: src/core/missense_evaluator.py
class MissenseEvaluator:
    def __init__(self):
        self.domain_regions = self._load_domain_regions()
        self.conservation_scores = self._load_conservation_data()
        self.structural_data = self._load_structural_data()
    
    def evaluate_missense_variant(self, variant_data):
        """Comprehensive missense variant evaluation"""
        scores = {
            'conservation': self._calculate_conservation_score(variant_data),
            'structural': self._calculate_structural_impact(variant_data),
            'functional': self._calculate_functional_score(variant_data),
            'domain': self._calculate_domain_impact(variant_data),
            'population': self._calculate_population_context(variant_data)
        }
        
        return self._generate_composite_score(scores)
    
    def _calculate_conservation_score(self, variant_data):
        """PhyloP, PhastCons, GERP++ integration"""
        phylop = getattr(variant_data, 'phylop_score', None)
        phastcons = getattr(variant_data, 'phastcons_score', None)
        gerp = getattr(variant_data, 'gerp_score', None)
        
        # Weighted conservation scoring
        conservation_weight = 0.0
        if phylop is not None:
            conservation_weight += 0.4 * self._normalize_phylop(phylop)
        if phastcons is not None:
            conservation_weight += 0.3 * self._normalize_phastcons(phastcons)
        if gerp is not None:
            conservation_weight += 0.3 * self._normalize_gerp(gerp)
            
        return conservation_weight
```

#### Required Data Files
```
data/
├── domain_regions.json         # Protein domain annotations
├── conservation_scores.json    # PhyloP, PhastCons, GERP++ data
├── structural_impact.json      # 3D structure impact predictions
└── functional_studies.json     # Curated functional study results
```

### 1.2 Inframe Deletion Logic (PM4)

**Problem**: Inframe deletions have 0% accuracy
**Impact**: Medium - affects 3/30 test cases
**Timeline**: 1-2 weeks

#### Implementation Strategy

```python
# Addition to src/core/evidence_evaluator.py
class InframeAnalyzer:
    def __init__(self):
        self.critical_regions = self._load_critical_regions()
        self.repeat_regions = self._load_repeat_regions()
        self.domain_boundaries = self._load_domain_boundaries()
    
    def evaluate_inframe_deletion(self, variant_data):
        """Enhanced inframe deletion evaluation"""
        # Check if deletion affects critical regions
        if self._affects_critical_region(variant_data):
            return 'PM4'
        
        # Check if deletion affects functional domains
        if self._affects_functional_domain(variant_data):
            return 'PM4'
        
        # Check if deletion affects structural integrity
        if self._affects_structural_integrity(variant_data):
            return 'PM4'
        
        # Check if deletion in repeat region (likely benign)
        if self._in_repeat_region(variant_data):
            return 'BP3'
        
        return None
    
    def _affects_critical_region(self, variant_data):
        """Check if deletion affects known critical regions"""
        gene = variant_data.gene
        position = self._extract_position(variant_data.hgvs_c)
        
        if gene in self.critical_regions:
            for region in self.critical_regions[gene]:
                if region['start'] <= position <= region['end']:
                    return True
        return False
```

### 1.3 Gene-Specific Hotspot Regions (PM1)

**Problem**: PM1 criteria has 0% accuracy
**Impact**: Medium - affects hotspot detection
**Timeline**: 1-2 weeks

#### Implementation Strategy

```python
# New file: src/core/gene_specific_rules.py
class GeneSpecificRules:
    def __init__(self):
        self.hotspot_regions = self._load_hotspot_regions()
        self.gene_specific_thresholds = self._load_gene_thresholds()
    
    def evaluate_hotspot_region(self, variant_data):
        """Gene-specific hotspot region evaluation"""
        gene = variant_data.gene
        position = self._extract_position(variant_data.hgvs_c)
        
        if gene in self.hotspot_regions:
            for hotspot in self.hotspot_regions[gene]:
                if hotspot['start'] <= position <= hotspot['end']:
                    # Check if this hotspot is well-established
                    if hotspot['evidence_level'] >= 3:
                        return 'PM1'
        
        return None
    
    def _load_hotspot_regions(self):
        """Load curated hotspot regions from literature"""
        return {
            'KRAS': [
                {'start': 34, 'end': 38, 'evidence_level': 5, 'description': 'Codon 12-13'},
                {'start': 175, 'end': 183, 'evidence_level': 5, 'description': 'Codon 59-61'},
                {'start': 436, 'end': 444, 'evidence_level': 4, 'description': 'Codon 146'}
            ],
            'TP53': [
                {'start': 730, 'end': 750, 'evidence_level': 5, 'description': 'DNA binding domain'},
                {'start': 818, 'end': 825, 'evidence_level': 4, 'description': 'Codon 273'}
            ],
            'BRAF': [
                {'start': 1798, 'end': 1800, 'evidence_level': 5, 'description': 'V600E hotspot'}
            ]
        }
```

## Phase 2: Functional Integration (Priority 2)

### 2.1 Functional Studies Integration (PS3/BS3)

**Problem**: PS3 criteria has 0% accuracy
**Impact**: High - affects experimental evidence
**Timeline**: 3-4 weeks

#### Implementation Strategy

```python
# New file: src/core/functional_studies_evaluator.py
class FunctionalStudiesEvaluator:
    def __init__(self):
        self.pubmed_client = PubMedClient()
        self.functional_db = FunctionalStudiesDatabase()
        self.study_quality_evaluator = StudyQualityEvaluator()
    
    def evaluate_functional_evidence(self, variant_data):
        """Evaluate functional studies evidence"""
        studies = self._search_functional_studies(variant_data)
        
        if not studies:
            return None
        
        # Evaluate study quality and consistency
        quality_scores = []
        pathogenicity_votes = []
        
        for study in studies:
            quality = self._evaluate_study_quality(study)
            if quality >= 3:  # Minimum quality threshold
                quality_scores.append(quality)
                pathogenicity_votes.append(study.pathogenicity_conclusion)
        
        if not quality_scores:
            return None
        
        # Determine consensus
        return self._determine_functional_consensus(quality_scores, pathogenicity_votes)
    
    def _search_functional_studies(self, variant_data):
        """Search for functional studies in multiple databases"""
        studies = []
        
        # Search PubMed
        pubmed_studies = self.pubmed_client.search_variant_studies(
            gene=variant_data.gene,
            variant=variant_data.hgvs_c
        )
        studies.extend(pubmed_studies)
        
        # Search ClinVar functional annotations
        clinvar_studies = self._search_clinvar_functional(variant_data)
        studies.extend(clinvar_studies)
        
        # Search LOVD
        lovd_studies = self._search_lovd_functional(variant_data)
        studies.extend(lovd_studies)
        
        return studies
```

### 2.2 Phenotype Matching (PP4)

**Problem**: PP4 criteria has 0% accuracy
**Impact**: Medium - affects phenotype-genotype correlation
**Timeline**: 2-3 weeks

#### Implementation Strategy

```python
# New file: src/core/phenotype_matcher.py
class PhenotypeMatcher:
    def __init__(self):
        self.hpo_client = HPOClient()
        self.gene_phenotype_db = GenePhenotypeDatabase()
        self.phenotype_similarity_calculator = PhenotypeSimilarityCalculator()
    
    def evaluate_phenotype_match(self, variant_data, patient_phenotypes):
        """Evaluate phenotype-genotype correlation"""
        gene = variant_data.gene
        
        # Get known phenotypes for this gene
        known_phenotypes = self.gene_phenotype_db.get_gene_phenotypes(gene)
        
        if not known_phenotypes or not patient_phenotypes:
            return None
        
        # Calculate phenotype similarity
        similarity_scores = []
        for known_phenotype in known_phenotypes:
            for patient_phenotype in patient_phenotypes:
                similarity = self.phenotype_similarity_calculator.calculate_similarity(
                    known_phenotype, patient_phenotype
                )
                similarity_scores.append(similarity)
        
        max_similarity = max(similarity_scores) if similarity_scores else 0
        
        # Determine PP4 applicability
        if max_similarity >= 0.8:
            return 'PP4'
        elif max_similarity <= 0.3:
            return 'BP5'  # Phenotype not consistent
        
        return None
```

## Phase 3: Population Genetics Enhancement (Priority 3)

### 3.1 Ethnicity-Specific Frequency Analysis

**Problem**: Population frequency analysis lacks ethnic stratification
**Impact**: Medium - affects BA1/BS1/PM2 accuracy
**Timeline**: 2-3 weeks

#### Implementation Strategy

```python
# Enhancement to src/core/population_analyzer.py
class EthnicityAwarePopulationAnalyzer:
    def __init__(self):
        self.gnomad_client = GnomADClient()
        self.ethnicity_thresholds = self._load_ethnicity_thresholds()
    
    def analyze_population_frequency(self, variant_data, patient_ethnicity=None):
        """Ethnicity-aware population frequency analysis"""
        frequencies = self.gnomad_client.get_population_frequencies(variant_data)
        
        if patient_ethnicity:
            # Use ethnicity-specific thresholds
            relevant_freq = frequencies.get(patient_ethnicity, frequencies.get('ALL'))
            threshold = self.ethnicity_thresholds.get(patient_ethnicity, 0.01)
        else:
            # Use overall frequency with conservative threshold
            relevant_freq = frequencies.get('ALL')
            threshold = 0.005  # More conservative when ethnicity unknown
        
        # BA1 - Stand-alone benign
        if relevant_freq >= threshold:
            return 'BA1'
        
        # BS1 - Strong benign
        if relevant_freq >= threshold * 0.1:
            return 'BS1'
        
        # PM2 - Moderate pathogenic
        if relevant_freq <= 0.0001:
            return 'PM2'
        
        return None
```

## Phase 4: Algorithm Architecture Enhancement (Priority 4)

### 4.1 Evidence Weighting System

**Problem**: Current system lacks evidence weighting
**Impact**: High - affects overall classification accuracy
**Timeline**: 3-4 weeks

#### Implementation Strategy

```python
# Enhancement to src/core/acmg_classifier.py
class WeightedACMGClassifier:
    def __init__(self):
        self.evidence_weights = self._load_evidence_weights()
        self.conflict_resolver = ConflictResolver()
    
    def classify_variant_weighted(self, evidence_list):
        """Weighted ACMG classification with conflict resolution"""
        # Separate pathogenic and benign evidence
        pathogenic_evidence = [e for e in evidence_list if e.startswith(('PVS', 'PS', 'PM', 'PP'))]
        benign_evidence = [e for e in evidence_list if e.startswith(('BA', 'BS', 'BP'))]
        
        # Calculate weighted scores
        pathogenic_score = self._calculate_weighted_score(pathogenic_evidence)
        benign_score = self._calculate_weighted_score(benign_evidence)
        
        # Resolve conflicts
        if pathogenic_score > 0 and benign_score > 0:
            return self.conflict_resolver.resolve_conflict(
                pathogenic_evidence, benign_evidence, pathogenic_score, benign_score
            )
        
        # Standard classification with weights
        return self._determine_classification(pathogenic_score, benign_score)
    
    def _calculate_weighted_score(self, evidence_list):
        """Calculate weighted evidence score"""
        score = 0
        for evidence in evidence_list:
            base_weight = self.evidence_weights.get(evidence[:3], 1.0)  # PVS, PS, PM, etc.
            confidence = self._get_evidence_confidence(evidence)
            score += base_weight * confidence
        return score
```

### 4.2 Gene-Specific Rules Engine

**Problem**: Lacks gene-specific classification rules
**Impact**: High - affects accuracy for specific genes
**Timeline**: 4-5 weeks

#### Implementation Strategy

```python
# New file: src/core/gene_rules_engine.py
class GeneRulesEngine:
    def __init__(self):
        self.gene_rules = self._load_gene_specific_rules()
        self.inheritance_patterns = self._load_inheritance_patterns()
    
    def apply_gene_specific_rules(self, variant_data, preliminary_evidence):
        """Apply gene-specific classification rules"""
        gene = variant_data.gene
        
        if gene not in self.gene_rules:
            return preliminary_evidence
        
        rules = self.gene_rules[gene]
        modified_evidence = preliminary_evidence.copy()
        
        for rule in rules:
            if self._rule_conditions_met(rule, variant_data):
                modified_evidence = self._apply_rule_modifications(rule, modified_evidence)
        
        return modified_evidence
    
    def _load_gene_specific_rules(self):
        """Load gene-specific classification rules"""
        return {
            'BRCA1': [
                {
                    'condition': 'missense_in_ring_domain',
                    'action': 'upgrade_PM1_to_PM2',
                    'evidence_level': 4
                },
                {
                    'condition': 'nonsense_before_exon_20',
                    'action': 'apply_PVS1',
                    'evidence_level': 5
                }
            ],
            'BRCA2': [
                {
                    'condition': 'missense_in_dna_binding_domain',
                    'action': 'apply_PM1',
                    'evidence_level': 4
                }
            ],
            'TP53': [
                {
                    'condition': 'missense_in_dna_binding_domain',
                    'action': 'upgrade_PM1_to_PS1',
                    'evidence_level': 5
                }
            ]
        }
```

## Phase 5: Performance Optimization (Priority 5)

### 5.1 Caching and Performance

**Timeline**: 1-2 weeks

#### Implementation Strategy

```python
# Enhancement to src/utils/api_client.py
class CachedAPIClient:
    def __init__(self):
        self.cache = APICache()
        self.rate_limiter = RateLimiter()
    
    @cached(ttl=3600)  # Cache for 1 hour
    def get_variant_annotations(self, variant_key):
        """Cached variant annotation retrieval"""
        return self._fetch_variant_annotations(variant_key)
    
    def batch_annotate_variants(self, variant_list):
        """Batch annotation for multiple variants"""
        cached_results = {}
        uncached_variants = []
        
        for variant in variant_list:
            cached_result = self.cache.get(variant.key)
            if cached_result:
                cached_results[variant.key] = cached_result
            else:
                uncached_variants.append(variant)
        
        # Batch process uncached variants
        if uncached_variants:
            batch_results = self._batch_fetch_annotations(uncached_variants)
            for variant, result in batch_results.items():
                self.cache.set(variant.key, result)
                cached_results[variant.key] = result
        
        return cached_results
```

## Implementation Timeline

### Week 1-2: Phase 1 Setup
- [ ] Create new file structure for enhanced modules
- [ ] Implement basic MissenseEvaluator class
- [ ] Set up InframeAnalyzer framework
- [ ] Create GeneSpecificRules foundation

### Week 3-4: Phase 1 Core Implementation
- [ ] Complete missense variant metascore engine
- [ ] Implement inframe deletion logic
- [ ] Add hotspot region detection
- [ ] Create comprehensive test cases

### Week 5-6: Phase 2 Functional Integration
- [ ] Implement FunctionalStudiesEvaluator
- [ ] Add PhenotypeMatcher
- [ ] Integrate with external APIs
- [ ] Add quality control mechanisms

### Week 7-8: Phase 3 Population Genetics
- [ ] Enhance population frequency analysis
- [ ] Add ethnicity-specific thresholds
- [ ] Implement stratified analysis
- [ ] Add population database integration

### Week 9-12: Phase 4 Architecture Enhancement
- [ ] Implement evidence weighting system
- [ ] Create gene-specific rules engine
- [ ] Add conflict resolution mechanisms
- [ ] Implement comprehensive logging

### Week 13-14: Phase 5 Optimization
- [ ] Add caching mechanisms
- [ ] Implement batch processing
- [ ] Optimize performance
- [ ] Add monitoring and metrics

## Data Requirements

### External Databases
- **gnomAD**: Population frequency data
- **ClinVar**: Clinical significance annotations
- **LOVD**: Locus-specific databases
- **PubMed**: Literature search API
- **HPO**: Human Phenotype Ontology
- **UniProt**: Protein domain annotations

### Internal Data Files
```
data/
├── gene_rules/
│   ├── brca1_rules.json
│   ├── brca2_rules.json
│   └── tp53_rules.json
├── domain_annotations/
│   ├── protein_domains.json
│   └── critical_regions.json
├── population_data/
│   ├── ethnicity_thresholds.json
│   └── population_frequencies.json
└── functional_studies/
    ├── curated_studies.json
    └── study_quality_metrics.json
```

## Success Metrics

### Target Improvements
- **Overall Accuracy**: 40% → 75% (87.5% improvement)
- **Missense Accuracy**: 12.5% → 60% (380% improvement)
- **Inframe Deletion**: 0% → 80% (new capability)
- **Gene-Specific Criteria**: 0% → 70% (new capability)

### Testing Strategy
- Expand test suite to 100+ variants
- Add edge case testing
- Implement regression testing
- Add performance benchmarking

## Risk Mitigation

### Technical Risks
1. **API Rate Limiting**: Implement caching and batch processing
2. **Data Quality**: Add data validation and quality checks
3. **Performance**: Optimize algorithms and add monitoring
4. **Complexity**: Maintain modular architecture

### Implementation Risks
1. **Timeline**: Prioritize critical fixes first
2. **Resource**: Implement incrementally
3. **Testing**: Comprehensive test coverage
4. **Documentation**: Maintain updated documentation

## Conclusion

This roadmap provides a systematic approach to improving ACMG Assistant accuracy from 40% to 75%+ through targeted algorithmic enhancements. The phased approach ensures critical fixes are implemented first while building a robust foundation for future improvements.

The implementation will be tracked through detailed logging, metrics, and regular testing to ensure progress towards accuracy targets.

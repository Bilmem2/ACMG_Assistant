"""
Tests for Multi-Source Predictor and Population API Clients
============================================================

Tests the PredictorAPIClient and PopulationAPIClient classes with
mocked API responses to verify:
- Multi-source data fetching with priority fallback
- Graceful degradation when sources are unavailable
- Proper handling of typed PredictorScore and PopulationStats objects
- Integration with MissenseEvaluator and PopulationAnalyzer
- Weight renormalization for partial data

Author: Can Sevilmiş
License: MIT License
"""

import pytest
from unittest.mock import Mock, patch, MagicMock
from dataclasses import asdict

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from config.predictors import (
    PredictorScore,
    PopulationStats,
    PREDICTOR_SOURCE_PRIORITY,
    POPULATION_SOURCE_PRIORITY,
    INVERTED_PREDICTORS,
    INSILICO_WEIGHTS,
    MIN_PREDICTORS_FOR_COMPOSITE,
)
from utils.predictor_api_client import PredictorAPIClient, PopulationAPIClient
from core.missense_evaluator import MissenseEvaluator
from core.population_analyzer import GnomADClient, EthnicityAwarePopulationAnalyzer
from core.variant_data import VariantData


# =============================================================================
# PredictorScore Dataclass Tests
# =============================================================================

class TestPredictorScore:
    """Test PredictorScore dataclass functionality."""
    
    def test_predictor_score_creation(self):
        """Test creating a PredictorScore with all fields."""
        score = PredictorScore(
            predictor='revel',
            value=0.85,
            source='dbNSFP',
            version='4.5',
            raw={'some': 'data'},
            is_inverted=False
        )
        
        assert score.predictor == 'revel'
        assert score.value == 0.85
        assert score.source == 'dbNSFP'
        assert score.version == '4.5'
        assert score.is_inverted is False
    
    def test_predictor_score_is_available(self):
        """Test is_available() method."""
        score_with_value = PredictorScore(predictor='revel', value=0.5)
        score_without_value = PredictorScore(predictor='revel', value=None)
        
        assert score_with_value.is_available() is True
        assert score_without_value.is_available() is False
    
    def test_predictor_score_normalized_value_normal(self):
        """Test get_normalized_value() for non-inverted predictors."""
        score = PredictorScore(predictor='revel', value=0.8, is_inverted=False)
        assert score.get_normalized_value() == 0.8
    
    def test_predictor_score_normalized_value_sift(self):
        """Test get_normalized_value() for SIFT (inverted)."""
        # SIFT 0.01 = deleterious → should normalize to ~0.99
        score = PredictorScore(predictor='sift', value=0.01, is_inverted=True)
        assert score.get_normalized_value() == pytest.approx(0.99, abs=0.01)
    
    def test_predictor_score_normalized_value_fathmm(self):
        """Test get_normalized_value() for FATHMM (inverted, negative = pathogenic)."""
        # FATHMM -3.0 = deleterious → should normalize to high value
        score = PredictorScore(predictor='fathmm', value=-3.0, is_inverted=True)
        normalized = score.get_normalized_value()
        assert 0.7 < normalized < 0.9  # Should be high (pathogenic)
    
    def test_predictor_score_none_value_normalized(self):
        """Test get_normalized_value() returns None when value is None."""
        score = PredictorScore(predictor='revel', value=None)
        assert score.get_normalized_value() is None


# =============================================================================
# PopulationStats Dataclass Tests
# =============================================================================

class TestPopulationStats:
    """Test PopulationStats dataclass functionality."""
    
    def test_population_stats_creation(self):
        """Test creating PopulationStats with all fields."""
        stats = PopulationStats(
            population='gnomad_v4',
            af=0.00015,
            an=152000,
            ac=23,
            homozygote_count=0,
            subpop={'nfe': {'af': 0.0002, 'ac': 15, 'an': 75000}},
            popmax_af=0.0003,
            popmax_population='nfe',
            source='gnomAD_GraphQL',
            version='v4.0'
        )
        
        assert stats.population == 'gnomad_v4'
        assert stats.af == 0.00015
        assert stats.an == 152000
        assert stats.ac == 23
    
    def test_population_stats_is_available(self):
        """Test is_available() method."""
        stats_with_af = PopulationStats(population='gnomad_v4', af=0.001)
        stats_with_ac = PopulationStats(population='gnomad_v4', ac=5)
        stats_empty = PopulationStats(population='gnomad_v4')
        
        assert stats_with_af.is_available() is True
        assert stats_with_ac.is_available() is True
        assert stats_empty.is_available() is False
    
    def test_population_stats_is_absent(self):
        """Test is_absent() method."""
        stats_absent_af = PopulationStats(population='gnomad_v4', af=0.0)
        stats_absent_ac = PopulationStats(population='gnomad_v4', ac=0)
        stats_present = PopulationStats(population='gnomad_v4', af=0.001)
        stats_no_data = PopulationStats(population='gnomad_v4')
        
        assert stats_absent_af.is_absent() is True
        assert stats_absent_ac.is_absent() is True
        assert stats_present.is_absent() is False
        assert stats_no_data.is_absent() is True
    
    def test_population_stats_get_max_subpop_af(self):
        """Test get_max_subpop_af() method."""
        stats = PopulationStats(
            population='gnomad_v4',
            af=0.0001,
            popmax_af=0.0003,
            subpop={
                'nfe': {'af': 0.0002},
                'afr': {'af': 0.0003}
            }
        )
        
        # Should return popmax_af when available
        assert stats.get_max_subpop_af() == 0.0003


# =============================================================================
# PredictorAPIClient Tests with Mocked Responses
# =============================================================================

class TestPredictorAPIClient:
    """Test PredictorAPIClient with mocked API responses."""
    
    def test_client_initialization(self):
        """Test client initialization with default parameters."""
        client = PredictorAPIClient(api_enabled=True, timeout=10)
        assert client.api_enabled is True
        assert client.timeout == 10
        assert client.test_mode is False
    
    def test_test_mode_returns_mock_scores(self):
        """Test that test_mode returns mock predictor scores."""
        client = PredictorAPIClient(api_enabled=True, test_mode=True)
        
        scores = client.get_predictor_scores(
            chrom='17', pos=43092919, ref='G', alt='A'
        )
        
        assert 'revel' in scores
        assert 'cadd_phred' in scores
        assert scores['revel'].value == 0.75
        assert scores['cadd_phred'].value == 25.0
        assert scores['revel'].source == 'mock'
    
    def test_api_disabled_returns_empty_scores(self):
        """Test that api_enabled=False returns empty scores."""
        client = PredictorAPIClient(api_enabled=False)
        
        scores = client.get_predictor_scores(
            chrom='17', pos=43092919, ref='G', alt='A'
        )
        
        # Should return PredictorScore objects with value=None
        assert all(s.value is None for s in scores.values())
    
    @patch('utils.predictor_api_client.requests.get')
    def test_myvariant_api_success(self, mock_get):
        """Test successful myvariant.info API call."""
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.json.return_value = {
            'dbnsfp': {
                'revel': {'score': 0.85},
                'cadd': {'phred': 28.5},
                'alphamissense': {'am_pathogenicity': 0.9},
                'sift': {'score': 0.02},
                'polyphen2': {'hdiv': {'score': 0.95}},
            }
        }
        mock_get.return_value = mock_response
        
        client = PredictorAPIClient(api_enabled=True)
        scores = client.get_predictor_scores(
            chrom='17', pos=43092919, ref='G', alt='A'
        )
        
        assert scores['revel'].value == 0.85
        assert scores['revel'].source == 'dbNSFP'
        assert scores['cadd_phred'].value == 28.5
        assert scores['alphamissense'].value == 0.9
        assert scores['sift'].value == 0.02
        assert scores['sift'].is_inverted is True
    
    @patch('utils.predictor_api_client.requests.get')
    def test_myvariant_api_404_variant_not_found(self, mock_get):
        """Test handling of variant not found in myvariant.info."""
        mock_response = Mock()
        mock_response.status_code = 404
        mock_get.return_value = mock_response
        
        client = PredictorAPIClient(api_enabled=True)
        scores = client.get_predictor_scores(
            chrom='17', pos=43092919, ref='G', alt='A'
        )
        
        # Should return scores with value=None
        assert scores['revel'].value is None
    
    @patch('utils.predictor_api_client.requests.get')
    def test_myvariant_api_timeout(self, mock_get):
        """Test handling of API timeout."""
        import requests
        mock_get.side_effect = requests.exceptions.Timeout()
        
        client = PredictorAPIClient(api_enabled=True)
        scores = client.get_predictor_scores(
            chrom='17', pos=43092919, ref='G', alt='A'
        )
        
        # Should handle timeout gracefully, return empty scores
        assert scores['revel'].value is None
    
    def test_get_available_predictor_count(self):
        """Test counting available predictors."""
        client = PredictorAPIClient(api_enabled=True, test_mode=True)
        scores = client.get_predictor_scores(
            chrom='17', pos=43092919, ref='G', alt='A'
        )
        
        count = client.get_available_predictor_count(scores)
        assert count >= 5  # Mock should have multiple scores
    
    def test_get_weighted_predictors(self):
        """Test getting weighted predictors sorted by weight."""
        client = PredictorAPIClient(api_enabled=True, test_mode=True)
        scores = client.get_predictor_scores(
            chrom='17', pos=43092919, ref='G', alt='A'
        )
        
        weighted = client.get_weighted_predictors(scores)
        
        # Should be sorted by weight descending
        weights = [w for _, _, w in weighted]
        assert weights == sorted(weights, reverse=True)
        
        # First should be REVEL (highest weight)
        assert weighted[0][0] == 'revel'


# =============================================================================
# PopulationAPIClient Tests with Mocked Responses
# =============================================================================

class TestPopulationAPIClient:
    """Test PopulationAPIClient with mocked API responses."""
    
    def test_client_initialization(self):
        """Test client initialization."""
        client = PopulationAPIClient(api_enabled=True, timeout=15)
        assert client.api_enabled is True
        assert client.timeout == 15
    
    def test_test_mode_returns_mock_stats(self):
        """Test that test_mode returns mock population stats."""
        client = PopulationAPIClient(api_enabled=True, test_mode=True)
        
        stats = client.get_population_stats(
            chrom='17', pos=43092919, ref='G', alt='A'
        )
        
        assert 'gnomad_v4' in stats
        assert stats['gnomad_v4'].af == 0.00001
        assert stats['gnomad_v4'].source == 'mock'
    
    def test_api_disabled_returns_empty(self):
        """Test that api_enabled=False returns empty stats."""
        client = PopulationAPIClient(api_enabled=False)
        
        stats = client.get_population_stats(
            chrom='17', pos=43092919, ref='G', alt='A'
        )
        
        assert stats == {}
    
    def test_get_max_frequency(self):
        """Test get_max_frequency() helper."""
        client = PopulationAPIClient(api_enabled=True, test_mode=True)
        stats = client.get_population_stats(
            chrom='17', pos=43092919, ref='G', alt='A'
        )
        
        max_af = client.get_max_frequency(stats)
        assert max_af is not None
        assert max_af >= 0
    
    def test_is_absent_from_all(self):
        """Test is_absent_from_all() helper."""
        client = PopulationAPIClient(api_enabled=True)
        
        # Empty stats = absent
        assert client.is_absent_from_all({}) is True
        
        # Stats with AF=0 = absent
        absent_stats = {
            'gnomad_v4': PopulationStats(population='gnomad_v4', af=0.0)
        }
        assert client.is_absent_from_all(absent_stats) is True
        
        # Stats with AF>0 = not absent
        present_stats = {
            'gnomad_v4': PopulationStats(population='gnomad_v4', af=0.001)
        }
        assert client.is_absent_from_all(present_stats) is False


# =============================================================================
# MissenseEvaluator Integration Tests
# =============================================================================

class TestMissenseEvaluatorIntegration:
    """Test MissenseEvaluator with typed PredictorScore objects."""
    
    def test_calculate_functional_score_with_typed_scores(self):
        """Test functional score calculation using typed predictor_scores."""
        evaluator = MissenseEvaluator()
        
        # Create VariantData with typed predictor_scores
        variant_data = VariantData()
        variant_data.predictor_scores = {
            'revel': PredictorScore(predictor='revel', value=0.85, source='dbNSFP'),
            'cadd_phred': PredictorScore(predictor='cadd_phred', value=28.0, source='dbNSFP'),
            'alphamissense': PredictorScore(predictor='alphamissense', value=0.9, source='dbNSFP'),
            'sift': PredictorScore(predictor='sift', value=0.01, is_inverted=True, source='dbNSFP'),
            'polyphen2': PredictorScore(predictor='polyphen2', value=0.95, source='dbNSFP'),
        }
        
        score = evaluator._calculate_functional_score(variant_data)
        
        # High pathogenic scores should result in high functional score
        assert score > 0.7
    
    def test_calculate_functional_score_fallback_to_insilico_data(self):
        """Test fallback to legacy insilico_data when predictor_scores is None."""
        evaluator = MissenseEvaluator()
        
        # Create VariantData with legacy insilico_data
        variant_data = VariantData()
        variant_data.insilico_data = {
            'revel': 0.85,
            'cadd_phred': 28.0,
            'sift': 0.01,
            'polyphen2': 0.95,
        }
        
        score = evaluator._calculate_functional_score(variant_data)
        
        # Should still calculate score from insilico_data
        assert score > 0.7
    
    def test_calculate_functional_score_partial_data(self):
        """Test functional score with only partial predictor data."""
        evaluator = MissenseEvaluator()
        
        # Create VariantData with only REVEL and CADD
        variant_data = VariantData()
        variant_data.predictor_scores = {
            'revel': PredictorScore(predictor='revel', value=0.80, source='dbNSFP'),
            'cadd_phred': PredictorScore(predictor='cadd_phred', value=25.0, source='dbNSFP'),
            # Other predictors have value=None
            'sift': PredictorScore(predictor='sift', value=None),
            'alphamissense': PredictorScore(predictor='alphamissense', value=None),
        }
        
        score = evaluator._calculate_functional_score(variant_data)
        
        # Should still produce a valid score with weight renormalization
        assert 0 <= score <= 1
    
    def test_calculate_functional_score_no_data(self):
        """Test functional score returns neutral (0.5) when no data available."""
        evaluator = MissenseEvaluator()
        
        variant_data = VariantData()
        # No predictor_scores and no insilico_data
        
        score = evaluator._calculate_functional_score(variant_data)
        
        # Should return neutral score
        assert score == 0.5
    
    def test_calculate_population_context_with_typed_stats(self):
        """Test population context using typed population_stats."""
        evaluator = MissenseEvaluator()
        
        # Rare variant
        variant_data = VariantData()
        variant_data.population_stats = {
            'gnomad_v4': PopulationStats(population='gnomad_v4', af=0.00001)
        }
        
        score = evaluator._calculate_population_context(variant_data)
        
        # Rare variant should have high population score (more likely pathogenic)
        assert score > 0.7
    
    def test_calculate_population_context_common_variant(self):
        """Test population context for common variant."""
        evaluator = MissenseEvaluator()
        
        # Common variant (>5%)
        variant_data = VariantData()
        variant_data.population_stats = {
            'gnomad_v4': PopulationStats(population='gnomad_v4', af=0.08)
        }
        
        score = evaluator._calculate_population_context(variant_data)
        
        # Common variant should have low population score (likely benign)
        assert score < 0.3
    
    def test_evaluate_missense_variant_full_pipeline(self):
        """Test full missense variant evaluation with typed data."""
        evaluator = MissenseEvaluator()
        
        variant_data = VariantData()
        variant_data.predictor_scores = {
            'revel': PredictorScore(predictor='revel', value=0.88, source='dbNSFP'),
            'cadd_phred': PredictorScore(predictor='cadd_phred', value=30.0, source='dbNSFP'),
            'alphamissense': PredictorScore(predictor='alphamissense', value=0.92, source='dbNSFP'),
            'sift': PredictorScore(predictor='sift', value=0.001, is_inverted=True, source='dbNSFP'),
            'polyphen2': PredictorScore(predictor='polyphen2', value=0.998, source='dbNSFP'),
        }
        variant_data.population_stats = {
            'gnomad_v4': PopulationStats(population='gnomad_v4', af=0.0)
        }
        variant_data.functional_data = {
            'in_functional_domain': True,
            'in_hotspot': False
        }
        
        result = evaluator.evaluate_missense_variant(variant_data)
        
        # Should have pathogenic direction with high scores
        # Using 0.65 threshold since weights average out to lower composite
        assert result['composite_score'] > 0.65
        assert 'PP3' in result['evidence_category']
        assert result['direction'] == 'pathogenic'


# =============================================================================
# PopulationAnalyzer Integration Tests
# =============================================================================

class TestPopulationAnalyzerIntegration:
    """Test PopulationAnalyzer with typed PopulationStats."""
    
    def test_gnomad_client_with_typed_stats(self):
        """Test GnomADClient uses pre-fetched population_stats."""
        client = GnomADClient()
        
        # Create VariantData with typed population_stats
        variant_data = VariantData()
        variant_data.population_stats = {
            'gnomad_v4': PopulationStats(
                population='gnomad_v4',
                af=0.0001,
                subpop={
                    'NFE': {'af': 0.00012},
                    'AFR': {'af': 0.00008}
                },
                popmax_af=0.00015
            )
        }
        
        frequencies = client.get_population_frequencies(variant_data)
        
        assert 'ALL' in frequencies
        assert frequencies['ALL'] == 0.0001
        assert 'POPMAX' in frequencies
        assert frequencies['POPMAX'] == 0.00015
    
    def test_gnomad_client_fallback_to_population_data(self):
        """Test GnomADClient falls back to legacy population_data."""
        client = GnomADClient()
        
        # Create VariantData with legacy population_data
        variant_data = VariantData()
        variant_data.population_data = {
            'gnomad_af': 0.0002,
            'gnomad_af_popmax': 0.0003
        }
        
        frequencies = client.get_population_frequencies(variant_data)
        
        assert 'ALL' in frequencies
        assert frequencies['ALL'] == 0.0002
        assert frequencies['POPMAX'] == 0.0003
    
    def test_ethnicity_aware_analyzer_with_typed_stats(self):
        """Test EthnicityAwarePopulationAnalyzer with typed population_stats."""
        analyzer = EthnicityAwarePopulationAnalyzer()
        
        # Rare variant (should trigger PM2)
        variant_data = VariantData()
        variant_data.population_stats = {
            'gnomad_v4': PopulationStats(
                population='gnomad_v4',
                af=0.00001
            )
        }
        
        result = analyzer.analyze_population_frequency(variant_data)
        assert result == 'PM2'
    
    def test_ethnicity_aware_analyzer_common_variant(self):
        """Test EthnicityAwarePopulationAnalyzer with common variant."""
        analyzer = EthnicityAwarePopulationAnalyzer()
        
        # Common variant (should trigger BA1 or BS1)
        variant_data = VariantData()
        variant_data.population_stats = {
            'gnomad_v4': PopulationStats(
                population='gnomad_v4',
                af=0.02  # 2%
            )
        }
        
        result = analyzer.analyze_population_frequency(variant_data)
        assert result in ('BA1', 'BS1')


# =============================================================================
# VariantData Integration Tests
# =============================================================================

class TestVariantDataIntegration:
    """Test VariantData with new typed fields."""
    
    def test_variant_data_with_predictor_scores(self):
        """Test VariantData accepts predictor_scores field."""
        variant_data = VariantData()
        variant_data.predictor_scores = {
            'revel': PredictorScore(predictor='revel', value=0.8)
        }
        
        assert variant_data.has_predictor_scores() is True
        assert variant_data.get_predictor_score('revel').value == 0.8
    
    def test_variant_data_get_available_predictors(self):
        """Test get_available_predictors() method."""
        variant_data = VariantData()
        variant_data.predictor_scores = {
            'revel': PredictorScore(predictor='revel', value=0.8),
            'cadd_phred': PredictorScore(predictor='cadd_phred', value=None),
            'sift': PredictorScore(predictor='sift', value=0.02),
        }
        
        available = variant_data.get_available_predictors()
        
        assert 'revel' in available
        assert 'sift' in available
        assert 'cadd_phred' not in available
    
    def test_variant_data_with_population_stats(self):
        """Test VariantData accepts population_stats field."""
        variant_data = VariantData()
        variant_data.population_stats = {
            'gnomad_v4': PopulationStats(population='gnomad_v4', af=0.0001)
        }
        
        assert variant_data.has_population_stats() is True
        assert variant_data.get_population_stat('gnomad_v4').af == 0.0001
    
    def test_variant_data_get_max_population_af(self):
        """Test get_max_population_af() with multiple sources."""
        variant_data = VariantData()
        variant_data.population_stats = {
            'gnomad_v4': PopulationStats(population='gnomad_v4', af=0.0001),
            'exac': PopulationStats(population='exac', af=0.0002),
        }
        
        max_af = variant_data.get_max_population_af()
        
        assert max_af == 0.0002  # Maximum across sources
    
    def test_variant_data_is_absent_from_populations(self):
        """Test is_absent_from_populations() method."""
        # Absent variant
        variant_data = VariantData()
        variant_data.population_stats = {
            'gnomad_v4': PopulationStats(population='gnomad_v4', af=0.0)
        }
        
        assert variant_data.is_absent_from_populations() is True
        
        # Present variant
        variant_data.population_stats = {
            'gnomad_v4': PopulationStats(population='gnomad_v4', af=0.001)
        }
        
        assert variant_data.is_absent_from_populations() is False
    
    def test_variant_data_to_dict_with_typed_fields(self):
        """Test to_dict() serializes typed fields correctly."""
        variant_data = VariantData()
        variant_data.basic_info = {'gene': 'BRCA1'}
        variant_data.predictor_scores = {
            'revel': PredictorScore(predictor='revel', value=0.8, source='dbNSFP')
        }
        variant_data.population_stats = {
            'gnomad_v4': PopulationStats(population='gnomad_v4', af=0.0001)
        }
        
        data_dict = variant_data.to_dict()
        
        assert 'predictor_scores' in data_dict
        assert 'population_stats' in data_dict
        assert data_dict['predictor_scores']['revel']['value'] == 0.8


# =============================================================================
# Run Tests
# =============================================================================

if __name__ == '__main__':
    pytest.main([__file__, '-v', '--tb=short'])

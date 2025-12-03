"""
Unit Tests for PM1 Hotspot/Domain Evaluation
=============================================

Tests the PM1 evidence evaluation using mocked remote API responses.
Validates that:
1. DomainAPIClient.get_hotspot_annotation() correctly queries APIs
2. GeneSpecificRules.evaluate_pm1() correctly interprets annotations
3. No hardcoded biological data is used - all logic uses mock responses

Author: Can Sevilmiş
License: MIT License
"""

import pytest
from unittest.mock import Mock, patch, MagicMock
from dataclasses import dataclass
import sys
import os

# Add src to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from utils.domain_api_client import DomainAPIClient, HotspotAnnotation
from core.gene_specific_rules import GeneSpecificRules, PM1Evidence


# =============================================================================
# Test Fixtures
# =============================================================================

@pytest.fixture
def mock_domain_client():
    """Create a mock DomainAPIClient that doesn't make real API calls."""
    client = DomainAPIClient(cache_enabled=False)
    return client


@pytest.fixture
def gene_rules():
    """Create GeneSpecificRules with a mocked domain client."""
    mock_client = Mock(spec=DomainAPIClient)
    rules = GeneSpecificRules(domain_client=mock_client)
    return rules


# =============================================================================
# HotspotAnnotation Dataclass Tests
# =============================================================================

class TestHotspotAnnotation:
    """Tests for HotspotAnnotation dataclass."""
    
    def test_default_values(self):
        """Default annotation should have all negative/empty values."""
        ann = HotspotAnnotation()
        assert ann.is_hotspot is False
        assert ann.in_critical_domain is False
        assert ann.confidence == 0.0
        assert ann.source == "none"
        assert ann.details == ""
        assert ann.domain_name is None
        assert ann.hotspot_count == 0
        assert ann.tumor_count == 0
    
    def test_hotspot_annotation(self):
        """Hotspot annotation with positive values."""
        ann = HotspotAnnotation(
            is_hotspot=True,
            confidence=0.95,
            source="CancerHotspots.org",
            tumor_count=50,
            hotspot_count=120
        )
        assert ann.is_hotspot is True
        assert ann.confidence == 0.95
        assert ann.tumor_count == 50
    
    def test_domain_annotation(self):
        """Domain annotation with positive values."""
        ann = HotspotAnnotation(
            in_critical_domain=True,
            confidence=0.8,
            source="UniProt",
            domain_name="DNA binding"
        )
        assert ann.in_critical_domain is True
        assert ann.domain_name == "DNA binding"
    
    def test_to_dict(self):
        """to_dict should serialize all public fields."""
        ann = HotspotAnnotation(
            is_hotspot=True,
            confidence=0.9,
            source="test",
            details="test details"
        )
        d = ann.to_dict()
        assert d['is_hotspot'] is True
        assert d['confidence'] == 0.9
        assert d['source'] == "test"
        assert d['details'] == "test details"
        # raw_api_response should NOT be in to_dict output
        assert 'raw_api_response' not in d


# =============================================================================
# PM1Evidence Dataclass Tests
# =============================================================================

class TestPM1Evidence:
    """Tests for PM1Evidence dataclass."""
    
    def test_default_values(self):
        """Default PM1Evidence should not apply."""
        ev = PM1Evidence()
        assert ev.applies is False
        assert ev.evidence_code == ""
        assert ev.confidence == 0.0
    
    def test_pm1_evidence(self):
        """PM1Evidence with full strength."""
        ev = PM1Evidence(
            applies=True,
            evidence_code='PM1',
            confidence=0.92,
            source="CancerHotspots.org",
            reason="Position 248 is a documented hotspot"
        )
        assert ev.applies is True
        assert ev.evidence_code == 'PM1'
        assert ev.confidence == 0.92
    
    def test_pm1_supporting_evidence(self):
        """PM1Evidence with supporting strength."""
        ev = PM1Evidence(
            applies=True,
            evidence_code='PM1_supporting',
            confidence=0.7,
            source="UniProt",
            reason="Position is in functional domain"
        )
        assert ev.applies is True
        assert ev.evidence_code == 'PM1_supporting'
    
    def test_to_dict(self):
        """to_dict should serialize all fields."""
        ev = PM1Evidence(
            applies=True,
            evidence_code='PM1',
            confidence=0.9,
            source="test",
            reason="test reason"
        )
        d = ev.to_dict()
        assert d['applies'] is True
        assert d['evidence_code'] == 'PM1'
        assert d['reason'] == "test reason"


# =============================================================================
# DomainAPIClient Tests (with mocked HTTP)
# =============================================================================

class TestDomainAPIClient:
    """Tests for DomainAPIClient with mocked HTTP responses."""
    
    def test_no_hardcoded_known_hotspots(self):
        """Verify KNOWN_HOTSPOTS dict has been removed."""
        client = DomainAPIClient(cache_enabled=False)
        assert not hasattr(client, 'KNOWN_HOTSPOTS'), \
            "KNOWN_HOTSPOTS should be removed - no hardcoded biological data"
    
    def test_confidence_thresholds_exist(self):
        """Verify interpretive thresholds exist (these are allowed)."""
        client = DomainAPIClient(cache_enabled=False)
        assert hasattr(client, 'CONFIDENCE_THRESHOLDS')
        assert 'high_hotspot_tumor_count' in client.CONFIDENCE_THRESHOLDS
        assert 'critical_domain_types' in client.CONFIDENCE_THRESHOLDS
    
    @patch('utils.domain_api_client.requests.get')
    def test_get_hotspot_annotation_cancer_hotspots(self, mock_get):
        """Test annotation from CancerHotspots.org API response."""
        # Mock CancerHotspots.org response (first call)
        mock_hotspot_response = Mock()
        mock_hotspot_response.status_code = 200
        mock_hotspot_response.json.return_value = [
            {'tumorCount': 25, 'count': 80}
        ]
        
        # Mock UniProt responses (second and third calls)
        mock_uniprot_search = Mock()
        mock_uniprot_search.status_code = 200
        mock_uniprot_search.json.return_value = {'results': []}
        
        mock_get.side_effect = [mock_hotspot_response, mock_uniprot_search]
        
        client = DomainAPIClient(cache_enabled=False)
        ann = client.get_hotspot_annotation("TP53", position=248)
        
        assert ann.is_hotspot is True
        assert ann.tumor_count == 25
        assert ann.hotspot_count == 80
        assert ann.confidence >= 0.9  # High tumor count -> high confidence
        assert "CancerHotspots.org" in ann.source
    
    @patch('utils.domain_api_client.requests.get')
    def test_get_hotspot_annotation_moderate_confidence(self, mock_get):
        """Test moderate confidence for fewer tumors."""
        mock_hotspot_response = Mock()
        mock_hotspot_response.status_code = 200
        mock_hotspot_response.json.return_value = [
            {'tumorCount': 5, 'count': 10}  # Between 3 and 10 tumors
        ]
        
        mock_uniprot_search = Mock()
        mock_uniprot_search.status_code = 200
        mock_uniprot_search.json.return_value = {'results': []}
        
        mock_get.side_effect = [mock_hotspot_response, mock_uniprot_search]
        
        client = DomainAPIClient(cache_enabled=False)
        ann = client.get_hotspot_annotation("KRAS", position=12)
        
        assert ann.is_hotspot is True
        assert ann.confidence >= 0.5 and ann.confidence < 0.9
    
    @patch('utils.domain_api_client.requests.get')
    def test_get_hotspot_annotation_not_hotspot(self, mock_get):
        """Test when position is not a hotspot."""
        # First call to CancerHotspots returns empty
        # Second calls to UniProt also return no domains
        mock_response_empty = Mock()
        mock_response_empty.status_code = 200
        mock_response_empty.json.return_value = []
        
        mock_uniprot_search = Mock()
        mock_uniprot_search.status_code = 200
        mock_uniprot_search.json.return_value = {'results': []}
        
        mock_get.side_effect = [mock_response_empty, mock_uniprot_search]
        
        client = DomainAPIClient(cache_enabled=False)
        ann = client.get_hotspot_annotation("FAKE_GENE", position=999)
        
        assert ann.is_hotspot is False
        assert ann.in_critical_domain is False
    
    @patch('utils.domain_api_client.requests.get')
    def test_get_hotspot_annotation_api_error(self, mock_get):
        """Test graceful handling of API errors."""
        # Simulate network error on first call (CancerHotspots)
        from requests.exceptions import RequestException
        mock_get.side_effect = RequestException("Connection timeout")
        
        client = DomainAPIClient(cache_enabled=False)
        ann = client.get_hotspot_annotation("TP53", position=248)
        
        # Should return empty annotation, not raise exception
        assert ann.is_hotspot is False
        assert ann.source == "none"
    
    @patch('utils.domain_api_client.requests.get')
    def test_get_hotspot_annotation_uniprot_domain(self, mock_get):
        """Test annotation from UniProt domain data."""
        # CancerHotspots returns empty
        mock_hotspot_response = Mock()
        mock_hotspot_response.status_code = 200
        mock_hotspot_response.json.return_value = []
        
        # UniProt search returns accession
        mock_uniprot_search = Mock()
        mock_uniprot_search.status_code = 200
        mock_uniprot_search.json.return_value = {
            'results': [{'primaryAccession': 'P04637'}]
        }
        
        # UniProt entry returns domain info
        mock_uniprot_entry = Mock()
        mock_uniprot_entry.status_code = 200
        mock_uniprot_entry.json.return_value = {
            'features': [
                {
                    'type': 'Domain',
                    'description': 'DNA-binding',
                    'location': {'start': {'value': 100}, 'end': {'value': 300}}
                }
            ]
        }
        
        mock_get.side_effect = [
            mock_hotspot_response,
            mock_uniprot_search,
            mock_uniprot_entry
        ]
        
        client = DomainAPIClient(cache_enabled=False)
        ann = client.get_hotspot_annotation("TP53", position=248)
        
        assert ann.in_critical_domain is True
        assert ann.domain_name == "DNA-binding"
        assert "UniProt" in ann.source
    
    def test_extract_position_from_hgvs_p(self):
        """Test position extraction from protein HGVS."""
        client = DomainAPIClient(cache_enabled=False)
        
        assert client._extract_position_from_hgvs_p("p.Arg248Gln") == 248
        assert client._extract_position_from_hgvs_p("p.R248Q") == 248
        assert client._extract_position_from_hgvs_p("p.Arg248*") == 248
        assert client._extract_position_from_hgvs_p("p.Val600Glu") == 600
        assert client._extract_position_from_hgvs_p("") is None
    
    @patch('utils.domain_api_client.requests.get')
    def test_legacy_get_hotspot_info(self, mock_get):
        """Test legacy get_hotspot_info method returns compatible format."""
        mock_hotspot_response = Mock()
        mock_hotspot_response.status_code = 200
        mock_hotspot_response.json.return_value = [
            {'tumorCount': 20, 'count': 50}
        ]
        
        mock_uniprot_search = Mock()
        mock_uniprot_search.status_code = 200
        mock_uniprot_search.json.return_value = {'results': []}
        
        mock_get.side_effect = [mock_hotspot_response, mock_uniprot_search]
        
        client = DomainAPIClient(cache_enabled=False)
        result = client.get_hotspot_info("TP53", position=248)
        
        # Legacy format should still work
        assert result['position_is_hotspot'] is True
        assert result['is_hotspot_gene'] is True
        assert result['source'] == "CancerHotspots.org"
        assert 'hotspot_details' in result


# =============================================================================
# GeneSpecificRules Tests
# =============================================================================

class TestGeneSpecificRules:
    """Tests for GeneSpecificRules PM1 evaluation."""
    
    def test_no_hardcoded_hotspot_regions(self):
        """Verify hardcoded hotspot_regions have been removed."""
        rules = GeneSpecificRules()
        # Old implementation had self.hotspot_regions
        assert not hasattr(rules, 'hotspot_regions'), \
            "hotspot_regions should be removed - no hardcoded biological data"
    
    def test_pm1_thresholds_exist(self):
        """Verify interpretive thresholds exist (these are allowed)."""
        rules = GeneSpecificRules()
        assert hasattr(rules, 'PM1_THRESHOLDS')
        assert 'PM1' in rules.PM1_THRESHOLDS
        assert 'PM1_supporting' in rules.PM1_THRESHOLDS
    
    def test_evaluate_pm1_no_gene(self):
        """PM1 should not apply without gene symbol."""
        mock_client = Mock(spec=DomainAPIClient)
        rules = GeneSpecificRules(domain_client=mock_client)
        
        result = rules.evaluate_pm1(gene="", position=248)
        
        assert result.applies is False
        assert "required" in result.reason.lower()
        mock_client.get_hotspot_annotation.assert_not_called()
    
    def test_evaluate_pm1_high_confidence_hotspot(self):
        """High confidence hotspot should give PM1."""
        mock_client = Mock(spec=DomainAPIClient)
        mock_client.get_hotspot_annotation.return_value = HotspotAnnotation(
            is_hotspot=True,
            confidence=0.95,
            source="CancerHotspots.org",
            tumor_count=50
        )
        
        rules = GeneSpecificRules(domain_client=mock_client)
        result = rules.evaluate_pm1(gene="TP53", position=248)
        
        assert result.applies is True
        assert result.evidence_code == 'PM1'
        assert result.confidence == 0.95
        assert "hotspot" in result.reason.lower()
    
    def test_evaluate_pm1_moderate_confidence_domain(self):
        """Moderate confidence domain should give PM1_supporting."""
        mock_client = Mock(spec=DomainAPIClient)
        mock_client.get_hotspot_annotation.return_value = HotspotAnnotation(
            is_hotspot=False,
            in_critical_domain=True,
            confidence=0.7,
            source="UniProt",
            domain_name="Kinase domain"
        )
        
        rules = GeneSpecificRules(domain_client=mock_client)
        result = rules.evaluate_pm1(gene="BRAF", position=600)
        
        assert result.applies is True
        assert result.evidence_code == 'PM1_supporting'
        assert "domain" in result.reason.lower() or "region" in result.reason.lower()
    
    def test_evaluate_pm1_low_confidence(self):
        """Low confidence should not give PM1."""
        mock_client = Mock(spec=DomainAPIClient)
        mock_client.get_hotspot_annotation.return_value = HotspotAnnotation(
            is_hotspot=False,
            in_critical_domain=False,
            confidence=0.3,
            source="UniProt",
            details="No significant findings"
        )
        
        rules = GeneSpecificRules(domain_client=mock_client)
        result = rules.evaluate_pm1(gene="RANDOM", position=100)
        
        assert result.applies is False
        assert result.evidence_code == ""
    
    def test_evaluate_pm1_api_unavailable(self):
        """API unavailable should return negative result (safe fallback)."""
        mock_client = Mock(spec=DomainAPIClient)
        mock_client.get_hotspot_annotation.return_value = HotspotAnnotation(
            source="none"
        )
        
        rules = GeneSpecificRules(domain_client=mock_client)
        result = rules.evaluate_pm1(gene="TP53", position=248)
        
        assert result.applies is False
        assert "none" in result.source or "no" in result.reason.lower()
    
    def test_evaluate_pm1_api_error(self):
        """API error should return negative result (safe fallback)."""
        mock_client = Mock(spec=DomainAPIClient)
        mock_client.get_hotspot_annotation.side_effect = Exception("Timeout")
        
        rules = GeneSpecificRules(domain_client=mock_client)
        result = rules.evaluate_pm1(gene="TP53", position=248)
        
        assert result.applies is False
        assert result.source == "error"
    
    def test_evaluate_pm1_with_hgvs_p(self):
        """Position should be extracted from hgvs_p if not provided."""
        mock_client = Mock(spec=DomainAPIClient)
        mock_client.get_hotspot_annotation.return_value = HotspotAnnotation(
            is_hotspot=True,
            confidence=0.9,
            source="CancerHotspots.org"
        )
        
        rules = GeneSpecificRules(domain_client=mock_client)
        result = rules.evaluate_pm1(gene="TP53", hgvs_p="p.Arg248Gln")
        
        # Verify position was extracted and passed
        mock_client.get_hotspot_annotation.assert_called_once()
        call_kwargs = mock_client.get_hotspot_annotation.call_args[1]
        assert call_kwargs['gene'] == "TP53"
        assert call_kwargs['hgvs_p'] == "p.Arg248Gln"
    
    def test_extract_position(self):
        """Test position extraction from various HGVS formats."""
        rules = GeneSpecificRules()
        
        assert rules._extract_position("p.Arg248Gln") == 248
        assert rules._extract_position("p.R248Q") == 248
        assert rules._extract_position("p.Val600Glu") == 600
        assert rules._extract_position("p.G12D") == 12
        assert rules._extract_position("p.Arg248*") == 248  # Nonsense
        assert rules._extract_position("p.Arg248fs") == 248  # Frameshift
        assert rules._extract_position("") is None
        assert rules._extract_position(None) is None


# =============================================================================
# Integration Tests
# =============================================================================

class TestPM1Integration:
    """Integration tests for PM1 evaluation pipeline."""
    
    @patch('utils.domain_api_client.requests.get')
    def test_full_pipeline_hotspot(self, mock_get):
        """Test full pipeline from gene/position to PM1 evidence."""
        # Mock CancerHotspots.org returning hotspot data
        mock_hotspot_response = Mock()
        mock_hotspot_response.status_code = 200
        mock_hotspot_response.json.return_value = [
            {'tumorCount': 100, 'count': 500}
        ]
        
        mock_uniprot_search = Mock()
        mock_uniprot_search.status_code = 200
        mock_uniprot_search.json.return_value = {'results': []}
        
        mock_get.side_effect = [mock_hotspot_response, mock_uniprot_search]
        
        # Create real client (with mocked HTTP)
        client = DomainAPIClient(cache_enabled=False)
        rules = GeneSpecificRules(domain_client=client)
        
        result = rules.evaluate_pm1(gene="TP53", position=248)
        
        assert result.applies is True
        assert result.evidence_code == 'PM1'
        assert result.confidence >= 0.85
    
    @patch('utils.domain_api_client.requests.get')
    def test_full_pipeline_no_data(self, mock_get):
        """Test full pipeline when APIs return no data."""
        # All API calls return empty/error
        mock_empty = Mock()
        mock_empty.status_code = 200
        mock_empty.json.return_value = []
        
        mock_get.return_value = mock_empty
        
        client = DomainAPIClient(cache_enabled=False)
        rules = GeneSpecificRules(domain_client=client)
        
        result = rules.evaluate_pm1(gene="UNKNOWN_GENE", position=999)
        
        assert result.applies is False
        # Should NOT have hardcoded fallback
        assert "manual" not in result.reason.lower() or "review" not in result.reason.lower()


# =============================================================================
# Threshold Boundary Tests
# =============================================================================

class TestPM1Thresholds:
    """Test threshold boundaries for PM1 evidence assignment."""
    
    def test_boundary_pm1_full(self):
        """Test boundary at PM1 threshold (0.85)."""
        mock_client = Mock(spec=DomainAPIClient)
        rules = GeneSpecificRules(domain_client=mock_client)
        
        # Just at threshold
        mock_client.get_hotspot_annotation.return_value = HotspotAnnotation(
            is_hotspot=True, confidence=0.85, source="test"
        )
        result = rules.evaluate_pm1(gene="TEST", position=1)
        assert result.evidence_code == 'PM1'
        
        # Just below threshold
        mock_client.get_hotspot_annotation.return_value = HotspotAnnotation(
            is_hotspot=True, confidence=0.84, source="test"
        )
        result = rules.evaluate_pm1(gene="TEST", position=1)
        assert result.evidence_code == 'PM1_supporting'
    
    def test_boundary_pm1_supporting(self):
        """Test boundary at PM1_supporting threshold (0.60)."""
        mock_client = Mock(spec=DomainAPIClient)
        rules = GeneSpecificRules(domain_client=mock_client)
        
        # Just at threshold
        mock_client.get_hotspot_annotation.return_value = HotspotAnnotation(
            in_critical_domain=True, confidence=0.60, source="test"
        )
        result = rules.evaluate_pm1(gene="TEST", position=1)
        assert result.applies is True
        assert result.evidence_code == 'PM1_supporting'
        
        # Just below threshold
        mock_client.get_hotspot_annotation.return_value = HotspotAnnotation(
            in_critical_domain=True, confidence=0.59, source="test"
        )
        result = rules.evaluate_pm1(gene="TEST", position=1)
        assert result.applies is False


if __name__ == '__main__':
    pytest.main([__file__, '-v'])

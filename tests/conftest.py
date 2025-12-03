"""
Pytest Configuration and Fixtures
=================================

Provides fixtures for offline testing of the ACMG classifier.
All tests run without network calls by mocking API clients.
"""

import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

# Add src to path for imports
src_path = Path(__file__).parent.parent / "src"
sys.path.insert(0, str(src_path))


@pytest.fixture(autouse=True)
def mock_api_client():
    """
    Mock the APIClient to prevent any network calls during tests.
    
    This fixture is auto-used for all tests to ensure offline execution.
    The APIClient is imported inside EvidenceEvaluator.__init__, so we
    patch it at the utils.api_client module level.
    """
    mock_client = MagicMock()
    
    # Mock gene constraint data (gnomAD-like response)
    mock_client.get_gene_constraint.return_value = {
        'is_lof_intolerant': True,
        'classification': 'LOF_intolerant',
        'source': 'gnomAD',
        'confidence': 'high',
        'pLI': 0.99,
        'LOEUF': 0.2,
        'oe_lof': 0.15
    }
    
    # Mock ClinGen gene validity
    mock_client.get_clingen_gene_validity.return_value = {
        'supports_lof_pathogenicity': True,
        'source': 'ClinGen',
        'confidence': 'high',
        'lof_diseases': ['Hereditary breast and ovarian cancer'],
        'primary_mechanism': 'loss_of_function'
    }
    
    # Mock ClinGen dosage sensitivity
    mock_client.get_clingen_dosage_sensitivity.return_value = {
        'haploinsufficiency_score': 3,
        'triplosensitivity_score': None,
        'pvs1_recommendation': 'very_strong',
        'source': 'ClinGen Dosage',
        'confidence': 'high'
    }
    
    # Patch the APIClient class at the module where it's imported from
    with patch('utils.api_client.APIClient', return_value=mock_client):
        yield mock_client


@pytest.fixture
def mock_api_client_lof_tolerant(mock_api_client):
    """
    Configure mock API client for LOF-tolerant gene scenarios.
    """
    mock_api_client.get_gene_constraint.return_value = {
        'is_lof_intolerant': False,
        'classification': 'LOF_tolerant',
        'source': 'gnomAD',
        'confidence': 'high',
        'pLI': 0.01,
        'LOEUF': 1.5,
        'oe_lof': 1.2
    }
    
    mock_api_client.get_clingen_gene_validity.return_value = {
        'supports_lof_pathogenicity': False,
        'source': 'ClinGen',
        'confidence': 'low',
        'lof_diseases': [],
        'primary_mechanism': None
    }
    
    mock_api_client.get_clingen_dosage_sensitivity.return_value = {
        'haploinsufficiency_score': 0,
        'triplosensitivity_score': None,
        'pvs1_recommendation': 'not_applicable',
        'source': 'ClinGen Dosage',
        'confidence': 'high'
    }
    
    return mock_api_client

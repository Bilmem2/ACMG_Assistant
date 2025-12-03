"""
Tests for Interactive Evidence Collector Module
================================================

This module tests the InteractiveEvidenceCollector and its pure functions
for mapping user input to ACMG evidence codes.

Test Categories:
1. Pure function tests (no I/O, deterministic)
2. ManualEvidence dataclass tests
3. MockInputProvider tests
4. InteractiveEvidenceCollector with mocked input
5. Integration tests with EvidenceEvaluator

Author: Can Sevilmiş
License: MIT
"""

import pytest
import sys
import os

# Add src to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from core.interactive_evidence import (
    ManualEvidence,
    InteractiveEvidenceCollector,
    InputProvider,
    MockInputProvider,
    EvidenceStrength,
    map_functional_studies_to_evidence,
    map_segregation_to_evidence,
    map_case_control_to_evidence,
)


# =============================================================================
# Pure Function Tests - Functional Studies (PS3/BS3)
# =============================================================================
class TestFunctionalStudiesMapping:
    """Tests for map_functional_studies_to_evidence pure function."""
    
    def test_no_studies_returns_none(self):
        """No functional studies → no evidence."""
        code, explanation = map_functional_studies_to_evidence(0, 0, 'high')
        assert code is None
        assert 'No functional studies' in explanation
    
    def test_single_high_quality_damaging_returns_ps3_moderate(self):
        """1 high-quality damaging study → PS3_moderate."""
        code, explanation = map_functional_studies_to_evidence(1, 0, 'high')
        assert code == 'PS3_moderate'
        assert 'high-quality' in explanation.lower()
    
    def test_multiple_high_quality_damaging_returns_ps3_strong(self):
        """≥2 high-quality damaging studies → PS3_strong."""
        code, explanation = map_functional_studies_to_evidence(2, 0, 'high')
        assert code == 'PS3_strong'
        
        code, explanation = map_functional_studies_to_evidence(5, 0, 'high')
        assert code == 'PS3_strong'
    
    def test_single_moderate_quality_damaging_returns_ps3_supporting(self):
        """1 moderate-quality damaging study → PS3_supporting."""
        code, explanation = map_functional_studies_to_evidence(1, 0, 'moderate')
        assert code == 'PS3_supporting'
    
    def test_multiple_moderate_quality_damaging_returns_ps3_moderate(self):
        """≥2 moderate-quality damaging studies → PS3_moderate."""
        code, explanation = map_functional_studies_to_evidence(2, 0, 'moderate')
        assert code == 'PS3_moderate'
    
    def test_single_low_quality_damaging_returns_none(self):
        """1 low-quality damaging study → insufficient."""
        code, explanation = map_functional_studies_to_evidence(1, 0, 'low')
        assert code is None
        assert 'insufficient' in explanation.lower()
    
    def test_multiple_low_quality_damaging_returns_ps3_supporting(self):
        """≥2 low-quality damaging studies → PS3_supporting."""
        code, explanation = map_functional_studies_to_evidence(2, 0, 'low')
        assert code == 'PS3_supporting'
    
    def test_single_high_quality_benign_returns_bs3_moderate(self):
        """1 high-quality benign study → BS3_moderate."""
        code, explanation = map_functional_studies_to_evidence(0, 1, 'high')
        assert code == 'BS3_moderate'
    
    def test_multiple_high_quality_benign_returns_bs3_strong(self):
        """≥2 high-quality benign studies → BS3_strong."""
        code, explanation = map_functional_studies_to_evidence(0, 2, 'high')
        assert code == 'BS3_strong'
    
    def test_conflicting_predominantly_damaging_returns_ps3_supporting(self):
        """Conflicting but >2x damaging → PS3_supporting."""
        code, explanation = map_functional_studies_to_evidence(5, 2, 'high')
        assert code == 'PS3_supporting'
        assert 'conflicting' in explanation.lower()
    
    def test_conflicting_predominantly_benign_returns_bs3_supporting(self):
        """Conflicting but >2x benign → BS3_supporting."""
        code, explanation = map_functional_studies_to_evidence(2, 5, 'high')
        assert code == 'BS3_supporting'
        assert 'conflicting' in explanation.lower()
    
    def test_conflicting_no_clear_direction_returns_none(self):
        """Conflicting with similar counts → no evidence."""
        code, explanation = map_functional_studies_to_evidence(3, 3, 'high')
        assert code is None
        assert 'conflicting' in explanation.lower()


# =============================================================================
# Pure Function Tests - Segregation (PP1/BS4)
# =============================================================================
class TestSegregationMapping:
    """Tests for map_segregation_to_evidence pure function."""
    
    def test_cosegregates_5_or_more_returns_pp1_strong(self):
        """≥5 affected carriers → PP1_strong."""
        code, explanation = map_segregation_to_evidence('cosegregates', affected_carriers=5)
        assert code == 'PP1_strong'
        
        code, explanation = map_segregation_to_evidence('cosegregates', affected_carriers=10)
        assert code == 'PP1_strong'
    
    def test_cosegregates_3_to_4_returns_pp1_moderate(self):
        """3-4 affected carriers → PP1_moderate."""
        code, explanation = map_segregation_to_evidence('cosegregates', affected_carriers=3)
        assert code == 'PP1_moderate'
        
        code, explanation = map_segregation_to_evidence('cosegregates', affected_carriers=4)
        assert code == 'PP1_moderate'
    
    def test_cosegregates_2_returns_pp1_supporting(self):
        """2 affected carriers → PP1_supporting."""
        code, explanation = map_segregation_to_evidence('cosegregates', affected_carriers=2)
        assert code == 'PP1_supporting'
    
    def test_cosegregates_1_returns_none(self):
        """1 affected carrier → insufficient."""
        code, explanation = map_segregation_to_evidence('cosegregates', affected_carriers=1)
        assert code is None
        assert 'insufficient' in explanation.lower()
    
    def test_does_not_segregate_3_or_more_returns_bs4_moderate(self):
        """≥3 affected non-carriers → BS4_moderate."""
        code, explanation = map_segregation_to_evidence('does_not_segregate', affected_noncarriers=3)
        assert code == 'BS4_moderate'
    
    def test_does_not_segregate_2_returns_bs4_supporting(self):
        """2 affected non-carriers → BS4_supporting."""
        code, explanation = map_segregation_to_evidence('does_not_segregate', affected_noncarriers=2)
        assert code == 'BS4_supporting'
    
    def test_does_not_segregate_1_returns_none(self):
        """1 affected non-carrier → insufficient."""
        code, explanation = map_segregation_to_evidence('does_not_segregate', affected_noncarriers=1)
        assert code is None
    
    def test_mixed_pattern_returns_none(self):
        """Mixed segregation pattern → no evidence."""
        code, explanation = map_segregation_to_evidence('mixed')
        assert code is None
        assert 'mixed' in explanation.lower()


# =============================================================================
# Pure Function Tests - Case-Control (PS4)
# =============================================================================
class TestCaseControlMapping:
    """Tests for map_case_control_to_evidence pure function."""
    
    def test_strong_enrichment_returns_ps4_strong(self):
        """≥5 cases with OR≥5 → PS4_strong."""
        # 10 cases with variant out of 100, 1 control out of 1000
        # OR = (10 * 999) / (90 * 1) ≈ 111
        code, explanation = map_case_control_to_evidence(10, 100, 1, 1000)
        assert code == 'PS4_strong'
    
    def test_moderate_enrichment_returns_ps4_moderate(self):
        """≥3 cases with OR≥3 → PS4_moderate."""
        # 3 cases out of 100, 0 controls out of 500
        code, explanation = map_case_control_to_evidence(3, 100, 0, 500)
        assert code == 'PS4_moderate'
    
    def test_supporting_enrichment_returns_ps4_supporting(self):
        """≥2 cases with OR≥2 → PS4_supporting."""
        # 2 cases out of 100, 0 controls out of 200
        code, explanation = map_case_control_to_evidence(2, 100, 0, 200)
        assert code == 'PS4_supporting'
    
    def test_no_cases_returns_none(self):
        """No cases → no evidence."""
        code, explanation = map_case_control_to_evidence(0, 100, 1, 1000)
        assert code is None
        assert 'insufficient' in explanation.lower()
    
    def test_insufficient_enrichment_returns_none(self):
        """Low OR → no evidence."""
        # Similar frequency in cases and controls
        code, explanation = map_case_control_to_evidence(10, 100, 10, 100)
        assert code is None
        assert 'not significant' in explanation.lower() or 'insufficient' in explanation.lower()


# =============================================================================
# ManualEvidence Dataclass Tests
# =============================================================================
class TestManualEvidence:
    """Tests for ManualEvidence dataclass."""
    
    def test_empty_initialization(self):
        """Empty ManualEvidence should have no codes."""
        ev = ManualEvidence()
        assert ev.codes == []
        assert ev.explanations == {}
        assert ev.has_evidence() is False
    
    def test_add_evidence(self):
        """Adding evidence should update codes and explanations."""
        ev = ManualEvidence()
        ev.add_evidence('PS3_strong', 'Functional studies show damage')
        
        assert 'PS3_strong' in ev.codes
        assert ev.explanations['PS3_strong'] == 'Functional studies show damage'
        assert ev.has_evidence() is True
    
    def test_add_duplicate_evidence(self):
        """Duplicate codes should not create duplicates in list."""
        ev = ManualEvidence()
        ev.add_evidence('PS3_strong', 'First explanation')
        ev.add_evidence('PS3_strong', 'Updated explanation')
        
        assert ev.codes.count('PS3_strong') == 1
        assert ev.explanations['PS3_strong'] == 'Updated explanation'
    
    def test_get_pathogenic_codes(self):
        """Should return only pathogenic codes."""
        ev = ManualEvidence()
        ev.add_evidence('PS3_strong', 'Damaging')
        ev.add_evidence('PP1_moderate', 'Segregation')
        ev.add_evidence('BS3_supporting', 'Benign')
        
        pathogenic = ev.get_pathogenic_codes()
        assert 'PS3_strong' in pathogenic
        assert 'PP1_moderate' in pathogenic
        assert 'BS3_supporting' not in pathogenic
    
    def test_get_benign_codes(self):
        """Should return only benign codes."""
        ev = ManualEvidence()
        ev.add_evidence('PS3_strong', 'Damaging')
        ev.add_evidence('BS3_supporting', 'Benign')
        ev.add_evidence('BP6_supporting', 'Benign source')
        
        benign = ev.get_benign_codes()
        assert 'BS3_supporting' in benign
        assert 'BP6_supporting' in benign
        assert 'PS3_strong' not in benign
    
    def test_merge_with(self):
        """Merging should combine codes and explanations."""
        ev1 = ManualEvidence()
        ev1.add_evidence('PS3_strong', 'Study 1')
        
        ev2 = ManualEvidence()
        ev2.add_evidence('PP1_moderate', 'Segregation')
        
        merged = ev1.merge_with(ev2)
        
        assert 'PS3_strong' in merged.codes
        assert 'PP1_moderate' in merged.codes
        assert merged.explanations['PS3_strong'] == 'Study 1'
        assert merged.explanations['PP1_moderate'] == 'Segregation'


# =============================================================================
# MockInputProvider Tests
# =============================================================================
class TestMockInputProvider:
    """Tests for MockInputProvider."""
    
    def test_scripted_responses(self):
        """Should return responses in order."""
        mock = MockInputProvider(['yes', 'no', '3'])
        
        assert mock.prompt('Q1: ') == 'yes'
        assert mock.prompt('Q2: ') == 'no'
        assert mock.prompt('Q3: ') == '3'
    
    def test_exhausted_responses(self):
        """Should return empty string when responses exhausted."""
        mock = MockInputProvider(['yes'])
        
        assert mock.prompt('Q1: ') == 'yes'
        assert mock.prompt('Q2: ') == ''
    
    def test_prompt_yes_no_yes(self):
        """Should parse 'y' as True."""
        mock = MockInputProvider(['y'])
        result = mock.prompt_yes_no('Continue?')
        assert result is True
    
    def test_prompt_yes_no_no(self):
        """Should parse 'n' as False."""
        mock = MockInputProvider(['n'])
        result = mock.prompt_yes_no('Continue?')
        assert result is False
    
    def test_prompt_yes_no_default(self):
        """Should use default on empty response."""
        mock = MockInputProvider([''])
        result = mock.prompt_yes_no('Continue?', default=True)
        assert result is True


# =============================================================================
# InteractiveEvidenceCollector Tests with Mock Input
# =============================================================================
class TestInteractiveEvidenceCollectorMocked:
    """Tests for InteractiveEvidenceCollector using mocked input."""
    
    def test_collect_ps3_bs3_damaging(self):
        """Collect PS3 with damaging functional studies."""
        # Responses: has_studies=yes, damaging=2, benign=0, quality=high
        mock = MockInputProvider(['y', '2', '0', 'high'])
        collector = InteractiveEvidenceCollector(input_provider=mock, show_prompts=False)
        
        evidence = collector.collect_ps3_bs3()
        
        assert evidence.has_evidence()
        assert any('PS3' in code for code in evidence.codes)
    
    def test_collect_ps3_bs3_benign(self):
        """Collect BS3 with benign functional studies."""
        # Responses: has_studies=yes, damaging=0, benign=2, quality=high
        mock = MockInputProvider(['y', '0', '2', 'high'])
        collector = InteractiveEvidenceCollector(input_provider=mock, show_prompts=False)
        
        evidence = collector.collect_ps3_bs3()
        
        assert evidence.has_evidence()
        assert any('BS3' in code for code in evidence.codes)
    
    def test_collect_ps3_bs3_no_studies(self):
        """No functional studies → no evidence."""
        mock = MockInputProvider(['n'])
        collector = InteractiveEvidenceCollector(input_provider=mock, show_prompts=False)
        
        evidence = collector.collect_ps3_bs3()
        
        assert not evidence.has_evidence()
    
    def test_collect_ps4_significant_enrichment(self):
        """Collect PS4 with significant case-control enrichment."""
        # has_data=yes, cases_with=5, total_cases=50, controls_with=0, total_controls=1000
        mock = MockInputProvider(['y', '5', '50', '0', '1000'])
        collector = InteractiveEvidenceCollector(input_provider=mock, show_prompts=False)
        
        evidence = collector.collect_ps4()
        
        assert evidence.has_evidence()
        assert any('PS4' in code for code in evidence.codes)
    
    def test_collect_ps4_no_data(self):
        """No case-control data → no evidence."""
        mock = MockInputProvider(['n'])
        collector = InteractiveEvidenceCollector(input_provider=mock, show_prompts=False)
        
        evidence = collector.collect_ps4()
        
        assert not evidence.has_evidence()
    
    def test_collect_pp1_cosegregation(self):
        """Collect PP1 with cosegregation."""
        # has_family=yes, pattern=cosegregates, affected_carriers=4, unaffected_noncarriers=2
        mock = MockInputProvider(['y', 'cosegregates', '4', '2'])
        collector = InteractiveEvidenceCollector(input_provider=mock, show_prompts=False)
        
        evidence = collector.collect_pp1_bs4()
        
        assert evidence.has_evidence()
        assert any('PP1' in code for code in evidence.codes)
    
    def test_collect_bs4_non_segregation(self):
        """Collect BS4 with non-segregation."""
        # has_family=yes, pattern=does_not_segregate, affected_noncarriers=3
        mock = MockInputProvider(['y', 'does_not_segregate', '3'])
        collector = InteractiveEvidenceCollector(input_provider=mock, show_prompts=False)
        
        evidence = collector.collect_pp1_bs4()
        
        assert evidence.has_evidence()
        assert any('BS4' in code for code in evidence.codes)
    
    def test_collect_ps1_same_aa_pathogenic(self):
        """Collect PS1 with same AA change pathogenic."""
        # has_same_aa=yes, splice_concern=no
        mock = MockInputProvider(['y', 'n'])
        collector = InteractiveEvidenceCollector(input_provider=mock, show_prompts=False)
        
        evidence = collector.collect_ps1_pm5()
        
        assert evidence.has_evidence()
        assert 'PS1_strong' in evidence.codes
    
    def test_collect_pm5_different_aa_same_position(self):
        """Collect PM5 with different AA at same position."""
        # has_same_aa=no, has_same_position=yes, classification=pathogenic
        mock = MockInputProvider(['n', 'y', 'pathogenic'])
        collector = InteractiveEvidenceCollector(input_provider=mock, show_prompts=False)
        
        evidence = collector.collect_ps1_pm5()
        
        assert evidence.has_evidence()
        assert 'PM5_moderate' in evidence.codes
    
    def test_collect_pp5_pathogenic_assertion(self):
        """Collect PP5 with reputable pathogenic assertion."""
        # use_criteria=yes, has_path_assertion=yes, has_conflicts=no, star_level=3_star
        mock = MockInputProvider(['y', 'y', 'n', '3_star'])
        collector = InteractiveEvidenceCollector(input_provider=mock, show_prompts=False)
        
        evidence = collector.collect_pp5_bp6()
        
        assert evidence.has_evidence()
        assert 'PP5_supporting' in evidence.codes
    
    def test_collect_bp6_benign_assertion(self):
        """Collect BP6 with reputable benign assertion."""
        # use_criteria=yes, has_path_assertion=no, has_benign_assertion=yes, has_conflicts=no
        mock = MockInputProvider(['y', 'n', 'y', 'n'])
        collector = InteractiveEvidenceCollector(input_provider=mock, show_prompts=False)
        
        evidence = collector.collect_pp5_bp6()
        
        assert evidence.has_evidence()
        assert 'BP6_supporting' in evidence.codes
    
    def test_collect_pp5_bp6_deprecated_not_used(self):
        """User chooses not to use deprecated PP5/BP6."""
        mock = MockInputProvider(['n'])
        collector = InteractiveEvidenceCollector(input_provider=mock, show_prompts=False)
        
        evidence = collector.collect_pp5_bp6()
        
        assert not evidence.has_evidence()


# =============================================================================
# Integration Tests
# =============================================================================
class TestEvidenceEvaluatorIntegration:
    """Tests for EvidenceEvaluator integration with manual evidence."""
    
    def test_set_and_get_manual_evidence(self):
        """Should store and retrieve manual evidence."""
        from core.evidence_evaluator import EvidenceEvaluator
        
        evaluator = EvidenceEvaluator(test_mode=True)
        manual_ev = ManualEvidence()
        manual_ev.add_evidence('PS3_strong', 'Test evidence')
        
        evaluator.set_manual_evidence(manual_ev)
        retrieved = evaluator.get_manual_evidence()
        
        assert retrieved is not None
        assert 'PS3_strong' in retrieved.codes
    
    def test_merge_manual_with_automated(self):
        """Should merge manual evidence into automated results."""
        from core.evidence_evaluator import EvidenceEvaluator
        
        evaluator = EvidenceEvaluator(test_mode=True)
        
        # Create manual evidence
        manual_ev = ManualEvidence()
        manual_ev.add_evidence('PS3_moderate', 'Functional studies show loss of function')
        manual_ev.add_evidence('PP1_supporting', 'Cosegregation in 2 affected family members')
        evaluator.set_manual_evidence(manual_ev)
        
        # Create mock automated results
        automated_results = {
            'pathogenic_criteria': {
                'PP3': {'applies': True, 'strength': 'Supporting', 'details': 'Computational evidence'}
            },
            'benign_criteria': {},
            'applied_criteria': {
                'PP3': {'applies': True, 'strength': 'Supporting', 'details': 'Computational evidence'}
            },
            'evidence_details': {},
            'vampp_score': 0.75,
            'statistical_tests': {}
        }
        
        # Merge
        merged = evaluator.merge_manual_with_automated(automated_results)
        
        # Check that manual evidence was added
        assert 'PS3' in merged['applied_criteria']
        assert merged['applied_criteria']['PS3']['applies'] is True
        assert merged['applied_criteria']['PS3']['data_source'] == 'manual_evidence'
        
        assert 'PP1' in merged['applied_criteria']
        assert merged['applied_criteria']['PP1']['strength'] == 'Supporting'
        
        # Check that automated evidence is preserved
        assert 'PP3' in merged['applied_criteria']
        
        # Check that manual_evidence summary is included
        assert 'manual_evidence' in merged
        assert len(merged['manual_evidence']['codes']) == 2
    
    def test_merge_with_no_manual_evidence(self):
        """Should return original results if no manual evidence."""
        from core.evidence_evaluator import EvidenceEvaluator
        
        evaluator = EvidenceEvaluator(test_mode=True)
        
        automated_results = {
            'pathogenic_criteria': {'PP3': {'applies': True}},
            'benign_criteria': {},
            'applied_criteria': {'PP3': {'applies': True}},
            'evidence_details': {},
            'vampp_score': 0.5,
            'statistical_tests': {}
        }
        
        merged = evaluator.merge_manual_with_automated(automated_results)
        
        # Should be essentially unchanged
        assert 'PP3' in merged['applied_criteria']
        assert 'manual_evidence' not in merged


# =============================================================================
# Edge Case Tests
# =============================================================================
class TestEdgeCases:
    """Tests for edge cases and boundary conditions."""
    
    def test_zero_affected_carriers(self):
        """Zero affected carriers should return None."""
        code, _ = map_segregation_to_evidence('cosegregates', affected_carriers=0)
        assert code is None
    
    def test_negative_study_count(self):
        """Negative counts should be handled gracefully."""
        # This shouldn't happen in practice, but function should handle it
        code, _ = map_functional_studies_to_evidence(-1, 0, 'high')
        assert code is None
    
    def test_empty_quality_string(self):
        """Empty quality should default to no evidence."""
        code, _ = map_functional_studies_to_evidence(1, 0, '')
        assert code is None
    
    def test_case_control_zero_total_cases(self):
        """Zero total cases should return None."""
        code, _ = map_case_control_to_evidence(0, 0, 0, 100)
        assert code is None
    
    def test_case_control_zero_total_controls(self):
        """Zero total controls should return None."""
        code, _ = map_case_control_to_evidence(5, 50, 0, 0)
        assert code is None


if __name__ == '__main__':
    pytest.main([__file__, '-v'])

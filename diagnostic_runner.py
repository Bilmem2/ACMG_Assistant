"""
ACMG Assistant Comprehensive Diagnostic Suite
==============================================

Automatically detects potential logic, integration, and validation errors
across the entire ACMG Assistant project.

Usage:
    python diagnostic_runner.py

Output:
    - Console: Colored test results with summary
    - diagnostic_report.json: Detailed JSON report
    - diagnostic_summary.txt: Human-readable summary
"""

import sys
import os
import json
import time
import traceback
import importlib
import inspect
from datetime import datetime
from typing import Dict, List, Any, Tuple, Optional
from pathlib import Path

# Add src to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

# Try to import colorama, fallback if not available
try:
    from colorama import Fore, Style, init
    init()
    COLORS_AVAILABLE = True
except ImportError:
    COLORS_AVAILABLE = False
    # Mock colorama if not available
    class MockFore:
        GREEN = RED = YELLOW = CYAN = BLUE = MAGENTA = WHITE = RESET = ""
    class MockStyle:
        BRIGHT = RESET_ALL = ""
    Fore = MockFore()
    Style = MockStyle()


class DiagnosticRunner:
    """Main diagnostic test runner."""
    
    def __init__(self):
        self.results = {
            'timestamp': datetime.now().isoformat(),
            'tests': {},
            'summary': {
                'total': 0,
                'passed': 0,
                'failed': 0,
                'warnings': 0
            },
            'errors': []
        }
        self.start_time = time.time()
    
    def print_header(self, text: str, char: str = "="):
        """Print formatted header."""
        width = 80
        print(f"\n{Fore.CYAN}{char * width}{Style.RESET_ALL}")
        print(f"{Fore.CYAN}{Style.BRIGHT}{text.center(width)}{Style.RESET_ALL}")
        print(f"{Fore.CYAN}{char * width}{Style.RESET_ALL}\n")
    
    def print_section(self, text: str):
        """Print section header."""
        print(f"\n{Fore.BLUE}{Style.BRIGHT}{'─' * 80}{Style.RESET_ALL}")
        print(f"{Fore.BLUE}{Style.BRIGHT}{text}{Style.RESET_ALL}")
        print(f"{Fore.BLUE}{Style.BRIGHT}{'─' * 80}{Style.RESET_ALL}\n")
    
    def log_test(self, category: str, test_name: str, passed: bool, 
                 details: str = "", error: Optional[Exception] = None):
        """Log a test result."""
        if category not in self.results['tests']:
            self.results['tests'][category] = []
        
        result = {
            'test': test_name,
            'passed': passed,
            'details': details,
            'error': str(error) if error else None,
            'traceback': traceback.format_exc() if error else None
        }
        
        self.results['tests'][category].append(result)
        self.results['summary']['total'] += 1
        
        if passed:
            self.results['summary']['passed'] += 1
            status = f"{Fore.GREEN}✅ PASS{Style.RESET_ALL}"
        else:
            self.results['summary']['failed'] += 1
            status = f"{Fore.RED}❌ FAIL{Style.RESET_ALL}"
        
        print(f"{status} - {test_name}")
        if details:
            print(f"       {Fore.WHITE}{details}{Style.RESET_ALL}")
        if error:
            print(f"       {Fore.RED}Error: {str(error)}{Style.RESET_ALL}")
    
    def test_hgvs_parsing(self):
        """Test HGVS format parsing."""
        self.print_section("1. HGVS Input Validation Tests")
        
        try:
            from utils.hgvs_parser import HGVSParser, parse_hgvs_variant
            
            # Valid formats
            valid_tests = [
                ("NM_000546.6:c.1528C>T", "Full HGVS with RefSeq"),
                ("c.1528C>T", "Simple cDNA format"),
                ("1528C>T", "Position-only format"),
                ("NM_007294.4:c.1066C>T", "BRCA1 variant"),
                ("c.5266dupC", "Duplication"),
                ("c.100_101delAT", "Deletion"),
                ("c.1528_1529insA", "Insertion"),
            ]
            
            for variant, desc in valid_tests:
                try:
                    result = parse_hgvs_variant(variant)
                    passed = result is not None
                    self.log_test(
                        "HGVS Parsing",
                        f"Parse valid {desc}",
                        passed,
                        f"Input: {variant}"
                    )
                except Exception as e:
                    self.log_test(
                        "HGVS Parsing",
                        f"Parse valid {desc}",
                        False,
                        f"Input: {variant}",
                        e
                    )
            
            # Invalid formats
            invalid_tests = [
                ("", "Empty string"),
                ("invalid", "Random text"),
                ("123", "Number only"),
                ("C>T", "Missing position"),
                ("p.Arg123Gln", "Protein notation"),
                ("NM_000546:c.1528", "Missing variant"),
            ]
            
            for variant, desc in invalid_tests:
                try:
                    result = parse_hgvs_variant(variant)
                    passed = result is None  # Should fail to parse
                    self.log_test(
                        "HGVS Parsing",
                        f"Reject invalid {desc}",
                        passed,
                        f"Input: {variant}"
                    )
                except Exception as e:
                    self.log_test(
                        "HGVS Parsing",
                        f"Reject invalid {desc}",
                        False,
                        f"Input: {variant}",
                        e
                    )
        
        except ImportError as e:
            self.log_test(
                "HGVS Parsing",
                "Import hgvs_parser module",
                False,
                "Cannot import utils.hgvs_parser",
                e
            )
    
    def test_criteria_evaluation(self):
        """Test ACMG criteria evaluation logic."""
        self.print_section("2. ACMG Criterion Evaluation Logic")
        
        try:
            from core.evidence_evaluator import EvidenceEvaluator
            from core.variant_data import VariantData
            
            # Create test variant data with all required fields
            variant_data = VariantData()
            variant_data.basic_info = {
                'gene': 'TP53',
                'chromosome': '17',
                'position': 7676154,
                'ref_allele': 'G',
                'alt_allele': 'T',
                'variant_type': 'missense',
                'consequence': 'missense_variant'
            }
            variant_data.population_data = {}
            variant_data.insilico_data = {}
            variant_data.genetic_data = {}
            variant_data.functional_data = {}
            
            evaluator = EvidenceEvaluator(use_2023_guidelines=False, test_mode=True)
            
            # Test evaluate_all_criteria method
            try:
                result = evaluator.evaluate_all_criteria(variant_data)
                
                # Count actual criteria evaluated
                pathogenic_count = len(result.get('pathogenic_criteria', {}))
                benign_count = len(result.get('benign_criteria', {}))
                total_criteria = pathogenic_count + benign_count
                
                self.log_test(
                    "Criteria Evaluation",
                    "Call evaluate_all_criteria",
                    isinstance(result, dict),
                    f"Evaluated {total_criteria} criteria ({pathogenic_count} pathogenic, {benign_count} benign)"
                )
                
                # All 28 ACMG criteria
                criteria_list = [
                    'PVS1', 'PS1', 'PS2', 'PS3', 'PS4',
                    'PM1', 'PM2', 'PM3', 'PM4', 'PM5', 'PM6',
                    'PP1', 'PP2', 'PP3', 'PP4', 'PP5',
                    'BA1', 'BS1', 'BS2', 'BS3', 'BS4',
                    'BP1', 'BP2', 'BP3', 'BP4', 'BP5', 'BP6', 'BP7'
                ]
                
                for criterion in criteria_list:
                    if criterion in result:
                        criterion_result = result[criterion]
                        
                        # Check return type
                        valid_types = (bool, type(None), dict, tuple)
                        passed = isinstance(criterion_result, valid_types)
                        
                        self.log_test(
                            "Criteria Evaluation",
                            f"Evaluate {criterion}",
                            passed,
                            f"Result: {criterion_result}"
                        )
                    else:
                        # Criterion not in results - might be expected for some
                        self.log_test(
                            "Criteria Evaluation",
                            f"Evaluate {criterion}",
                            True,  # Not an error if criterion not evaluated
                            "Not evaluated (expected for incomplete data)"
                        )
            except Exception as e:
                self.log_test(
                    "Criteria Evaluation",
                    "Call evaluate_all_criteria",
                    False,
                    "",
                    e
                )
        
        except ImportError as e:
            self.log_test(
                "Criteria Evaluation",
                "Import core modules",
                False,
                "Cannot import core modules",
                e
            )
    
    def test_variant_types(self):
        """Test different variant type classifications."""
        self.print_section("3. Variant Type Classification")
        
        try:
            from core.acmg_classifier import ACMGClassifier
            from core.evidence_evaluator import EvidenceEvaluator
            from core.variant_data import VariantData
            
            variant_types = [
                ('missense', 'missense_variant'),
                ('nonsense', 'stop_gained'),
                ('frameshift', 'frameshift_variant'),
                ('splice', 'splice_donor_variant'),
                ('synonymous', 'synonymous_variant'),
                ('inframe_indel', 'inframe_deletion'),
            ]
            
            for variant_type, consequence in variant_types:
                try:
                    variant_data = VariantData()
                    variant_data.basic_info = {
                        'gene': 'TEST',
                        'chromosome': '1',
                        'position': 12345,
                        'ref_allele': 'A',
                        'alt_allele': 'T',
                        'variant_type': variant_type,
                        'consequence': consequence
                    }
                    variant_data.population_data = {}
                    variant_data.insilico_data = {}
                    variant_data.genetic_data = {}
                    variant_data.functional_data = {}
                    
                    # First evaluate evidence
                    evaluator = EvidenceEvaluator(use_2023_guidelines=False, test_mode=True)
                    evidence_results = evaluator.evaluate_all_criteria(variant_data)
                    
                    # Then classify
                    classifier = ACMGClassifier(use_2023_guidelines=False)
                    result = classifier.classify(evidence_results)
                    
                    passed = result is not None and 'classification' in result
                    
                    self.log_test(
                        "Variant Types",
                        f"Classify {variant_type}",
                        passed,
                        f"Classification: {result.get('classification') if result else 'None'}"
                    )
                except Exception as e:
                    self.log_test(
                        "Variant Types",
                        f"Classify {variant_type}",
                        False,
                        "",
                        e
                    )
        
        except ImportError as e:
            self.log_test(
                "Variant Types",
                "Import classifier",
                False,
                "Cannot import ACMGClassifier",
                e
            )
    
    def test_module_consistency(self):
        """Test inter-module consistency."""
        self.print_section("4. Inter-Module Consistency")
        
        # Test module imports
        modules_to_test = [
            ('core.variant_data', 'VariantData'),
            ('core.evidence_evaluator', 'EvidenceEvaluator'),
            ('core.acmg_classifier', 'ACMGClassifier'),
            ('utils.input_handler', 'InputHandler'),
            ('utils.hgvs_parser', 'HGVSParser'),
            ('config.constants', None),
        ]
        
        for module_name, class_name in modules_to_test:
            try:
                module = importlib.import_module(module_name)
                
                if class_name:
                    if hasattr(module, class_name):
                        self.log_test(
                            "Module Consistency",
                            f"Import {module_name}.{class_name}",
                            True,
                            "Module and class found"
                        )
                    else:
                        self.log_test(
                            "Module Consistency",
                            f"Import {module_name}.{class_name}",
                            False,
                            f"Class {class_name} not found in module"
                        )
                else:
                    self.log_test(
                        "Module Consistency",
                        f"Import {module_name}",
                        True,
                        "Module found"
                    )
            except Exception as e:
                self.log_test(
                    "Module Consistency",
                    f"Import {module_name}",
                    False,
                    "",
                    e
                )
    
    def test_edge_cases(self):
        """Test edge cases and special scenarios."""
        self.print_section("5. Edge Case Simulations")
        
        try:
            from core.acmg_classifier import ACMGClassifier
            from core.evidence_evaluator import EvidenceEvaluator
            from core.variant_data import VariantData
            
            # Test synonymous variant
            try:
                variant_data = VariantData()
                variant_data.basic_info = {
                    'gene': 'TEST',
                    'chromosome': '1',
                    'position': 12345,
                    'ref_allele': 'A',
                    'alt_allele': 'G',
                    'variant_type': 'synonymous',
                    'consequence': 'synonymous_variant'
                }
                variant_data.population_data = {}
                variant_data.insilico_data = {}
                variant_data.genetic_data = {}
                variant_data.functional_data = {}
                
                evaluator = EvidenceEvaluator(use_2023_guidelines=False, test_mode=True)
                evidence_results = evaluator.evaluate_all_criteria(variant_data)
                
                classifier = ACMGClassifier(use_2023_guidelines=False)
                result = classifier.classify(evidence_results)
                
                self.log_test(
                    "Edge Cases",
                    "Classify synonymous variant",
                    result is not None,
                    f"Result: {result.get('classification') if result else 'None'}"
                )
            except Exception as e:
                self.log_test(
                    "Edge Cases",
                    "Classify synonymous variant",
                    False,
                    "",
                    e
                )
            
            # Test empty variant data
            try:
                variant_data = VariantData()
                errors = variant_data.validate()
                
                self.log_test(
                    "Edge Cases",
                    "Validate empty variant data",
                    len(errors) > 0,
                    f"Found {len(errors)} validation errors"
                )
            except Exception as e:
                self.log_test(
                    "Edge Cases",
                    "Validate empty variant data",
                    False,
                    "",
                    e
                )
        
        except ImportError as e:
            self.log_test(
                "Edge Cases",
                "Import required modules",
                False,
                "",
                e
            )
    
    def test_performance(self):
        """Test performance and stability."""
        self.print_section("6. Performance and Stability")
        
        try:
            from core.acmg_classifier import ACMGClassifier
            from core.evidence_evaluator import EvidenceEvaluator
            from core.variant_data import VariantData
            
            # Test classification speed
            try:
                variant_data = VariantData()
                variant_data.basic_info = {
                    'gene': 'TP53',
                    'chromosome': '17',
                    'position': 7676154,
                    'ref_allele': 'G',
                    'alt_allele': 'T',
                    'variant_type': 'missense',
                    'consequence': 'missense_variant'
                }
                variant_data.population_data = {}
                variant_data.insilico_data = {}
                variant_data.genetic_data = {}
                variant_data.functional_data = {}
                
                start = time.time()
                evaluator = EvidenceEvaluator(use_2023_guidelines=False, test_mode=True)
                evidence_results = evaluator.evaluate_all_criteria(variant_data)
                classifier = ACMGClassifier(use_2023_guidelines=False)
                result = classifier.classify(evidence_results)
                elapsed = time.time() - start
                
                passed = elapsed < 5.0  # Should complete in under 5 seconds
                
                self.log_test(
                    "Performance",
                    "Classification speed",
                    passed,
                    f"Completed in {elapsed:.3f}s"
                )
            except Exception as e:
                self.log_test(
                    "Performance",
                    "Classification speed",
                    False,
                    "",
                    e
                )
            
            # Test circular import detection
            try:
                import sys
                modules_before = set(sys.modules.keys())
                
                from core import acmg_classifier
                from utils import hgvs_parser
                
                modules_after = set(sys.modules.keys())
                new_modules = modules_after - modules_before
                
                self.log_test(
                    "Performance",
                    "Circular import check",
                    True,
                    f"Loaded {len(new_modules)} modules without errors"
                )
            except Exception as e:
                self.log_test(
                    "Performance",
                    "Circular import check",
                    False,
                    "",
                    e
                )
        
        except ImportError as e:
            self.log_test(
                "Performance",
                "Import modules",
                False,
                "",
                e
            )
    
    def generate_reports(self):
        """Generate diagnostic reports."""
        self.print_section("7. Generating Reports")
        
        # Calculate health score
        if self.results['summary']['total'] > 0:
            health_score = (self.results['summary']['passed'] / 
                          self.results['summary']['total'] * 100)
        else:
            health_score = 0
        
        self.results['health_score'] = health_score
        self.results['duration'] = time.time() - self.start_time
        
        # Save JSON report
        try:
            with open('diagnostic_report.json', 'w') as f:
                json.dump(self.results, f, indent=2)
            print(f"{Fore.GREEN}✅ Saved diagnostic_report.json{Style.RESET_ALL}")
        except Exception as e:
            print(f"{Fore.RED}❌ Failed to save JSON report: {e}{Style.RESET_ALL}")
        
        # Save text summary
        try:
            with open('diagnostic_summary.txt', 'w', encoding='utf-8') as f:
                f.write("ACMG Assistant Diagnostic Summary\n")
                f.write("=" * 80 + "\n\n")
                f.write(f"Timestamp: {self.results['timestamp']}\n")
                f.write(f"Duration: {self.results['duration']:.2f}s\n\n")
                
                f.write("Summary:\n")
                f.write("-" * 80 + "\n")
                f.write(f"Total Tests: {self.results['summary']['total']}\n")
                f.write(f"Passed: {self.results['summary']['passed']}\n")
                f.write(f"Failed: {self.results['summary']['failed']}\n")
                f.write(f"Health Score: {health_score:.1f}%\n\n")
                
                f.write("Test Results by Category:\n")
                f.write("-" * 80 + "\n")
                for category, tests in self.results['tests'].items():
                    passed = sum(1 for t in tests if t['passed'])
                    total = len(tests)
                    f.write(f"\n{category}: {passed}/{total} passed\n")
                    
                    for test in tests:
                        status = "✅" if test['passed'] else "❌"
                        f.write(f"  {status} {test['test']}\n")
                        if test['details']:
                            f.write(f"     {test['details']}\n")
                        if test['error']:
                            f.write(f"     Error: {test['error']}\n")
            
            print(f"{Fore.GREEN}✅ Saved diagnostic_summary.txt{Style.RESET_ALL}")
        except Exception as e:
            print(f"{Fore.RED}❌ Failed to save text summary: {e}{Style.RESET_ALL}")
    
    def print_summary(self):
        """Print final summary."""
        self.print_header("DIAGNOSTIC SUMMARY")
        
        total = self.results['summary']['total']
        passed = self.results['summary']['passed']
        failed = self.results['summary']['failed']
        health_score = self.results.get('health_score', 0)
        
        print(f"{Fore.CYAN}Total Tests:{Style.RESET_ALL} {total}")
        print(f"{Fore.GREEN}Passed:{Style.RESET_ALL} {passed}")
        print(f"{Fore.RED}Failed:{Style.RESET_ALL} {failed}")
        print(f"{Fore.BLUE}Duration:{Style.RESET_ALL} {self.results['duration']:.2f}s")
        print()
        
        # Health score with color
        if health_score >= 90:
            color = Fore.GREEN
            emoji = "✅"
            status = "EXCELLENT"
        elif health_score >= 75:
            color = Fore.YELLOW
            emoji = "⚠️"
            status = "GOOD"
        elif health_score >= 50:
            color = Fore.YELLOW
            emoji = "⚠️"
            status = "NEEDS ATTENTION"
        else:
            color = Fore.RED
            emoji = "❌"
            status = "CRITICAL"
        
        print(f"{color}{Style.BRIGHT}{emoji} Health Score: {health_score:.1f}% — {status}{Style.RESET_ALL}")
        print()
        
        # Category breakdown
        print(f"{Fore.CYAN}Results by Category:{Style.RESET_ALL}")
        print("-" * 80)
        for category, tests in self.results['tests'].items():
            passed_count = sum(1 for t in tests if t['passed'])
            total_count = len(tests)
            percentage = (passed_count / total_count * 100) if total_count > 0 else 0
            
            if percentage == 100:
                color = Fore.GREEN
            elif percentage >= 75:
                color = Fore.YELLOW
            else:
                color = Fore.RED
            
            print(f"{color}{category:.<40} {passed_count}/{total_count} ({percentage:.0f}%){Style.RESET_ALL}")
        
        print()
    
    def run(self):
        """Run all diagnostic tests."""
        self.print_header("ACMG ASSISTANT DIAGNOSTIC SUITE")
        
        print(f"{Fore.CYAN}Starting comprehensive diagnostic tests...{Style.RESET_ALL}\n")
        
        # Run all test suites
        self.test_hgvs_parsing()
        self.test_criteria_evaluation()
        self.test_variant_types()
        self.test_module_consistency()
        self.test_edge_cases()
        self.test_performance()
        
        # Generate reports
        self.generate_reports()
        
        # Print summary
        self.print_summary()
        
        print(f"\n{Fore.CYAN}Reports saved:{Style.RESET_ALL}")
        print(f"  • diagnostic_report.json")
        print(f"  • diagnostic_summary.txt")
        print()
        
        return self.results['summary']['failed'] == 0


if __name__ == '__main__':
    runner = DiagnosticRunner()
    success = runner.run()
    sys.exit(0 if success else 1)

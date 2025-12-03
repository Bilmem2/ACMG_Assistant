"""Core modules for ACMG Assistant."""
# Makes this folder a Python package

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

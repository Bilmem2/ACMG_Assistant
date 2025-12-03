"""
Gene-Specific Rules for ACMG Criteria Evaluation
=================================================

This module provides gene-specific rule evaluation, particularly for PM1
(mutational hotspot / critical functional domain) evidence.

DESIGN PRINCIPLE:
    This module MUST NOT contain hardcoded biological data (gene names,
    hotspot positions, domain coordinates). All biological facts are
    fetched from remote APIs via DomainAPIClient. Local configuration
    is limited to interpretive thresholds.

Author: Can Sevilmiş
License: MIT License
"""

from dataclasses import dataclass
from typing import Dict, Any, Optional
import re


@dataclass
class PM1Evidence:
    """
    Result of PM1 evidence evaluation.
    
    Attributes:
        applies: Whether PM1 evidence applies
        evidence_code: 'PM1' for full strength, 'PM1_supporting' for reduced
        confidence: Float 0.0-1.0 indicating evidence strength
        source: Data source(s) used for evaluation
        reason: Human-readable explanation
    """
    applies: bool = False
    evidence_code: str = ""
    confidence: float = 0.0
    source: str = ""
    reason: str = ""
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for serialization."""
        return {
            'applies': self.applies,
            'evidence_code': self.evidence_code,
            'confidence': self.confidence,
            'source': self.source,
            'reason': self.reason,
        }


class GeneSpecificRules:
    """
    Evaluates gene-specific ACMG criteria using remote API data.
    
    This class acts as an adapter between the raw API responses from
    DomainAPIClient and the ACMG evidence system. It applies interpretive
    thresholds to convert hotspot/domain annotations into PM1 evidence.
    
    IMPORTANT: No hardcoded biological data. All gene/position-specific
    information comes from DomainAPIClient which fetches from remote APIs.
    
    Usage:
        rules = GeneSpecificRules()
        pm1 = rules.evaluate_pm1(gene="TP53", position=248)
        if pm1.applies:
            print(f"{pm1.evidence_code}: {pm1.reason}")
    """
    
    # Interpretive thresholds (not biological data)
    # These define how to convert API confidence scores to ACMG evidence
    PM1_THRESHOLDS = {
        'PM1': 0.85,           # High confidence -> full PM1
        'PM1_supporting': 0.60  # Moderate confidence -> PM1_supporting
    }
    
    def __init__(self, domain_client=None):
        """
        Initialize gene-specific rules evaluator.
        
        Args:
            domain_client: Optional DomainAPIClient instance. If not provided,
                          one will be created on first use.
        """
        self._domain_client = domain_client
    
    @property
    def domain_client(self):
        """Lazy-load DomainAPIClient to avoid circular imports."""
        if self._domain_client is None:
            from utils.domain_api_client import DomainAPIClient
            self._domain_client = DomainAPIClient(cache_enabled=True)
        return self._domain_client
    
    def evaluate_pm1(
        self,
        gene: str,
        position: Optional[int] = None,
        hgvs_p: Optional[str] = None
    ) -> PM1Evidence:
        """
        Evaluate PM1 evidence for a variant position.
        
        PM1: Located in a mutational hot spot and/or critical and well-
        established functional domain (e.g., active site of an enzyme)
        without benign variation.
        
        This method queries remote APIs via DomainAPIClient and applies
        local thresholds to determine if PM1 evidence applies.
        
        Args:
            gene: Gene symbol (e.g., 'TP53')
            position: Amino acid position (e.g., 248)
            hgvs_p: Protein HGVS notation (e.g., 'p.Arg248Gln')
        
        Returns:
            PM1Evidence with applies, evidence_code, confidence, and reason
        """
        # Extract position from hgvs_p if not provided
        if position is None and hgvs_p:
            position = self._extract_position(hgvs_p)
        
        # Cannot evaluate PM1 without gene
        if not gene:
            return PM1Evidence(
                applies=False,
                reason="Gene symbol required for PM1 evaluation"
            )
        
        # Query remote APIs for hotspot/domain annotation
        try:
            annotation = self.domain_client.get_hotspot_annotation(
                gene=gene,
                position=position,
                hgvs_p=hgvs_p
            )
        except Exception as e:
            # API error - return negative result (safe fallback)
            return PM1Evidence(
                applies=False,
                source="error",
                reason=f"API query failed: {str(e)[:100]}"
            )
        
        # No remote data available
        if annotation.source == "none":
            return PM1Evidence(
                applies=False,
                source="none",
                reason="No hotspot/domain data available from remote APIs"
            )
        
        # Apply thresholds to determine evidence level
        return self._apply_pm1_thresholds(annotation, gene, position)
    
    def _apply_pm1_thresholds(
        self,
        annotation,  # HotspotAnnotation
        gene: str,
        position: Optional[int]
    ) -> PM1Evidence:
        """
        Apply interpretive thresholds to convert annotation to PM1 evidence.
        
        Args:
            annotation: HotspotAnnotation from DomainAPIClient
            gene: Gene symbol
            position: Amino acid position
        
        Returns:
            PM1Evidence based on threshold application
        """
        confidence = annotation.confidence
        
        # Hotspot or critical domain with high confidence -> PM1
        if confidence >= self.PM1_THRESHOLDS['PM1']:
            if annotation.is_hotspot:
                reason = (
                    f"Position {position} in {gene} is a well-documented "
                    f"mutational hotspot ({annotation.tumor_count} tumors)"
                )
            elif annotation.in_critical_domain:
                reason = (
                    f"Position {position} in {gene} is within critical "
                    f"functional domain: {annotation.domain_name}"
                )
            else:
                reason = f"Strong hotspot/domain evidence for {gene} position {position}"
            
            return PM1Evidence(
                applies=True,
                evidence_code='PM1',
                confidence=confidence,
                source=annotation.source,
                reason=reason
            )
        
        # Moderate confidence -> PM1_supporting
        if confidence >= self.PM1_THRESHOLDS['PM1_supporting']:
            if annotation.is_hotspot:
                reason = (
                    f"Position {position} in {gene} shows moderate hotspot "
                    f"evidence ({annotation.tumor_count} tumors)"
                )
            elif annotation.in_critical_domain:
                reason = (
                    f"Position {position} in {gene} is within functional "
                    f"region: {annotation.domain_name}"
                )
            else:
                reason = f"Moderate hotspot/domain evidence for {gene} position {position}"
            
            return PM1Evidence(
                applies=True,
                evidence_code='PM1_supporting',
                confidence=confidence,
                source=annotation.source,
                reason=reason
            )
        
        # Low confidence or no evidence -> PM1 does not apply
        return PM1Evidence(
            applies=False,
            confidence=confidence,
            source=annotation.source,
            reason=annotation.details or f"Insufficient evidence for PM1 in {gene}"
        )
    
    def _extract_position(self, hgvs_p: Optional[str]) -> Optional[int]:
        """
        Extract amino acid position from protein HGVS notation.
        
        Args:
            hgvs_p: Protein HGVS string (e.g., 'p.Arg248Gln', 'p.R248Q')
        
        Returns:
            Integer position or None if extraction fails
        """
        if not hgvs_p:
            return None
        
        # Match patterns like:
        # - p.Arg248Gln (three-letter code)
        # - p.R248Q (single-letter code)
        # - p.Arg248* (nonsense)
        match = re.search(r'p\.(?:[A-Za-z]{1,3})?(\d+)', hgvs_p)
        if match:
            return int(match.group(1))
        return None

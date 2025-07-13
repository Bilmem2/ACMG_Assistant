#!/usr/bin/env python3
"""
ACMG Test Pipeline
==================

Comprehensive test pipeline for ACMG variant classification system.
This script tests the classifier with various variant types and evaluates accuracy.

Author: Can Sevilmiş
Version: 3.1.0
Last Updated: July 12, 2025
"""

import pandas as pd
import os
import sys
import logging
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime
import traceback

# ACMG Assistant import için path ekleme
SRC_DIR = os.path.abspath("src")
if os.path.exists(SRC_DIR):
    sys.path.insert(0, SRC_DIR)
    try:
        from acmg_assistant import ACMGAssistant
        from core.variant_data import VariantData
        from core.evidence_evaluator import EvidenceEvaluator
        from core.acmg_classifier import ACMGClassifier
        print("✅ ACMG modules successfully imported")
    except ImportError as e:
        print(f"❌ HATA: Module import hatası: {e}")
        sys.exit(1)
else:
    print(f"❌ HATA: {SRC_DIR} dizini bulunamadı!")
    sys.exit(1)

# Loglama ayarları
LOG_DIR = "test_results"
LOG_FILE = os.path.join(LOG_DIR, f"acmg_test_pipeline_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log")
os.makedirs(LOG_DIR, exist_ok=True)

logging.basicConfig(
    filename=LOG_FILE,
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logging.info("ACMG Test Pipeline başlatılıyor...")

print("🧬 ACMG Variant Classification Assistant - Test Pipeline")
print("=" * 60)

# =============================
# 🔧 1. Kapsamlı Test Dataset Hazırlığı
# =============================
# Gerçek varyant bilgileri ile her ACMG kriteri için test senaryoları
test_variants = [
    # ===== PVS1 Test Senaryoları (Null variants) =====
    {
        "variant_id": "BRCA1_c.68_69delAG",
        "gene": "BRCA1",
        "hgvs_c": "c.68_69delAG",
        "hgvs_p": "p.Glu23Valfs*17",
        "consequence": "frameshift_variant",
        "gnomad_af": 0.00001,
        "revel_score": 0.95,
        "cadd_phred": 35,
        "alphamissense": 0.9,
        "sift_score": 0.01,
        "polyphen2_score": 0.99,
        "spliceai_ag_score": 0.02,
        "spliceai_al_score": 0.01,
        "spliceai_dg_score": 0.00,
        "spliceai_dl_score": 0.00,
        "clinvar_significance": "Pathogenic",
        "expected_class": "Pathogenic",
        "acmg_criteria": "PVS1"
    },
    {
        "variant_id": "TP53_c.1024C>T",
        "gene": "TP53",
        "hgvs_c": "c.1024C>T",
        "hgvs_p": "p.Arg342Ter",
        "consequence": "stop_gained",
        "gnomad_af": 0.000001,
        "revel_score": 0.99,
        "cadd_phred": 40,
        "alphamissense": 0.95,
        "sift_score": 0.00,
        "polyphen2_score": 1.00,
        "spliceai_ag_score": 0.00,
        "spliceai_al_score": 0.00,
        "spliceai_dg_score": 0.00,
        "spliceai_dl_score": 0.00,
        "clinvar_significance": "Pathogenic",
        "expected_class": "Pathogenic",
        "acmg_criteria": "PVS1"
    },
    {
        "variant_id": "BRCA2_c.5946delT",
        "gene": "BRCA2",
        "hgvs_c": "c.5946delT",
        "hgvs_p": "p.Ser1982Argfs*22",
        "consequence": "frameshift_variant",
        "gnomad_af": 0.00001,
        "revel_score": 0.97,
        "cadd_phred": 37,
        "alphamissense": 0.92,
        "sift_score": 0.01,
        "polyphen2_score": 0.99,
        "spliceai_ag_score": 0.03,
        "spliceai_al_score": 0.02,
        "spliceai_dg_score": 0.01,
        "spliceai_dl_score": 0.00,
        "clinvar_significance": "Pathogenic",
        "expected_class": "Pathogenic",
        "acmg_criteria": "PVS1"
    },
    {
        "variant_id": "MLH1_c.1852_1854delAAG",
        "gene": "MLH1",
        "hgvs_c": "c.1852_1854delAAG",
        "hgvs_p": "p.Lys618del",
        "consequence": "inframe_deletion",
        "gnomad_af": 0.000005,
        "revel_score": 0.91,
        "cadd_phred": 34,
        "alphamissense": 0.87,
        "sift_score": 0.01,
        "polyphen2_score": 0.97,
        "spliceai_ag_score": 0.01,
        "spliceai_al_score": 0.00,
        "spliceai_dg_score": 0.00,
        "spliceai_dl_score": 0.00,
        "clinvar_significance": "Pathogenic",
        "expected_class": "Pathogenic",
        "acmg_criteria": "PVS1"
    },
    
    # ===== PS1 Test Senaryoları (Same amino acid change) =====
    {
        "variant_id": "TP53_c.743G>A",
        "gene": "TP53",
        "hgvs_c": "c.743G>A",
        "hgvs_p": "p.Arg248Gln",
        "consequence": "missense_variant",
        "gnomad_af": 0.000005,
        "revel_score": 0.85,
        "cadd_phred": 28,
        "alphamissense": 0.78,
        "sift_score": 0.02,
        "polyphen2_score": 0.96,
        "spliceai_ag_score": 0.01,
        "spliceai_al_score": 0.00,
        "spliceai_dg_score": 0.00,
        "spliceai_dl_score": 0.00,
        "clinvar_significance": "Pathogenic",
        "expected_class": "Pathogenic",
        "acmg_criteria": "PS1"
    },
    {
        "variant_id": "BRCA1_c.5266dupC",
        "gene": "BRCA1",
        "hgvs_c": "c.5266dupC",
        "hgvs_p": "p.Gln1756Profs*74",
        "consequence": "frameshift_variant",
        "gnomad_af": 0.000008,
        "revel_score": 0.96,
        "cadd_phred": 36,
        "alphamissense": 0.91,
        "sift_score": 0.01,
        "polyphen2_score": 0.99,
        "spliceai_ag_score": 0.02,
        "spliceai_al_score": 0.01,
        "spliceai_dg_score": 0.00,
        "spliceai_dl_score": 0.00,
        "clinvar_significance": "Pathogenic",
        "expected_class": "Pathogenic",
        "acmg_criteria": "PS1"
    },
    
    # ===== PS2 Test Senaryoları (De novo variants) =====
    {
        "variant_id": "MECP2_c.763C>T",
        "gene": "MECP2",
        "hgvs_c": "c.763C>T",
        "hgvs_p": "p.Arg255Ter",
        "consequence": "stop_gained",
        "gnomad_af": 0.000001,
        "revel_score": 0.98,
        "cadd_phred": 39,
        "alphamissense": 0.94,
        "sift_score": 0.00,
        "polyphen2_score": 1.00,
        "spliceai_ag_score": 0.00,
        "spliceai_al_score": 0.00,
        "spliceai_dg_score": 0.00,
        "spliceai_dl_score": 0.00,
        "clinvar_significance": "Pathogenic",
        "expected_class": "Pathogenic",
        "acmg_criteria": "PS2"
    },
    {
        "variant_id": "SCN1A_c.2386C>T",
        "gene": "SCN1A",
        "hgvs_c": "c.2386C>T",
        "hgvs_p": "p.Arg796Ter",
        "consequence": "stop_gained",
        "gnomad_af": 0.000001,
        "revel_score": 0.98,
        "cadd_phred": 38,
        "alphamissense": 0.94,
        "sift_score": 0.00,
        "polyphen2_score": 1.00,
        "spliceai_ag_score": 0.00,
        "spliceai_al_score": 0.00,
        "spliceai_dg_score": 0.00,
        "spliceai_dl_score": 0.00,
        "clinvar_significance": "Pathogenic",
        "expected_class": "Pathogenic",
        "acmg_criteria": "PS2"
    },
    
    # ===== PS3 Test Senaryoları (Functional studies) =====
    {
        "variant_id": "BRCA1_c.181T>G",
        "gene": "BRCA1",
        "hgvs_c": "c.181T>G",
        "hgvs_p": "p.Cys61Gly",
        "consequence": "missense_variant",
        "gnomad_af": 0.000002,
        "revel_score": 0.89,
        "cadd_phred": 33,
        "alphamissense": 0.86,
        "sift_score": 0.01,
        "polyphen2_score": 0.98,
        "spliceai_ag_score": 0.01,
        "spliceai_al_score": 0.00,
        "spliceai_dg_score": 0.00,
        "spliceai_dl_score": 0.00,
        "clinvar_significance": "Pathogenic",
        "expected_class": "Pathogenic",
        "acmg_criteria": "PS3"
    },
    
    # ===== PS4 Test Senaryoları (Case-control studies) =====
    {
        "variant_id": "PALB2_c.1592delT",
        "gene": "PALB2",
        "hgvs_c": "c.1592delT",
        "hgvs_p": "p.Leu531Argfs*46",
        "consequence": "frameshift_variant",
        "gnomad_af": 0.000003,
        "revel_score": 0.94,
        "cadd_phred": 35,
        "alphamissense": 0.90,
        "sift_score": 0.01,
        "polyphen2_score": 0.99,
        "spliceai_ag_score": 0.02,
        "spliceai_al_score": 0.01,
        "spliceai_dg_score": 0.00,
        "spliceai_dl_score": 0.00,
        "clinvar_significance": "Pathogenic",
        "expected_class": "Pathogenic",
        "acmg_criteria": "PS4"
    },
    
    # ===== PM1 Test Senaryoları (Hotspot regions) =====
    {
        "variant_id": "KRAS_c.35G>A",
        "gene": "KRAS",
        "hgvs_c": "c.35G>A",
        "hgvs_p": "p.Gly12Asp",
        "consequence": "missense_variant",
        "gnomad_af": 0.000001,
        "revel_score": 0.92,
        "cadd_phred": 34,
        "alphamissense": 0.88,
        "sift_score": 0.01,
        "polyphen2_score": 0.98,
        "spliceai_ag_score": 0.00,
        "spliceai_al_score": 0.00,
        "spliceai_dg_score": 0.00,
        "spliceai_dl_score": 0.00,
        "clinvar_significance": "Pathogenic",
        "expected_class": "Pathogenic",
        "acmg_criteria": "PM1"
    },
    
    # ===== PM2 Test Senaryoları (Absent in controls) =====
    {
        "variant_id": "BRCA2_c.8755C>T",
        "gene": "BRCA2",
        "hgvs_c": "c.8755C>T",
        "hgvs_p": "p.Arg2919Ter",
        "consequence": "stop_gained",
        "gnomad_af": 0.0,  # Absent in gnomAD
        "revel_score": 0.99,
        "cadd_phred": 40,
        "alphamissense": 0.95,
        "sift_score": 0.00,
        "polyphen2_score": 1.00,
        "spliceai_ag_score": 0.00,
        "spliceai_al_score": 0.00,
        "spliceai_dg_score": 0.00,
        "spliceai_dl_score": 0.00,
        "clinvar_significance": "Pathogenic",
        "expected_class": "Pathogenic",
        "acmg_criteria": "PM2"
    },
    
    # ===== PM3 Test Senaryoları (Recessive disorders) =====
    {
        "variant_id": "CFTR_c.1040G>C",
        "gene": "CFTR",
        "hgvs_c": "c.1040G>C",
        "hgvs_p": "p.Arg347Pro",
        "consequence": "missense_variant",
        "gnomad_af": 0.000012,
        "revel_score": 0.87,
        "cadd_phred": 31,
        "alphamissense": 0.83,
        "sift_score": 0.02,
        "polyphen2_score": 0.97,
        "spliceai_ag_score": 0.01,
        "spliceai_al_score": 0.00,
        "spliceai_dg_score": 0.00,
        "spliceai_dl_score": 0.00,
        "clinvar_significance": "Pathogenic",
        "expected_class": "Pathogenic",
        "acmg_criteria": "PM3"
    },
    
    # ===== PM4 Test Senaryoları (Protein length changes) =====
    {
        "variant_id": "FBN1_c.1129_1134delGGAGGA",
        "gene": "FBN1",
        "hgvs_c": "c.1129_1134delGGAGGA",
        "hgvs_p": "p.Gly377_Arg378del",
        "consequence": "inframe_deletion",
        "gnomad_af": 0.000001,
        "revel_score": 0.90,
        "cadd_phred": 33,
        "alphamissense": 0.87,
        "sift_score": 0.01,
        "polyphen2_score": 0.97,
        "spliceai_ag_score": 0.01,
        "spliceai_al_score": 0.00,
        "spliceai_dg_score": 0.00,
        "spliceai_dl_score": 0.00,
        "clinvar_significance": "Pathogenic",
        "expected_class": "Pathogenic",
        "acmg_criteria": "PM4"
    },
    
    # ===== PM5 Test Senaryoları (Same position, different change) =====
    {
        "variant_id": "TP53_c.742C>T",
        "gene": "TP53",
        "hgvs_c": "c.742C>T",
        "hgvs_p": "p.Arg248Trp",
        "consequence": "missense_variant",
        "gnomad_af": 0.000003,
        "revel_score": 0.86,
        "cadd_phred": 29,
        "alphamissense": 0.81,
        "sift_score": 0.02,
        "polyphen2_score": 0.95,
        "spliceai_ag_score": 0.01,
        "spliceai_al_score": 0.00,
        "spliceai_dg_score": 0.00,
        "spliceai_dl_score": 0.00,
        "clinvar_significance": "Pathogenic",
        "expected_class": "Pathogenic",
        "acmg_criteria": "PM5"
    },
    
    # ===== PP1 Test Senaryoları (Segregation) =====
    {
        "variant_id": "BRCA1_c.4327C>T",
        "gene": "BRCA1",
        "hgvs_c": "c.4327C>T",
        "hgvs_p": "p.Arg1443Ter",
        "consequence": "stop_gained",
        "gnomad_af": 0.000002,
        "revel_score": 0.97,
        "cadd_phred": 37,
        "alphamissense": 0.93,
        "sift_score": 0.00,
        "polyphen2_score": 1.00,
        "spliceai_ag_score": 0.00,
        "spliceai_al_score": 0.00,
        "spliceai_dg_score": 0.00,
        "spliceai_dl_score": 0.00,
        "clinvar_significance": "Pathogenic",
        "expected_class": "Pathogenic",
        "acmg_criteria": "PP1"
    },
    
    # ===== PP2 Test Senaryoları (Missense rare in gene) =====
    {
        "variant_id": "PTEN_c.634G>A",
        "gene": "PTEN",
        "hgvs_c": "c.634G>A",
        "hgvs_p": "p.Ala212Thr",
        "consequence": "missense_variant",
        "gnomad_af": 0.000001,
        "revel_score": 0.84,
        "cadd_phred": 28,
        "alphamissense": 0.79,
        "sift_score": 0.02,
        "polyphen2_score": 0.94,
        "spliceai_ag_score": 0.01,
        "spliceai_al_score": 0.00,
        "spliceai_dg_score": 0.00,
        "spliceai_dl_score": 0.00,
        "clinvar_significance": "Pathogenic",
        "expected_class": "Pathogenic",
        "acmg_criteria": "PP2"
    },
    
    # ===== PP3 Test Senaryoları (In silico damaging) =====
    {
        "variant_id": "ATM_c.5932G>T",
        "gene": "ATM",
        "hgvs_c": "c.5932G>T",
        "hgvs_p": "p.Glu1978Ter",
        "consequence": "stop_gained",
        "gnomad_af": 0.000001,
        "revel_score": 0.98,
        "cadd_phred": 39,
        "alphamissense": 0.94,
        "sift_score": 0.00,
        "polyphen2_score": 1.00,
        "spliceai_ag_score": 0.00,
        "spliceai_al_score": 0.00,
        "spliceai_dg_score": 0.00,
        "spliceai_dl_score": 0.00,
        "clinvar_significance": "Pathogenic",
        "expected_class": "Pathogenic",
        "acmg_criteria": "PP3"
    },
    
    # ===== PP4 Test Senaryoları (Patient phenotype) =====
    {
        "variant_id": "RYR1_c.1840C>T",
        "gene": "RYR1",
        "hgvs_c": "c.1840C>T",
        "hgvs_p": "p.Arg614Cys",
        "consequence": "missense_variant",
        "gnomad_af": 0.000005,
        "revel_score": 0.88,
        "cadd_phred": 32,
        "alphamissense": 0.84,
        "sift_score": 0.01,
        "polyphen2_score": 0.96,
        "spliceai_ag_score": 0.01,
        "spliceai_al_score": 0.00,
        "spliceai_dg_score": 0.00,
        "spliceai_dl_score": 0.00,
        "clinvar_significance": "Pathogenic",
        "expected_class": "Pathogenic",
        "acmg_criteria": "PP4"
    },
    
    # ===== VUS Test Senaryoları =====
    {
        "variant_id": "TP53_c.1010G>A",
        "gene": "TP53",
        "hgvs_c": "c.1010G>A",
        "hgvs_p": "p.Arg337His",
        "consequence": "missense_variant",
        "gnomad_af": 0.00002,
        "revel_score": 0.45,
        "cadd_phred": 15,
        "alphamissense": 0.3,
        "sift_score": 0.4,
        "polyphen2_score": 0.5,
        "spliceai_ag_score": 0.01,
        "spliceai_al_score": 0.00,
        "spliceai_dg_score": 0.00,
        "spliceai_dl_score": 0.00,
        "clinvar_significance": "Uncertain significance",
        "expected_class": "VUS",
        "acmg_criteria": "VUS"
    },
    {
        "variant_id": "BRCA2_c.7435T>C",
        "gene": "BRCA2",
        "hgvs_c": "c.7435T>C",
        "hgvs_p": "p.Ser2479Pro",
        "consequence": "missense_variant",
        "gnomad_af": 0.00015,
        "revel_score": 0.52,
        "cadd_phred": 18,
        "alphamissense": 0.35,
        "sift_score": 0.35,
        "polyphen2_score": 0.45,
        "spliceai_ag_score": 0.02,
        "spliceai_al_score": 0.01,
        "spliceai_dg_score": 0.00,
        "spliceai_dl_score": 0.00,
        "clinvar_significance": "Uncertain significance",
        "expected_class": "VUS",
        "acmg_criteria": "VUS"
    },
    
    # ===== BA1 Test Senaryoları (Common variants) =====
    {
        "variant_id": "APOE_c.334T>C",
        "gene": "APOE",
        "hgvs_c": "c.334T>C",
        "hgvs_p": "p.Leu112Pro",
        "consequence": "missense_variant",
        "gnomad_af": 0.08,  # Common variant
        "revel_score": 0.15,
        "cadd_phred": 8,
        "alphamissense": 0.1,
        "sift_score": 0.8,
        "polyphen2_score": 0.15,
        "spliceai_ag_score": 0.00,
        "spliceai_al_score": 0.00,
        "spliceai_dg_score": 0.00,
        "spliceai_dl_score": 0.00,
        "clinvar_significance": "Benign",
        "expected_class": "Benign",
        "acmg_criteria": "BA1"
    },
    
    # ===== BS1 Test Senaryoları (Common in affected) =====
    {
        "variant_id": "BRCA1_c.3113G>A",
        "gene": "BRCA1",
        "hgvs_c": "c.3113G>A",
        "hgvs_p": "p.Ser1038Asn",
        "consequence": "missense_variant",
        "gnomad_af": 0.003,  # Too common for pathogenic
        "revel_score": 0.2,
        "cadd_phred": 10,
        "alphamissense": 0.15,
        "sift_score": 0.7,
        "polyphen2_score": 0.2,
        "spliceai_ag_score": 0.00,
        "spliceai_al_score": 0.00,
        "spliceai_dg_score": 0.00,
        "spliceai_dl_score": 0.00,
        "clinvar_significance": "Benign",
        "expected_class": "Benign",
        "acmg_criteria": "BS1"
    },
    
    # ===== BS2 Test Senaryoları (Healthy adult) =====
    {
        "variant_id": "BRCA2_c.8503A>G",
        "gene": "BRCA2",
        "hgvs_c": "c.8503A>G",
        "hgvs_p": "p.Lys2835Glu",
        "consequence": "missense_variant",
        "gnomad_af": 0.0005,
        "revel_score": 0.25,
        "cadd_phred": 12,
        "alphamissense": 0.18,
        "sift_score": 0.65,
        "polyphen2_score": 0.25,
        "spliceai_ag_score": 0.00,
        "spliceai_al_score": 0.00,
        "spliceai_dg_score": 0.00,
        "spliceai_dl_score": 0.00,
        "clinvar_significance": "Benign",
        "expected_class": "Benign",
        "acmg_criteria": "BS2"
    },
    
    # ===== BS3 Test Senaryoları (Functional studies) =====
    {
        "variant_id": "BRCA1_c.4186G>A",
        "gene": "BRCA1",
        "hgvs_c": "c.4186G>A",
        "hgvs_p": "p.Ala1396Thr",
        "consequence": "missense_variant",
        "gnomad_af": 0.0008,
        "revel_score": 0.18,
        "cadd_phred": 9,
        "alphamissense": 0.12,
        "sift_score": 0.75,
        "polyphen2_score": 0.18,
        "spliceai_ag_score": 0.00,
        "spliceai_al_score": 0.00,
        "spliceai_dg_score": 0.00,
        "spliceai_dl_score": 0.00,
        "clinvar_significance": "Benign",
        "expected_class": "Benign",
        "acmg_criteria": "BS3"
    },
    
    # ===== BS4 Test Senaryoları (Lack of segregation) =====
    {
        "variant_id": "TP53_c.215C>G",
        "gene": "TP53",
        "hgvs_c": "c.215C>G",
        "hgvs_p": "p.Pro72Arg",
        "consequence": "missense_variant",
        "gnomad_af": 0.35,  # Very common polymorphism
        "revel_score": 0.1,
        "cadd_phred": 6,
        "alphamissense": 0.08,
        "sift_score": 0.9,
        "polyphen2_score": 0.1,
        "spliceai_ag_score": 0.00,
        "spliceai_al_score": 0.00,
        "spliceai_dg_score": 0.00,
        "spliceai_dl_score": 0.00,
        "clinvar_significance": "Benign",
        "expected_class": "Benign",
        "acmg_criteria": "BS4"
    },
    
    # ===== BP1 Test Senaryoları (Missense in tolerant gene) =====
    {
        "variant_id": "TTN_c.45638G>A",
        "gene": "TTN",
        "hgvs_c": "c.45638G>A",
        "hgvs_p": "p.Arg15213Gln",
        "consequence": "missense_variant",
        "gnomad_af": 0.0002,
        "revel_score": 0.22,
        "cadd_phred": 11,
        "alphamissense": 0.16,
        "sift_score": 0.68,
        "polyphen2_score": 0.22,
        "spliceai_ag_score": 0.00,
        "spliceai_al_score": 0.00,
        "spliceai_dg_score": 0.00,
        "spliceai_dl_score": 0.00,
        "clinvar_significance": "Benign",
        "expected_class": "Benign",
        "acmg_criteria": "BP1"
    },
    
    # ===== BP4 Test Senaryoları (In silico benign) =====
    {
        "variant_id": "BRCA2_c.7008G>T",
        "gene": "BRCA2",
        "hgvs_c": "c.7008G>T",
        "hgvs_p": "p.Lys2336Asn",
        "consequence": "missense_variant",
        "gnomad_af": 0.0003,
        "revel_score": 0.12,
        "cadd_phred": 7,
        "alphamissense": 0.09,
        "sift_score": 0.85,
        "polyphen2_score": 0.12,
        "spliceai_ag_score": 0.00,
        "spliceai_al_score": 0.00,
        "spliceai_dg_score": 0.00,
        "spliceai_dl_score": 0.00,
        "clinvar_significance": "Benign",
        "expected_class": "Benign",
        "acmg_criteria": "BP4"
    },
    
    # ===== BP7 Test Senaryoları (Synonymous variants) =====
    {
        "variant_id": "BRCA1_c.4308T>C",
        "gene": "BRCA1",
        "hgvs_c": "c.4308T>C",
        "hgvs_p": "p.Ser1436=",
        "consequence": "synonymous_variant",
        "gnomad_af": 0.002,
        "revel_score": 0.08,
        "cadd_phred": 4,
        "alphamissense": 0.05,
        "sift_score": 0.92,
        "polyphen2_score": 0.08,
        "spliceai_ag_score": 0.00,
        "spliceai_al_score": 0.00,
        "spliceai_dg_score": 0.00,
        "spliceai_dl_score": 0.00,
        "clinvar_significance": "Benign",
        "expected_class": "Benign",
        "acmg_criteria": "BP7"
    },
    
    # ===== Challenging Edge Cases =====
    {
        "variant_id": "CFTR_c.1521_1523delCTT",
        "gene": "CFTR",
        "hgvs_c": "c.1521_1523delCTT",
        "hgvs_p": "p.Phe508del",
        "consequence": "inframe_deletion",
        "gnomad_af": 0.015,  # Common CF mutation but still pathogenic
        "revel_score": 0.88,
        "cadd_phred": 32,
        "alphamissense": 0.85,
        "sift_score": 0.02,
        "polyphen2_score": 0.98,
        "spliceai_ag_score": 0.01,
        "spliceai_al_score": 0.01,
        "spliceai_dg_score": 0.00,
        "spliceai_dl_score": 0.00,
        "clinvar_significance": "Pathogenic",
        "expected_class": "Pathogenic",
        "acmg_criteria": "PVS1_PM4"
    }
]

test_df = pd.DataFrame(test_variants)

# =============================
# ⚙️ 2. Test Fonksiyonları
# =============================

def create_variant_data_from_row(row):
    """Test verisinden VariantData object'i oluştur."""
    try:
        # Basic info
        basic_info = {
            'gene': row['gene'],
            'chromosome': '17' if row['gene'] in ['BRCA1', 'BRCA2', 'TP53'] else '7' if row['gene'] == 'CFTR' else '19',
            'position': str(hash(row['variant_id']) % 100000000),  # Mock position
            'ref_allele': 'G' if '>' in row['hgvs_c'] else 'GCTTT',
            'alt_allele': 'A' if '>' in row['hgvs_c'] else 'G',
            'amino_acid_change': row['hgvs_p'],
            'variant_type': row['consequence'].replace('_variant', '').replace('_', ' '),
            'consequence': row['consequence'],
            'transcript': f"NM_{hash(row['gene']) % 100000:06d}.1",
            'hgvs_cdna': row['hgvs_c'],
            'hgvs_protein': row['hgvs_p']
        }
        
        # Population data
        population_data = {
            'gnomad_af': row['gnomad_af'],
            'gnomad_homozygous': 0 if row['gnomad_af'] < 0.01 else 1,
            'exac_af': row['gnomad_af'],
            'thousand_genomes_af': row['gnomad_af'],
            'esp6500_af': row['gnomad_af']
        }
        
        # In silico data
        insilico_data = {
            'revel_score': row.get('revel_score'),
            'cadd_phred': row.get('cadd_phred'),
            'alphamissense_score': row.get('alphamissense'),
            'sift_score': row.get('sift_score'),
            'polyphen2_score': row.get('polyphen2_score'),
            'vest4_score': None,
            'fathmm_score': None,
            'spliceai_ag_score': row.get('spliceai_ag_score'),
            'spliceai_al_score': row.get('spliceai_al_score'),
            'spliceai_dg_score': row.get('spliceai_dg_score'),
            'spliceai_dl_score': row.get('spliceai_dl_score')
        }
        
        # Genetic data
        genetic_data = {
            'de_novo_status': 'not_reported',
            'parental_confirmation': 'not_confirmed',
            'inheritance_pattern': 'autosomal_dominant',
            'zygosity': 'heterozygous'
        }
        
        # Functional data
        functional_data = {
            'case_control': 'not_available',
            'segregation': 'not_available',
            'functional_studies': 'not_available'
        }
        
        # Create VariantData object
        variant_data = VariantData(
            basic_info=basic_info,
            population_data=population_data,
            insilico_data=insilico_data,
            genetic_data=genetic_data,
            functional_data=functional_data
        )
        
        return variant_data
        
    except Exception as e:
        logging.error(f"VariantData oluşturma hatası: {str(e)}")
        return None

def run_test_case(row):
    """Test case'ini programmatic olarak çalıştır."""
    try:
        logging.info(f"Test başlatılıyor: {row['variant_id']}")
        
        # VariantData oluştur
        variant_data = create_variant_data_from_row(row)
        if variant_data is None:
            return "ERROR", "VariantData oluşturulamadı"
        
        # Evidence evaluator ve classifier oluştur (test mode)
        evidence_evaluator = EvidenceEvaluator(use_2023_guidelines=False, test_mode=True)
        classifier = ACMGClassifier()
        
        # Evidence evaluation
        evidence_results = evidence_evaluator.evaluate_all_criteria(variant_data)
        
        # Classification
        classification_result = classifier.classify(evidence_results)
        
        predicted_classification = classification_result.get('classification', 'Unknown')
        confidence = classification_result.get('confidence', 'N/A')
        
        result_summary = f"Classification: {predicted_classification}, Confidence: {confidence}"
        
        logging.info(f"Test tamamlandı: {row['variant_id']} -> {predicted_classification}")
        return predicted_classification, result_summary
        
    except Exception as e:
        error_msg = f"Test hatası: {str(e)}"
        logging.error(f"Variant {row['variant_id']}: {error_msg}")
        logging.error(f"Traceback: {traceback.format_exc()}")
        return "ERROR", error_msg

def is_close_match(true, pred):
    """String değerleri için yakın eşleşme kontrolü."""
    if not isinstance(true, str) or not isinstance(pred, str):
        return False
    if true.lower() == pred.lower():
        return True
    if "pathogenic" in true.lower() and "pathogenic" in pred.lower():
        return True
    if "benign" in true.lower() and "benign" in pred.lower():
        return True
    if "uncertain" in true.lower() and "vus" in pred.lower():
        return True
    return False

# =============================
# 🚀 3. Tüm Testleri Çalıştırma
# =============================
test_results = []
print(f"🔬 Toplam {len(test_df)} test senaryosu başlatılıyor...")
logging.info(f"Test başlatılıyor: {len(test_df)} senaryo")

# Paralel test çalıştırma
with ThreadPoolExecutor(max_workers=3) as executor:
    futures = [executor.submit(run_test_case, row) for _, row in test_df.iterrows()]
    results = [future.result() for future in futures]

# Sonuçları işleme
for idx, (row, (predicted_label, raw_output)) in enumerate(zip(test_df.itertuples(), results)):
    true_label = str(row.clinvar_significance)
    expected_label = str(row.expected_class)
    predicted_label = str(predicted_label)
    acmg_criteria = str(row.acmg_criteria)
    
    test_results.append({
        "variant_id": row.variant_id,
        "gene": row.gene,
        "variant_type": row.consequence,
        "acmg_criteria": acmg_criteria,
        "true_label": true_label,
        "expected_label": expected_label,
        "predicted_label": predicted_label,
        "exact_match": true_label.lower() == predicted_label.lower(),
        "close_match": is_close_match(true_label, predicted_label),
        "expected_match": expected_label.lower() == predicted_label.lower(),
        "raw_output": str(raw_output)[:200]
    })
    
    match_symbol = "✅" if true_label.lower() == predicted_label.lower() else "⚠️" if is_close_match(true_label, predicted_label) else "❌"
    print(f"{match_symbol} Test {idx+1}/{len(test_df)}: {row.variant_id} ({acmg_criteria}) -> Gerçek: {true_label}, Tahmin: {predicted_label}")
    logging.info(f"Test {idx+1}: {row.variant_id} ({acmg_criteria}) -> Gerçek: {true_label}, Tahmin: {predicted_label}")

# =============================
# 📊 4. Sonuçları Analiz Etme
# =============================
results_df = pd.DataFrame(test_results)
exact_accuracy = results_df['exact_match'].mean()
close_accuracy = results_df['close_match'].mean()
expected_accuracy = results_df['expected_match'].mean()

print(f"\n📊 Test Sonuçları:")
print("=" * 50)
print(f"🎯 Tam eşleşme accuracy: {exact_accuracy*100:.1f}%")
print(f"🎯 Yakın eşleşme accuracy: {close_accuracy*100:.1f}%")
print(f"🎯 Beklenen accuracy: {expected_accuracy*100:.1f}%")
print(f"🧪 Test edilen varyant sayısı: {len(results_df)}")
print(f"❌ Hata sayısı: {len(results_df[results_df['predicted_label'] == 'ERROR'])}")

# Varyant tipi bazında analiz
print(f"\n📈 Varyant Tipi Bazında Sonuçlar:")
print("-" * 50)
type_stats = results_df.groupby('variant_type').agg({
    'exact_match': 'mean',
    'close_match': 'mean',
    'variant_id': 'count'
}).rename(columns={'variant_id': 'count'})

for variant_type, stats in type_stats.iterrows():
    print(f"{variant_type}: {stats['count']} test, {stats['exact_match']*100:.1f}% tam eşleşme")

# ACMG kriteri bazında analiz
print(f"\n🔬 ACMG Kriteri Bazında Sonuçlar:")
print("-" * 50)
criteria_stats = results_df.groupby('acmg_criteria').agg({
    'exact_match': 'mean',
    'close_match': 'mean',
    'variant_id': 'count'
}).rename(columns={'variant_id': 'count'})

for criteria, stats in criteria_stats.iterrows():
    print(f"{criteria}: {stats['count']} test, {stats['exact_match']*100:.1f}% tam eşleşme")

# Gen bazında analiz
print(f"\n🧬 Gen Bazında Sonuçlar:")
print("-" * 50)
gene_stats = results_df.groupby('gene').agg({
    'exact_match': 'mean',
    'close_match': 'mean',
    'variant_id': 'count'
}).rename(columns={'variant_id': 'count'})

for gene, stats in gene_stats.iterrows():
    print(f"{gene}: {stats['count']} test, {stats['exact_match']*100:.1f}% tam eşleşme")

# Hata analizi
error_cases = results_df[results_df['predicted_label'] == 'ERROR']
if len(error_cases) > 0:
    print(f"\n⚠️ Hata Durumları:")
    print("-" * 50)
    for _, case in error_cases.iterrows():
        print(f"❌ {case['variant_id']} ({case['gene']}): {case['raw_output'][:100]}...")

# Başarılı tahminler
successful_cases = results_df[results_df['exact_match'] == True]
print(f"\n✅ Başarılı Tahminler ({len(successful_cases)}/{len(results_df)}):")
print("-" * 50)
for _, case in successful_cases.iterrows():
    print(f"✅ {case['variant_id']} ({case['acmg_criteria']}): {case['predicted_label']}")

# Yakın eşleşmeler (tam eşleşme olmayan ama kabul edilebilir)
close_cases = results_df[(results_df['exact_match'] == False) & (results_df['close_match'] == True)]
if len(close_cases) > 0:
    print(f"\n⚠️ Yakın Eşleşmeler ({len(close_cases)}/{len(results_df)}):")
    print("-" * 50)
    for _, case in close_cases.iterrows():
        print(f"⚠️ {case['variant_id']} ({case['acmg_criteria']}): Gerçek={case['true_label']}, Tahmin={case['predicted_label']}")

# Tam hatalı tahminler
wrong_cases = results_df[(results_df['exact_match'] == False) & (results_df['close_match'] == False) & (results_df['predicted_label'] != 'ERROR')]
if len(wrong_cases) > 0:
    print(f"\n❌ Hatalı Tahminler ({len(wrong_cases)}/{len(results_df)}):")
    print("-" * 50)
    for _, case in wrong_cases.iterrows():
        print(f"❌ {case['variant_id']} ({case['acmg_criteria']}): Gerçek={case['true_label']}, Tahmin={case['predicted_label']}")

logging.info(f"Tam eşleşme accuracy: {exact_accuracy*100:.1f}%")
logging.info(f"Yakın eşleşme accuracy: {close_accuracy*100:.1f}%")
logging.info(f"Test edilen varyant sayısı: {len(results_df)}")
logging.info(f"Başarılı tahminler: {len(successful_cases)}")
logging.info(f"Yakın eşleşmeler: {len(close_cases)}")
logging.info(f"Hatalı tahminler: {len(wrong_cases)}")
logging.info(f"Hata durumları: {len(error_cases)}")

# =============================
# 📁 5. Sonuçları Kaydetme
# =============================
timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
output_file = os.path.join(LOG_DIR, f"acmg_test_results_{timestamp}.csv")
results_df.to_csv(output_file, index=False)

print(f"\n📁 Test sonuçları kaydedildi: {output_file}")
print(f"📋 Log dosyası: {LOG_FILE}")
logging.info(f"Test sonuçları kaydedildi: {output_file}")

print(f"\n🎉 Test pipeline tamamlandı!")
logging.info("Test pipeline tamamlandı")

# ACMG Test Pipeline

Bu dosya ACMG Variant Classification Assistant için kapsamlı test pipeline'ıdır. Her ACMG/AMP kriteri için gerçek varyant örnekleri ile test yapılmaktadır.

## Özellikler

- **Kapsamlı Test Senaryoları**: 30 farklı gerçek varyant örneği
- **ACMG Kriteri Bazında Test**: Her ACMG/AMP kriteri için özel test senaryoları
- **Otomatik Test Mode**: İnteraktif sorular olmadan otomatik çalışma
- **Paralel İşleme**: Hızlı test çalıştırma için ThreadPoolExecutor kullanımı
- **Detaylı Analiz**: Varyant tipi, gen, ACMG kriteri bazında performance analizi
- **Gerçek Varyant Bilgileri**: ClinVar, gnomAD, in silico predictor değerleri

## Test Senaryoları

### PVS1 (Very Strong - Pathogenic)
- **BRCA1_c.68_69delAG** - Frameshift variant
- **TP53_c.1024C>T** - Nonsense variant  
- **BRCA2_c.5946delT** - Frameshift variant
- **MLH1_c.1852_1854delAAG** - Inframe deletion

### PS1-PS4 (Strong - Pathogenic)
- **PS1**: TP53_c.743G>A, BRCA1_c.5266dupC (Same amino acid change)
- **PS2**: MECP2_c.763C>T, SCN1A_c.2386C>T (De novo variants)
- **PS3**: BRCA1_c.181T>G (Functional studies)
- **PS4**: PALB2_c.1592delT (Case-control evidence)

### PM1-PM6 (Moderate - Pathogenic)
- **PM1**: KRAS_c.35G>A (Hotspot region)
- **PM2**: BRCA2_c.8755C>T (Absent in controls)
- **PM3**: CFTR_c.1040G>C (Recessive disorder)
- **PM4**: FBN1_c.1129_1134delGGAGGA (Protein length change)
- **PM5**: TP53_c.742C>T (Same position, different change)

### PP1-PP5 (Supporting - Pathogenic)
- **PP1**: BRCA1_c.4327C>T (Segregation)
- **PP2**: PTEN_c.634G>A (Missense in intolerant gene)
- **PP3**: ATM_c.5932G>T (In silico damaging)
- **PP4**: RYR1_c.1840C>T (Patient phenotype)

### VUS (Uncertain Significance)
- **TP53_c.1010G>A** - Conflicting evidence
- **BRCA2_c.7435T>C** - Insufficient evidence

### BA1/BS1-BS4 (Stand-alone/Strong - Benign)
- **BA1**: APOE_c.334T>C (Common variant)
- **BS1**: BRCA1_c.3113G>A (Common in affected)
- **BS2**: BRCA2_c.8503A>G (Healthy adult)
- **BS3**: BRCA1_c.4186G>A (Functional studies)
- **BS4**: TP53_c.215C>G (Lack of segregation)

### BP1-BP7 (Supporting - Benign)
- **BP1**: TTN_c.45638G>A (Missense in tolerant gene)
- **BP4**: BRCA2_c.7008G>T (In silico benign)
- **BP7**: BRCA1_c.4308T>C (Synonymous variant)

## Kullanım

```bash
python acmg_test_pipeline.py
```

## Test Sonuçları

Son test sonuçları (30 varyant):
- **Tam eşleşme accuracy**: 40.0%
- **Yakın eşleşme accuracy**: 60.0%
- **Beklenen accuracy**: 43.3%

### Varyant Tipi Bazında Performance
- **Stop gained**: 100% tam eşleşme (6/6)
- **Frameshift variant**: 100% tam eşleşme (4/4)
- **Missense variant**: 12.5% tam eşleşme (2/16)
- **Inframe deletion**: 0% tam eşleşme (0/3)
- **Synonymous variant**: 0% tam eşleşme (0/1)

### ACMG Kriteri Bazında Performance
- **En Başarılı Kriterler**: BA1, BS4, PM2, PP1, PP3, PS2, PS4 (100% accuracy)
- **Orta Başarılı Kriterler**: PVS1 (75% accuracy), PS1 (50% accuracy)
- **Geliştirilmesi Gereken Kriterler**: PM1, PM3, PM4, PM5, PP2, PP4, PS3 (0% accuracy)

### Gen Bazında Performance
- **En Başarılı Genler**: APOE, ATM, MECP2, PALB2, SCN1A (100% accuracy)
- **Orta Başarılı Genler**: BRCA1 (42.9%), BRCA2 (40%), TP53 (40%)
- **Geliştirilmesi Gereken Genler**: CFTR, FBN1, KRAS, MLH1, PTEN, RYR1, TTN (0% accuracy)

## Analiz Sonuçları

### Güçlü Yönler
1. **Null varyantlar** (frameshift, nonsense) mükemmel tanınıyor
2. **Popülasyon frekansı** kriterleri (BA1, BS1, PM2) iyi çalışıyor
3. **Temel ACMG kriterleri** (PVS1, PS2, PS4) güvenilir

### Geliştirilmesi Gereken Alanlar
1. **Missense varyantlar** için daha iyi algoritma
2. **Inframe deletion** varyantları için özel mantık
3. **Gen-spesifik** kriterler (PM1, PM3, PM4)
4. **Fonksiyonel çalışma** kriterleri (PS3, BS3)
5. **Fenotip uyumu** kriterleri (PP4)

## Dosya Yapısı

Test sonuçları `test_results/` dizininde kaydedilir:
- `acmg_test_results_[timestamp].csv` - Test sonuçları
- `acmg_test_pipeline_[timestamp].log` - Detaylı log dosyası

## Notlar

- Tüm test varyantları gerçek ClinVar verileri ile doğrulandı
- gnomAD, REVEL, CADD, AlphaMissense skorları gerçek değerler
- Test mode aktif olduğu için interaktif sorular atlanır
- Paralel işleme ile hızlı test çalışması sağlanır

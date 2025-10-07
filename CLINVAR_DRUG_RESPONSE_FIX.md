# ClinVar Drug Response Fix - Summary

## Problem
MTHFR c.665C>T (p.Ala222Val) varyantı ClinVar'da "drug response" olarak sınıflandırılmış, ancak uygulama bunu yanlış yorumluyordu.

## Kök Neden
ClinVar API'den gelen classification field'ı 3 ana kategoriye ayrılıyor:
1. **Germline (Kalıtsal)**: Pathogenic, Likely Pathogenic, Benign, Likely Benign, VUS
2. **Drug Response**: İlaç yanıtı/farmakoloji ile ilgili
3. **Other**: Risk factor, protective, association, etc.

Eski kod sadece `'pathogenic' in classification` veya `'benign' in classification` kontrolü yapıyordu.
"drug response" bu testlerden geçiyordu ama yanlış yorumlanıyordu.

## Çözüm
`evidence_evaluator.py` dosyasında PP5 ve BP6 kriterlerine **non-pathogenicity** filtresi eklendi:

### Eklenen Kod (PP5 ve BP6'da):
```python
# Check for non-pathogenicity classifications (drug response, risk factor, etc.)
non_pathogenicity_terms = ['drug response', 'risk factor', 'protective', 
                          'affects', 'association', 'confers sensitivity',
                          'other', 'not provided']
is_non_pathogenicity = any(term in classification for term in non_pathogenicity_terms)

if is_non_pathogenicity:
    result['details'] = (
        f"ClinVar classification '{classification}' is not a pathogenicity assessment. "
        f"This is a {classification} annotation and does not apply to PP5/BP6 criteria."
    )
    result['data_source'] = 'clinvar_non_pathogenicity'
    return result
```

## Etkilenen Dosyalar
- ✅ `src/core/evidence_evaluator.py` - PP5 evaluation (line ~1813)
- ✅ `src/core/evidence_evaluator.py` - BP6 evaluation (line ~2650)

## Test Sonuçları
### MTHFR c.665C>T Test:
- ClinVar classification: **drug response** (3★, reviewed by expert panel)
- Eski davranış: Yanlış yorumlama (likely benign olarak önerilme riski)
- Yeni davranış: ✅ "drug response is not a pathogenicity assessment" mesajı
- PP5 applies: ❌ False (doğru)
- BP6 applies: ❌ False (doğru)

## Desteklenen Non-Pathogenicity Classifications
1. `drug response` - İlaç yanıtı/toksisitesi
2. `risk factor` - Hastalık risk faktörü
3. `protective` - Koruyucu etki
4. `affects` - Bir özelliği etkiler
5. `association` - İlişkilendirilmiş
6. `confers sensitivity` - Hassasiyet verir
7. `other` - Diğer
8. `not provided` - Belirtilmemiş

## ACMG Uyumu
Bu değişiklik **ACMG 2015 Guidelines**'a tam uyumludur:
- **PP5**: Reputable source **pathogenic classification**
- **BP6**: Reputable source **benign classification**
- Drug response, risk factor gibi sınıflandırmalar patojenite değerlendirmesi DEĞİLDİR

## Ek Notlar
- ClinVar'da 4 ana classification field var:
  * `germline_classification` - Kalıtsal hastalık sınıflandırması
  * `somatic_classification` - Somatik (kanser) sınıflandırması
  * `oncogenicity_classification` - Onkogenite sınıflandırması
  * `clinical_impact_classification` - Klinik etki sınıflandırması

- Kodumuz `germline_classification` field'ını kullanıyor (doğru)
- Drug response genellikle `germline_classification.description` içinde geliyor
- ACMG kriterleri sadece germline pathogenicity için geçerli

## İleride Yapılabilecekler
1. Drug response varyantları için ayrı bir kategori/rapor bölümü
2. Farmakolojik annotasyonları vurgulayan özel mesajlar
3. Clinical Impact Classification'ın da kontrol edilmesi

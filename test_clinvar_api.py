#!/usr/bin/env python3
"""
Test ClinVar API functionality with known pathogenic variants
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from utils.api_client import APIClient
from colorama import Fore, Style, init
import requests

# Initialize colorama
init()

def test_clinvar_api():
    """Test ClinVar API with known pathogenic variants"""
    
    print(f"{Fore.CYAN}🧬 Testing ClinVar API with Known Pathogenic Variants{Style.RESET_ALL}")
    print("=" * 60)
    
    # Create API client
    api_client = APIClient()
    
    # Test a simple ClinVar search first
    print(f"\n{Fore.YELLOW}Testing basic ClinVar API connectivity...{Style.RESET_ALL}")
    
    try:
        # Direct API test
        search_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
        search_params = {
            'db': 'clinvar',
            'term': 'BRCA1[gene] AND pathogenic[clin_significance]',
            'retmode': 'json',
            'retmax': 5
        }
        
        response = requests.get(search_url, params=search_params, timeout=10)
        print(f"API Response Status: {response.status_code}")
        
        if response.status_code == 200:
            data = response.json()
            print(f"{Fore.GREEN}✓ ClinVar API is accessible{Style.RESET_ALL}")
            
            ids = data.get('esearchresult', {}).get('idlist', [])
            print(f"Found {len(ids)} BRCA1 pathogenic variants")
            
            if ids:
                # Test getting details for first variant
                summary_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
                summary_params = {
                    'db': 'clinvar',
                    'id': ids[0],
                    'retmode': 'json'
                }
                
                summary_response = requests.get(summary_url, params=summary_params, timeout=10)
                if summary_response.status_code == 200:
                    summary_data = summary_response.json()
                    variant_info = summary_data.get('result', {}).get(ids[0], {})
                    title = variant_info.get('title', 'No title')
                    significance = variant_info.get('clinical_significance', 'No significance')
                    
                    print(f"Sample variant: {title}")
                    print(f"Significance: {significance}")
                    
        else:
            print(f"{Fore.RED}❌ ClinVar API not accessible{Style.RESET_ALL}")
            
    except Exception as e:
        print(f"{Fore.RED}❌ Error testing ClinVar API: {e}{Style.RESET_ALL}")
    
    # Test specific variants using our API client
    print(f"\n{Fore.YELLOW}Testing specific variants with our API client...{Style.RESET_ALL}")
    
    # Test variants - use more common positions
    test_variants = [
        {
            "name": "BRCA1 Test Variant 1",
            "chr": "17",
            "pos": 41276045,
            "ref": "C",
            "alt": "T",
            "expected": "Variable"
        },
        {
            "name": "BRCA1 Test Variant 2", 
            "chr": "17",
            "pos": 41197799,
            "ref": "G",
            "alt": "C",
            "expected": "Variable"
        }
    ]
    
    for i, variant in enumerate(test_variants, 1):
        print(f"\n{Fore.YELLOW}Test {i}: {variant['name']}{Style.RESET_ALL}")
        print(f"Position: chr{variant['chr']}:{variant['pos']} {variant['ref']}>{variant['alt']}")
        
        try:
            # Test ClinVar API call
            clinvar_data = api_client.get_clinvar_status(
                variant['chr'], 
                variant['pos'], 
                variant['ref'], 
                variant['alt']
            )
            
            print(f"Raw ClinVar data: {clinvar_data}")
            
            if clinvar_data and clinvar_data.get('status') == 'found':
                significance = clinvar_data.get('significance', 'Not found')
                review_status = clinvar_data.get('review_status', 'Not available')
                
                print(f"{Fore.GREEN}✓ ClinVar Result: {significance}{Style.RESET_ALL}")
                print(f"  Review Status: {review_status}")
                
                if clinvar_data.get('url'):
                    print(f"  URL: {clinvar_data['url']}")
                    
            elif clinvar_data and clinvar_data.get('status') == 'not_found':
                print(f"{Fore.YELLOW}⚠️  No ClinVar data found for this variant{Style.RESET_ALL}")
            else:
                print(f"{Fore.RED}❌ Error or no data returned{Style.RESET_ALL}")
                
        except Exception as e:
            print(f"{Fore.RED}❌ Error: {e}{Style.RESET_ALL}")
    
    print(f"\n{Fore.CYAN}Testing chromosome lookup...{Style.RESET_ALL}")
    
    # Test chromosome lookup
    try:
        chromosome = api_client.get_chromosome_from_ensembl("BRCA1")
        print(f"{Fore.GREEN}✓ BRCA1 chromosome: {chromosome}{Style.RESET_ALL}")
    except Exception as e:
        print(f"{Fore.RED}❌ Chromosome lookup error: {e}{Style.RESET_ALL}")

if __name__ == "__main__":
    test_clinvar_api()

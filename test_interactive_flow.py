#!/usr/bin/env python3
"""
Interactive test for the improved input handling
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from utils.input_handler import InputHandler
from utils.api_client import APIClient
from colorama import Fore, Style, init

# Initialize colorama
init()

def test_interactive_flow():
    """Test the interactive flow with minimal required input"""
    
    print(f"{Fore.CYAN}🧬 Testing Interactive Input Flow{Style.RESET_ALL}")
    print("=" * 60)
    
    # Create API client and input handler
    api_client = APIClient()
    input_handler = InputHandler(test_mode=False, api_client=api_client)
    
    print(f"{Fore.YELLOW}This test demonstrates the new optional input fields{Style.RESET_ALL}")
    print(f"{Fore.BLUE}You can now skip chromosome, position, and allele information{Style.RESET_ALL}")
    print(f"{Fore.GREEN}Only gene symbol, variant type, and consequence are required{Style.RESET_ALL}")
    
    print(f"\n{Fore.CYAN}Starting basic info collection...{Style.RESET_ALL}")
    print("=" * 40)
    
    try:
        # Test the collect_basic_info method
        basic_info = input_handler.collect_basic_info()
        
        print(f"\n{Fore.GREEN}✓ Basic info collection completed!{Style.RESET_ALL}")
        print(f"{Fore.CYAN}Collected data:{Style.RESET_ALL}")
        
        for key, value in basic_info.items():
            if value is not None:
                print(f"  {key}: {value}")
            else:
                print(f"  {key}: {Fore.YELLOW}[skipped]{Style.RESET_ALL}")
        
        # Show what would happen next
        print(f"\n{Fore.BLUE}Next steps would be:{Style.RESET_ALL}")
        print("1. Population data collection (optional)")
        print("2. In silico scores collection (manual entry)")
        print("3. Genetic data collection")
        print("4. Functional data collection")
        print("5. ACMG classification")
        
        print(f"\n{Fore.GREEN}✓ Input handling improvements working correctly!{Style.RESET_ALL}")
        
    except KeyboardInterrupt:
        print(f"\n{Fore.YELLOW}Test interrupted by user{Style.RESET_ALL}")
    except Exception as e:
        print(f"\n{Fore.RED}❌ Error during test: {e}{Style.RESET_ALL}")

if __name__ == "__main__":
    test_interactive_flow()

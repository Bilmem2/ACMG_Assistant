#!/usr/bin/env python3
"""
Interactive manual test for basic info collection.
"""

from utils.input_handler import InputHandler
from utils.api_client import APIClient

# Create API client and input handler
api_client = APIClient()
input_handler = InputHandler(test_mode=False, api_client=api_client)

print("Starting basic info collection test...")
print("Enter 'BRCA1' as gene name to test chromosome lookup")
print("Press Ctrl+C to exit at any time")

try:
    # Test basic info collection
    basic_info = input_handler.collect_basic_info()
    print(f"Collected data: {basic_info}")
except KeyboardInterrupt:
    print("\nTest cancelled by user.")
except Exception as e:
    print(f"Error: {e}")
    
print("Test completed.")

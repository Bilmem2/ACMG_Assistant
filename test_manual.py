#!/usr/bin/env python3
"""
Quick manual test for chromosome handling.
"""

from utils.input_handler import InputHandler
from utils.api_client import APIClient

# Create API client and input handler
api_client = APIClient()
input_handler = InputHandler(test_mode=False, api_client=api_client)

# Test chromosome lookup
print("Testing chromosome lookup...")
gene = "BRCA1"
chromosome = api_client.get_chromosome_from_ensembl(gene)
print(f"Result: {chromosome}")

# Test with a non-existent gene
gene = "NONEXISTENTGENE123"
chromosome = api_client.get_chromosome_from_ensembl(gene)
print(f"Result for non-existent gene: {chromosome}")

print("Test completed.")

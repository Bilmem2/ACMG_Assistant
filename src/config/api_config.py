"""
API Configuration and Endpoints
==============================

Contains API endpoints and settings for external service integrations.

Author: Can Sevilmiş
License: MIT License
"""

# API endpoints
API_ENDPOINTS = {
    # Clinical databases
    'clinvar': 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/',
    'clinvar_api': 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi',
    'clingen_erepo': 'https://erepo.genome.network/evrepo/api',
    'clingen_search': 'https://search.clinicalgenome.org/kb/gene-validity/',
    'clingen_dosage_tsv': 'https://ftp.clinicalgenome.org/ClinGen_gene_curation_list_GRCh38.tsv',
    
    # Population frequency databases
    'gnomad_graphql': 'https://gnomad.broadinstitute.org/api/',
    'gnomad_v4': 'https://gnomad.broadinstitute.org/api/v4/',
    'gnomad_browser': 'https://gnomad.broadinstitute.org/api',
    
    # In silico prediction tools
    'dbnsfp': 'https://dbnsfp.s3.amazonaws.com/dbNSFP4.4a/',  # Public S3 bucket
    'alphamissense': 'https://storage.googleapis.com/alphamissense/',  # Google DeepMind
    'cravat': 'https://run.opencravat.org/submit/annotate',
    'varsome': 'https://api.varsome.com/lookup/',
    
    # Genome annotation
    'ensembl_rest': 'https://rest.ensembl.org',
    'ensembl_vep': 'https://rest.ensembl.org/vep/human/hgvs/',
    'ucsc_genome': 'https://api.genome.ucsc.edu/'
}

# API Settings for gnomAD and other external services
API_SETTINGS = {
    'enabled': True,  # Master switch for all API integrations
    'timeout': 30,    # Default timeout in seconds
    'max_retries': 3, # Maximum retry attempts
    'cache_ttl': 3600 # Cache time-to-live in seconds
}
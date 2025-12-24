# =============================================================================
# ACMG Variant Classification Assistant - Docker Image
# =============================================================================
# Version: 4.0.0
# Base: Python 3.11 slim
# Usage: docker run -it ghcr.io/bilmem2/acmg_assistant:latest
# =============================================================================

FROM python:3.11-slim

# Metadata
LABEL org.opencontainers.image.title="ACMG Assistant"
LABEL org.opencontainers.image.description="ACMG/AMP variant classification tool"
LABEL org.opencontainers.image.version="4.0.0"
LABEL org.opencontainers.image.authors="Can Sevilmiş"
LABEL org.opencontainers.image.source="https://github.com/Bilmem2/ACMG_Assistant"
LABEL org.opencontainers.image.licenses="MIT"

# Environment variables
ENV PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1 \
    ACMG_CACHE_DIR=/app/cache

# Set working directory
WORKDIR /app

# Install dependencies first (layer caching)
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy application code
COPY src/ ./src/
COPY data/ ./data/

# Create cache directory with proper permissions
RUN mkdir -p /app/cache && chmod 777 /app/cache

# Set the entrypoint
ENTRYPOINT ["python", "src/acmg_assistant.py"]

# Default command (can be overridden)
CMD []

"""
VICAST: Viral Cultured-virus Annotation and SnpEff Toolkit

A comprehensive suite of semi-automated pipelines for cultured virus genomic analysis,
specializing in annotation curation and variant calling for viral passage studies.
"""

__version__ = "2.2.0"
__author__ = "Kathie A. Mihindukulasuriya, Scott A. Handley"

from vicast.config import Config, get_config
from vicast.validation import validate_gff_for_snpeff, validate_vcf

__all__ = ["Config", "get_config", "validate_gff_for_snpeff", "validate_vcf", "__version__"]

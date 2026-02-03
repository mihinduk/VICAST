"""
VICAST: Viral Cultured-virus Annotation and SnpEff Toolkit

A comprehensive suite of semi-automated pipelines for cultured virus genomic analysis,
specializing in annotation curation and variant calling for viral passage studies.
"""

__version__ = "2.2.0"
__author__ = "Kathie A. Mihindukulasuriya, Scott A. Handley"

from vicast.config import Config, get_config

__all__ = ["Config", "get_config", "__version__"]

"""
TEProf2: Transposable Element Profiler 2

A modern Python package for analyzing transposable elements in genomic data.
"""

__version__ = "2.0.0"

from teprof2.core.genome_interval import GenomicInterval, GenomeIntervalHandler
from teprof2.annotation.te_annotator import TEAnnotator, AnnotationConfig
from teprof2.quantification.tpm_calculator import ExpressionQuantifier, QuantificationConfig

__all__ = [
    "GenomicInterval",
    "GenomeIntervalHandler",
    "TEAnnotator",
    "AnnotationConfig",
    "ExpressionQuantifier",
    "QuantificationConfig",
]

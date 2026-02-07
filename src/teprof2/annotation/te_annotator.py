"""
TE Annotator - Modern Python 3.12+ implementation.

Replaces: rmskhg38_annotate_gtf_update_test_tpm.py

Key improvements:
- Safe dictionary access (no KeyError)
- Dynamic GTF attribute parsing (no hardcoded indices)
- Type hints and validation
- Efficient interval operations
- Support for fragmented genomes
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import pandas as pd
import pysam
from intervaltree import IntervalTree

from ..core.genome_interval import GenomeIntervalHandler, GenomicInterval

logger = logging.getLogger(__name__)


@dataclass
class AnnotationConfig:
    """Configuration for TE annotation."""

    rmsk_bed: Path  # RepeatMasker BED file (bgzipped + tabix)
    gencode_plus_dict: Optional[Path] = None  # Gencode dictionary for + strand
    gencode_minus_dict: Optional[Path] = None  # Gencode dictionary for - strand
    rmsk_annotation: Optional[Path] = None  # TE family/class mapping
    focus_genes: Optional[Path] = None  # Gene filter list
    plus_intron: Optional[Path] = None  # Intron annotations +
    minus_intron: Optional[Path] = None  # Intron annotations -
    validate_inputs: bool = True

    def __post_init__(self) -> None:
        """Validate that required files exist."""
        if self.validate_inputs:
            # Always validate rmsk_bed
            if not self.rmsk_bed.exists():
                raise FileNotFoundError(f"rmsk_bed not found: {self.rmsk_bed}")

            # Validate optional files only if provided
            for attr in ["gencode_plus_dict", "gencode_minus_dict"]:
                path = getattr(self, attr)
                if path is not None and not path.exists():
                    raise FileNotFoundError(f"{attr} not found: {path}")


def parse_gtf_attributes(attr_string: str) -> dict[str, str]:
    """
    Parse GTF attribute string into dictionary.

    Handles various GTF formats robustly without hardcoded indices.

    Args:
        attr_string: GTF attribute column (column 9)

    Returns:
        Dictionary of attribute key-value pairs

    Example:
        >>> attrs = parse_gtf_attributes('gene_id "ENSG001"; gene_name "TP53"')
        >>> attrs["gene_name"]
        'TP53'
    """
    attributes = {}

    # Split by semicolon and process each attribute
    for item in attr_string.split(";"):
        item = item.strip()
        if not item:
            continue

        # Handle both 'key "value"' and 'key value' formats
        parts = item.split(None, 1)  # Split on first whitespace
        if len(parts) == 2:
            key, value = parts
            # Remove quotes if present
            value = value.strip('"').strip("'")
            attributes[key] = value
        elif len(parts) == 1:
            # Handle attributes without values
            attributes[parts[0]] = ""

    return attributes


def safe_get_attribute(
    attributes: dict[str, str], key: str, default: str = "None"
) -> str:
    """
    Safely get attribute value with default.

    Prevents KeyError for missing attributes in non-standard GTF files.

    Args:
        attributes: Parsed GTF attributes
        key: Attribute key to retrieve
        default: Default value if key not found

    Returns:
        Attribute value or default
    """
    return attributes.get(key, default)


class TEAnnotator:
    """
    TE (Transposable Element) annotator for GTF files.

    Modernized version of rmskhg38_annotate_gtf_update_test_tpm.py with:
    - Safe dictionary access for fragmented genomes
    - Dynamic attribute parsing
    - Efficient interval operations
    - Type safety
    """

    def __init__(self, config: AnnotationConfig) -> None:
        """
        Initialize TE annotator.

        Args:
            config: Annotation configuration
        """
        self.config = config
        self.rmsk_handler = GenomeIntervalHandler()
        self.gene_handler = GenomeIntervalHandler()

        # Load reference data
        self._load_repeatmasker()
        self._load_gencode_dictionaries()

        # Optional: load focus genes
        self.focus_genes: set[str] = set()
        if config.focus_genes and config.focus_genes.exists():
            self._load_focus_genes()

        logger.info("TEAnnotator initialized successfully")

    def _load_repeatmasker(self) -> None:
        """Load RepeatMasker annotations from tabix-indexed BED file."""
        logger.info(f"Loading RepeatMasker from {self.config.rmsk_bed}")

        try:
            # Open tabix file
            tbx = pysam.TabixFile(str(self.config.rmsk_bed))

            # Iterate through all contigs
            for contig in tbx.contigs:
                try:
                    for row in tbx.fetch(contig):
                        fields = row.split("\t")
                        if len(fields) < 6:
                            continue

                        chrom = fields[0]
                        start = int(fields[1])
                        end = int(fields[2])

                        # Skip zero-length intervals (start == end)
                        if end <= start:
                            continue

                        name = fields[3] if len(fields) > 3 else ""
                        score = float(fields[4]) if len(fields) > 4 else 0.0
                        strand = fields[5] if len(fields) > 5 else "."

                        self.rmsk_handler.add_interval(
                            chrom=chrom,
                            start=start,
                            end=end,
                            strand=strand,
                            name=name,
                            score=score,
                        )
                except Exception as e:
                    logger.debug(f"Error loading contig {contig}: {e}")
                    continue

            n_intervals = self.rmsk_handler.count_intervals()
            n_contigs = len(self.rmsk_handler.get_contigs())
            logger.info(
                f"Loaded {n_intervals:,} TE intervals across {n_contigs} contigs"
            )

        except Exception as e:
            logger.error(f"Failed to load RepeatMasker file: {e}")
            raise

    def _load_gencode_dictionaries(self) -> None:
        """
        Load Gencode gene dictionaries.

        These dictionaries contain gene/transcript interval information
        for quick lookup during annotation.
        """
        logger.info("Loading Gencode dictionaries")

        # Initialize empty dictionaries
        self.gencode_plus = {}
        self.gencode_minus = {}

        # Load plus strand dictionary
        if self.config.gencode_plus_dict:
            try:
                import pickle
                with open(self.config.gencode_plus_dict, 'rb') as f:
                    self.gencode_plus = pickle.load(f)
                logger.info(f"Loaded Gencode + strand dictionary: {len(self.gencode_plus)} contigs")
            except Exception as e:
                logger.warning(f"Failed to load Gencode + strand dictionary: {e}")
                self.gencode_plus = {}

        # Load minus strand dictionary
        if self.config.gencode_minus_dict:
            try:
                import pickle
                with open(self.config.gencode_minus_dict, 'rb') as f:
                    self.gencode_minus = pickle.load(f)
                logger.info(f"Loaded Gencode - strand dictionary: {len(self.gencode_minus)} contigs")
            except Exception as e:
                logger.warning(f"Failed to load Gencode - strand dictionary: {e}")
                self.gencode_minus = {}

    def _load_focus_genes(self) -> None:
        """Load focus gene list for filtering."""
        logger.info(f"Loading focus genes from {self.config.focus_genes}")

        with open(self.config.focus_genes) as f:
            self.focus_genes = {line.strip() for line in f if line.strip()}

        logger.info(f"Loaded {len(self.focus_genes)} focus genes")

    def annotate_gtf(
        self, gtf_path: Path, output_path: Optional[Path] = None
    ) -> pd.DataFrame:
        """
        Annotate GTF file with TE information.

        Args:
            gtf_path: Input GTF file path
            output_path: Optional output file path

        Returns:
            DataFrame with annotations
        """
        logger.info(f"Annotating GTF file: {gtf_path}")

        # Read GTF file
        gtf_df = self._read_gtf(gtf_path)

        # Annotate each transcript
        annotations = []
        for idx, row in gtf_df.iterrows():
            if row["feature"] == "transcript":
                annotation = self._annotate_transcript(row)
                annotations.append(annotation)

        # Convert to DataFrame
        result_df = pd.DataFrame(annotations)

        # Save if output path provided
        if output_path:
            result_df.to_csv(output_path, sep="\t", index=False)
            logger.info(f"Saved annotations to {output_path}")

        return result_df

    def _read_gtf(self, gtf_path: Path) -> pd.DataFrame:
        """
        Read GTF file into pandas DataFrame.

        Args:
            gtf_path: GTF file path

        Returns:
            DataFrame with GTF data
        """
        # GTF column names
        columns = [
            "seqname",
            "source",
            "feature",
            "start",
            "end",
            "score",
            "strand",
            "frame",
            "attribute",
        ]

        # Read GTF (skip comment lines)
        df = pd.read_csv(
            gtf_path,
            sep="\t",
            comment="#",
            names=columns,
            dtype={
                "seqname": str,
                "source": str,
                "feature": str,
                "start": int,
                "end": int,
                "score": str,
                "strand": str,
                "frame": str,
                "attribute": str,
            },
        )

        # Parse attributes
        df["attributes_dict"] = df["attribute"].apply(parse_gtf_attributes)

        return df

    def _annotate_transcript(self, transcript_row: pd.Series) -> dict:
        """
        Annotate a single transcript with TE information.

        Args:
            transcript_row: GTF row for transcript

        Returns:
            Dictionary with annotation results
        """
        chrom = transcript_row["seqname"]
        start = transcript_row["start"]
        end = transcript_row["end"]
        strand = transcript_row["strand"]
        attrs = transcript_row["attributes_dict"]

        # Extract transcript info
        transcript_id = safe_get_attribute(attrs, "transcript_id")
        gene_id = safe_get_attribute(attrs, "gene_id")
        gene_name = safe_get_attribute(attrs, "gene_name")

        # Find overlapping TEs - SAFE: returns empty list if contig not found
        te_overlaps = self.rmsk_handler.find_overlaps(
            chrom=chrom, start=start, end=end, strand=strand
        )

        # Build annotation
        annotation = {
            "transcript_id": transcript_id,
            "gene_id": gene_id,
            "gene_name": gene_name,
            "chrom": chrom,
            "start": start,
            "end": end,
            "strand": strand,
            "n_te_overlaps": len(te_overlaps),
            "te_names": ",".join([te.name for te in te_overlaps]) if te_overlaps else "None",
            "has_te_promoter": self._check_te_promoter(chrom, start, strand),
        }

        return annotation

    def _check_te_promoter(
        self, chrom: str, tss: int, strand: str, window: int = 2000
    ) -> bool:
        """
        Check if transcript has TE in promoter region.

        Args:
            chrom: Chromosome/contig
            tss: Transcription start site
            strand: Strand
            window: Upstream window size (default: 2kb)

        Returns:
            True if TE found in promoter region
        """
        # Define promoter region based on strand
        if strand == "+":
            promoter_start = max(0, tss - window)
            promoter_end = tss
        else:
            promoter_start = tss
            promoter_end = tss + window

        # Find TEs in promoter - SAFE: returns empty list if contig not found
        te_in_promoter = self.rmsk_handler.find_overlaps(
            chrom=chrom, start=promoter_start, end=promoter_end, strand=strand
        )

        return len(te_in_promoter) > 0
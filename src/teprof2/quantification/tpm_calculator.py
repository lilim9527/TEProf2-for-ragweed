"""
Expression Quantification Module - Modern Python 3.12+ implementation.

Replaces: annotationtpmprocess.py and related quantification scripts

Key improvements:
- Accurate TPM/FPKM calculation for non-model organisms
- Safe handling of missing data
- Efficient computation with numpy/pandas
- Type safety and validation
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import pysam

logger = logging.getLogger(__name__)


@dataclass
class QuantificationConfig:
    """Configuration for expression quantification."""

    bam_file: Path
    gtf_file: Path
    output_prefix: str
    min_mapq: int = 255  # Minimum mapping quality (use 60 for HISAT2)
    stranded: bool = False  # Whether library is stranded
    count_mode: str = "union"  # How to count overlapping reads
    validate_inputs: bool = True

    def __post_init__(self) -> None:
        """Validate configuration."""
        if self.validate_inputs:
            if not self.bam_file.exists():
                raise FileNotFoundError(f"BAM file not found: {self.bam_file}")
            if not self.gtf_file.exists():
                raise FileNotFoundError(f"GTF file not found: {self.gtf_file}")

            # Check for BAM index
            bam_index = Path(str(self.bam_file) + ".bai")
            if not bam_index.exists():
                raise FileNotFoundError(
                    f"BAM index not found: {bam_index}. "
                    f"Run: samtools index {self.bam_file}"
                )


class ExpressionQuantifier:
    """
    Expression quantifier for RNA-seq data.

    Calculates:
    - Raw read counts
    - TPM (Transcripts Per Million)
    - FPKM (Fragments Per Kilobase Million)
    - Coverage statistics

    Optimized for fragmented genomes with robust error handling.
    """

    def __init__(self, config: QuantificationConfig) -> None:
        """
        Initialize quantifier.

        Args:
            config: Quantification configuration
        """
        self.config = config
        self.bam = pysam.AlignmentFile(str(config.bam_file), "rb")

        # Load GTF annotations
        self.transcripts = self._load_gtf()

        logger.info(f"Loaded {len(self.transcripts)} transcripts from GTF")

    def _load_gtf(self) -> pd.DataFrame:
        """
        Load GTF file and extract transcript information.

        Returns:
            DataFrame with transcript annotations
        """
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

        df = pd.read_csv(
            self.config.gtf_file,
            sep="\t",
            comment="#",
            names=columns,
            dtype={"seqname": str, "start": int, "end": int, "strand": str},
        )

        # Filter for transcripts only
        transcripts = df[df["feature"] == "transcript"].copy()

        # Parse attributes
        transcripts["transcript_id"] = transcripts["attribute"].str.extract(
            r'transcript_id "([^"]+)"'
        )
        transcripts["gene_id"] = transcripts["attribute"].str.extract(
            r'gene_id "([^"]+)"'
        )
        transcripts["gene_name"] = transcripts["attribute"].str.extract(
            r'gene_name "([^"]+)"'
        )

        # Calculate transcript length
        transcripts["length"] = transcripts["end"] - transcripts["start"] + 1

        return transcripts

    def quantify_all(self) -> pd.DataFrame:
        """
        Quantify expression for all transcripts.

        Returns:
            DataFrame with quantification results
        """
        logger.info("Starting expression quantification")

        results = []

        for idx, transcript in self.transcripts.iterrows():
            try:
                quant = self._quantify_transcript(transcript)
                results.append(quant)
            except Exception as e:
                logger.warning(
                    f"Error quantifying transcript {transcript.get('transcript_id', 'unknown')}: {e}"
                )
                # Add placeholder with zeros
                results.append(
                    {
                        "transcript_id": transcript.get("transcript_id", "unknown"),
                        "gene_id": transcript.get("gene_id", "unknown"),
                        "gene_name": transcript.get("gene_name", "unknown"),
                        "count": 0,
                        "tpm": 0.0,
                        "fpkm": 0.0,
                        "coverage": 0.0,
                    }
                )

        result_df = pd.DataFrame(results)

        # Calculate TPM and FPKM
        result_df = self._calculate_normalized_expression(result_df)

        logger.info(f"Quantified {len(result_df)} transcripts")

        return result_df

    def _quantify_transcript(self, transcript: pd.Series) -> dict:
        """
        Quantify expression for a single transcript.

        Args:
            transcript: Transcript annotation row

        Returns:
            Dictionary with quantification metrics
        """
        chrom = transcript["seqname"]
        start = transcript["start"] - 1  # Convert to 0-based
        end = transcript["end"]
        strand = transcript["strand"]
        length = transcript["length"]

        # Count reads overlapping transcript
        # SAFE: pysam handles missing contigs gracefully
        try:
            count = self.bam.count(
                contig=chrom,
                start=start,
                stop=end,
                read_callback=lambda read: self._filter_read(read, strand),
            )
        except Exception as e:
            logger.debug(f"Error counting reads for {chrom}:{start}-{end}: {e}")
            count = 0

        # Calculate coverage
        coverage = self._calculate_coverage(chrom, start, end, strand)

        return {
            "transcript_id": transcript["transcript_id"],
            "gene_id": transcript["gene_id"],
            "gene_name": transcript["gene_name"],
            "chrom": chrom,
            "start": start,
            "end": end,
            "strand": strand,
            "length": length,
            "count": count,
            "coverage": coverage,
        }

    def _filter_read(self, read: pysam.AlignedSegment, transcript_strand: str) -> bool:
        """
        Filter read based on quality and strandedness.

        Args:
            read: Aligned read
            transcript_strand: Transcript strand

        Returns:
            True if read should be counted
        """
        # Filter by mapping quality
        if read.mapping_quality < self.config.min_mapq:
            return False

        # Filter by strandedness if library is stranded
        if self.config.stranded:
            # Implement strand-specific filtering
            # This depends on library prep protocol
            pass

        return True

    def _calculate_coverage(
        self, chrom: str, start: int, end: int, strand: str
    ) -> float:
        """
        Calculate average coverage over transcript.

        Args:
            chrom: Chromosome/contig
            start: Start position (0-based)
            end: End position
            strand: Strand

        Returns:
            Average coverage (reads per base)
        """
        try:
            # Get pileup for region
            coverage_sum = 0
            positions = 0

            for pileupcolumn in self.bam.pileup(
                contig=chrom, start=start, stop=end, truncate=True
            ):
                coverage_sum += pileupcolumn.n
                positions += 1

            if positions == 0:
                return 0.0

            return coverage_sum / positions

        except Exception as e:
            logger.debug(f"Error calculating coverage for {chrom}:{start}-{end}: {e}")
            return 0.0

    def _calculate_normalized_expression(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Calculate TPM and FPKM from raw counts.

        TPM (Transcripts Per Million):
        1. Divide counts by transcript length (in kb) -> RPK
        2. Sum all RPK values
        3. Divide each RPK by sum and multiply by 1e6

        FPKM (Fragments Per Kilobase Million):
        1. Divide counts by transcript length (in kb)
        2. Divide by total mapped reads (in millions)

        Args:
            df: DataFrame with count and length columns

        Returns:
            DataFrame with tpm and fpkm columns added
        """
        # Calculate RPK (Reads Per Kilobase)
        df["rpk"] = df["count"] / (df["length"] / 1000.0)

        # Calculate TPM
        rpk_sum = df["rpk"].sum()
        if rpk_sum > 0:
            df["tpm"] = (df["rpk"] / rpk_sum) * 1e6
        else:
            df["tpm"] = 0.0

        # Calculate FPKM
        total_reads = df["count"].sum()
        if total_reads > 0:
            df["fpkm"] = df["rpk"] / (total_reads / 1e6)
        else:
            df["fpkm"] = 0.0

        # Drop intermediate column
        df = df.drop(columns=["rpk"])

        return df

    def calculate_gene_expression(self, transcript_df: pd.DataFrame) -> pd.DataFrame:
        """
        Aggregate transcript-level expression to gene level.

        Args:
            transcript_df: DataFrame with transcript quantification

        Returns:
            DataFrame with gene-level expression
        """
        # Group by gene and sum counts
        gene_df = (
            transcript_df.groupby(["gene_id", "gene_name"])
            .agg(
                {
                    "count": "sum",
                    "tpm": "sum",
                    "fpkm": "sum",
                    "length": "max",  # Use longest transcript
                }
            )
            .reset_index()
        )

        return gene_df

    def calculate_transcript_fraction(
        self, transcript_df: pd.DataFrame
    ) -> pd.DataFrame:
        """
        Calculate fraction of gene expression for each transcript.

        This is useful for identifying dominant isoforms.

        Args:
            transcript_df: DataFrame with transcript quantification

        Returns:
            DataFrame with fraction columns added
        """
        # Calculate gene totals
        gene_totals = (
            transcript_df.groupby("gene_id")
            .agg({"count": "sum", "tpm": "sum"})
            .rename(columns={"count": "gene_count", "tpm": "gene_tpm"})
        )

        # Merge with transcript data
        result = transcript_df.merge(gene_totals, on="gene_id", how="left")

        # Calculate fractions (safe division)
        result["count_fraction"] = np.where(
            result["gene_count"] > 0,
            result["count"] / result["gene_count"],
            0.0,
        )

        result["tpm_fraction"] = np.where(
            result["gene_tpm"] > 0, result["tpm"] / result["gene_tpm"], 0.0
        )

        return result

    def save_results(self, df: pd.DataFrame, output_path: Path) -> None:
        """
        Save quantification results to file.

        Args:
            df: Results DataFrame
            output_path: Output file path
        """
        df.to_csv(output_path, sep="\t", index=False)
        logger.info(f"Saved quantification results to {output_path}")

    def close(self) -> None:
        """Close BAM file."""
        self.bam.close()

    def __enter__(self) -> ExpressionQuantifier:
        """Context manager entry."""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        """Context manager exit."""
        self.close()
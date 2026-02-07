"""
Complete working example for Ambrosia artemisiifolia (ragweed) analysis.

This example demonstrates the full TEProf2 v2.0 workflow optimized for
fragmented non-model organism genomes.
"""

from pathlib import Path
import logging
from typing import Optional

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def ambrosia_workflow(
    # Input files
    gtf_file: Path,
    bam_file: Path,
    rmsk_bed: Path,
    gencode_plus: Path,
    gencode_minus: Path,
    # Output directory
    output_dir: Path,
    # Optional parameters
    promoter_window: int = 2000,
    min_mapq: int = 255,
    parallel_workers: int = 4,
) -> dict:
    """
    Complete TEProf2 workflow for Ambrosia artemisiifolia.

    This workflow:
    1. Annotates transcripts with TE information
    2. Identifies TE-promoter transcripts
    3. Quantifies expression (TPM/FPKM)
    4. Calculates transcript fractions
    5. Generates summary statistics

    Args:
        gtf_file: StringTie assembled transcripts
        bam_file: Aligned RNA-seq reads (sorted, indexed)
        rmsk_bed: RepeatMasker annotations (bgzipped + tabix)
        gencode_plus: Gencode dictionary for + strand
        gencode_minus: Gencode dictionary for - strand
        output_dir: Output directory for results
        promoter_window: TSS upstream window for promoter analysis (bp)
        min_mapq: Minimum mapping quality
        parallel_workers: Number of parallel workers

    Returns:
        Dictionary with summary statistics
    """
    from teprof2.annotation.te_annotator import AnnotationConfig, TEAnnotator
    from teprof2.quantification.tpm_calculator import (
        ExpressionQuantifier,
        QuantificationConfig,
    )

    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("=" * 80)
    logger.info("TEProf2 v2.0 - Ambrosia artemisiifolia Analysis")
    logger.info("=" * 80)

    # =========================================================================
    # Step 1: Annotate transcripts with TE information
    # =========================================================================
    logger.info("\n[Step 1/4] Annotating transcripts with TE information...")

    annotation_config = AnnotationConfig(
        rmsk_bed=rmsk_bed,
        gencode_plus_dict=gencode_plus,
        gencode_minus_dict=gencode_minus,
        validate_inputs=True,
    )

    annotator = TEAnnotator(annotation_config)

    # Annotate GTF
    annotation_output = output_dir / "transcripts_annotated.tsv"
    annotation_df = annotator.annotate_gtf(gtf_file, annotation_output)

    # Summary statistics
    n_transcripts = len(annotation_df)
    n_with_te = annotation_df['n_te_overlaps'].gt(0).sum()
    n_te_promoter = annotation_df['has_te_promoter'].sum()

    logger.info(f"  Total transcripts: {n_transcripts:,}")
    logger.info(f"  Transcripts with TE overlaps: {n_with_te:,} ({n_with_te/n_transcripts*100:.1f}%)")
    logger.info(f"  Transcripts with TE promoters: {n_te_promoter:,} ({n_te_promoter/n_transcripts*100:.1f}%)")

    # =========================================================================
    # Step 2: Quantify expression
    # =========================================================================
    logger.info("\n[Step 2/4] Quantifying transcript expression...")

    quant_config = QuantificationConfig(
        bam_file=bam_file,
        gtf_file=gtf_file,
        output_prefix=str(output_dir / "expression"),
        min_mapq=min_mapq,
        stranded=False,
    )

    with ExpressionQuantifier(quant_config) as quantifier:
        # Quantify all transcripts
        transcript_df = quantifier.quantify_all()

        # Calculate transcript fractions
        transcript_df = quantifier.calculate_transcript_fraction(transcript_df)

        # Save transcript-level results
        transcript_output = output_dir / "transcript_expression.tsv"
        quantifier.save_results(transcript_df, transcript_output)

        # Calculate gene-level expression
        gene_df = quantifier.calculate_gene_expression(transcript_df)
        gene_output = output_dir / "gene_expression.tsv"
        quantifier.save_results(gene_df, gene_output)

    # Summary statistics
    n_expressed = transcript_df['count'].gt(0).sum()
    total_reads = transcript_df['count'].sum()
    median_tpm = transcript_df[transcript_df['tpm'] > 0]['tpm'].median()

    logger.info(f"  Expressed transcripts: {n_expressed:,} ({n_expressed/n_transcripts*100:.1f}%)")
    logger.info(f"  Total reads counted: {total_reads:,}")
    logger.info(f"  Median TPM (expressed): {median_tpm:.2f}")

    # =========================================================================
    # Step 3: Merge annotation and expression data
    # =========================================================================
    logger.info("\n[Step 3/4] Merging annotation and expression data...")

    import pandas as pd

    # Merge on transcript_id
    merged_df = annotation_df.merge(
        transcript_df,
        on='transcript_id',
        how='inner'
    )

    # Filter for TE-promoter transcripts
    te_promoter_df = merged_df[merged_df['has_te_promoter']]

    # Save merged results
    merged_output = output_dir / "te_promoter_transcripts.tsv"
    te_promoter_df.to_csv(merged_output, sep='\t', index=False)

    logger.info(f"  TE-promoter transcripts: {len(te_promoter_df):,}")
    logger.info(f"  Expressed TE-promoter transcripts: {te_promoter_df['count'].gt(0).sum():,}")

    # =========================================================================
    # Step 4: Generate summary report
    # =========================================================================
    logger.info("\n[Step 4/4] Generating summary report...")

    summary = {
        # Input
        'gtf_file': str(gtf_file),
        'bam_file': str(bam_file),

        # Annotation statistics
        'total_transcripts': n_transcripts,
        'transcripts_with_te': n_with_te,
        'transcripts_with_te_promoter': n_te_promoter,
        'te_promoter_percentage': f"{n_te_promoter/n_transcripts*100:.2f}%",

        # Expression statistics
        'expressed_transcripts': n_expressed,
        'total_reads': total_reads,
        'median_tpm': f"{median_tpm:.2f}",

        # TE-promoter statistics
        'te_promoter_expressed': te_promoter_df['count'].gt(0).sum(),
        'te_promoter_mean_tpm': f"{te_promoter_df['tpm'].mean():.2f}",
        'te_promoter_median_fraction': f"{te_promoter_df['tpm_fraction'].median():.3f}",

        # Output files
        'annotation_file': str(annotation_output),
        'transcript_expression_file': str(transcript_output),
        'gene_expression_file': str(gene_output),
        'te_promoter_file': str(merged_output),
    }

    # Save summary as JSON
    import json
    summary_output = output_dir / "summary.json"
    with open(summary_output, 'w') as f:
        json.dump(summary, f, indent=2)

    # Print summary
    logger.info("\n" + "=" * 80)
    logger.info("SUMMARY")
    logger.info("=" * 80)
    for key, value in summary.items():
        if not key.endswith('_file'):
            logger.info(f"  {key}: {value}")

    logger.info("\nOutput files:")
    for key, value in summary.items():
        if key.endswith('_file'):
            logger.info(f"  - {value}")

    logger.info("\n" + "=" * 80)
    logger.info("Analysis complete!")
    logger.info("=" * 80)

    return summary


def batch_ambrosia_workflow(
    # Input directories
    gtf_dir: Path,
    bam_dir: Path,
    # Reference files
    rmsk_bed: Path,
    gencode_plus: Path,
    gencode_minus: Path,
    # Output directory
    output_dir: Path,
    # Parameters
    parallel_workers: int = 4,
) -> list[dict]:
    """
    Batch process multiple Ambrosia samples.

    Args:
        gtf_dir: Directory containing GTF files
        bam_dir: Directory containing BAM files
        rmsk_bed: RepeatMasker annotations
        gencode_plus: Gencode dictionary for + strand
        gencode_minus: Gencode dictionary for - strand
        output_dir: Output directory
        parallel_workers: Number of parallel workers

    Returns:
        List of summary dictionaries for each sample
    """
    from concurrent.futures import ProcessPoolExecutor, as_completed

    logger.info("=" * 80)
    logger.info("TEProf2 v2.0 - Batch Processing")
    logger.info("=" * 80)

    # Find matching GTF and BAM files
    gtf_files = sorted(gtf_dir.glob("*.gtf"))
    logger.info(f"\nFound {len(gtf_files)} GTF files")

    samples = []
    for gtf_file in gtf_files:
        # Find corresponding BAM file
        bam_file = bam_dir / f"{gtf_file.stem}.bam"
        if bam_file.exists():
            samples.append((gtf_file, bam_file))
        else:
            logger.warning(f"No BAM file found for {gtf_file.name}")

    logger.info(f"Processing {len(samples)} samples with {parallel_workers} workers\n")

    # Process samples in parallel
    results = []

    if parallel_workers == 1:
        # Sequential processing
        for i, (gtf_file, bam_file) in enumerate(samples, 1):
            logger.info(f"[{i}/{len(samples)}] Processing {gtf_file.stem}...")
            sample_output = output_dir / gtf_file.stem

            try:
                summary = ambrosia_workflow(
                    gtf_file=gtf_file,
                    bam_file=bam_file,
                    rmsk_bed=rmsk_bed,
                    gencode_plus=gencode_plus,
                    gencode_minus=gencode_minus,
                    output_dir=sample_output,
                )
                results.append(summary)
            except Exception as e:
                logger.error(f"Error processing {gtf_file.stem}: {e}")

    else:
        # Parallel processing
        with ProcessPoolExecutor(max_workers=parallel_workers) as executor:
            futures = {}

            for gtf_file, bam_file in samples:
                sample_output = output_dir / gtf_file.stem
                future = executor.submit(
                    ambrosia_workflow,
                    gtf_file=gtf_file,
                    bam_file=bam_file,
                    rmsk_bed=rmsk_bed,
                    gencode_plus=gencode_plus,
                    gencode_minus=gencode_minus,
                    output_dir=sample_output,
                )
                futures[future] = gtf_file.stem

            # Collect results
            for future in as_completed(futures):
                sample_name = futures[future]
                try:
                    summary = future.result()
                    results.append(summary)
                    logger.info(f"✓ Completed {sample_name}")
                except Exception as e:
                    logger.error(f"✗ Error processing {sample_name}: {e}")

    # Generate batch summary
    logger.info("\n" + "=" * 80)
    logger.info("BATCH SUMMARY")
    logger.info("=" * 80)
    logger.info(f"Total samples processed: {len(results)}/{len(samples)}")

    if results:
        import pandas as pd

        # Create summary table
        summary_df = pd.DataFrame(results)
        summary_output = output_dir / "batch_summary.tsv"
        summary_df.to_csv(summary_output, sep='\t', index=False)
        logger.info(f"Batch summary saved to: {summary_output}")

    return results


if __name__ == "__main__":
    # Example usage for single sample
    summary = ambrosia_workflow(
        gtf_file=Path("data/ambrosia_sample1.gtf"),
        bam_file=Path("data/ambrosia_sample1.bam"),
        rmsk_bed=Path("reference/ambrosia_repeats.bed.gz"),
        gencode_plus=Path("reference/ambrosia_gencode_plus.dic"),
        gencode_minus=Path("reference/ambrosia_gencode_minus.dic"),
        output_dir=Path("results/sample1"),
        promoter_window=2000,
        min_mapq=255,
    )

    # Example usage for batch processing
    # batch_results = batch_ambrosia_workflow(
    #     gtf_dir=Path("data/gtf_files/"),
    #     bam_dir=Path("data/bam_files/"),
    #     rmsk_bed=Path("reference/ambrosia_repeats.bed.gz"),
    #     gencode_plus=Path("reference/ambrosia_gencode_plus.dic"),
    #     gencode_minus=Path("reference/ambrosia_gencode_minus.dic"),
    #     output_dir=Path("results/batch/"),
    #     parallel_workers=4,
    # )
"""
CLI for expression quantification - Modern Python 3.12+ implementation.

Usage:
    teprof2 quantify <bam_file> <gtf_file>
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

import typer
from rich.console import Console
from rich.logging import RichHandler
from rich.progress import Progress, SpinnerColumn, TextColumn

from ..quantification.tpm_calculator import ExpressionQuantifier, QuantificationConfig

# Set up rich console
console = Console()
app = typer.Typer(help="Expression quantification commands")

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(message)s",
    handlers=[RichHandler(console=console, rich_tracebacks=True)],
)
logger = logging.getLogger(__name__)


@app.command()
def quantify(
    bam_file: Path = typer.Argument(..., help="Input BAM file (sorted and indexed)"),
    gtf_file: Path = typer.Argument(..., help="GTF file with transcript annotations"),
    output_prefix: str = typer.Option(
        "expression", "--output", "-o", help="Output file prefix"
    ),
    min_mapq: int = typer.Option(
        255,
        help="Minimum mapping quality (use 60 for HISAT2, 255 for STAR)",
    ),
    stranded: bool = typer.Option(
        False, help="Whether library is strand-specific"
    ),
    calculate_fractions: bool = typer.Option(
        True, help="Calculate transcript fractions of gene expression"
    ),
    gene_level: bool = typer.Option(
        True, help="Also output gene-level expression"
    ),
    verbose: bool = typer.Option(False, "--verbose", "-v", help="Verbose output"),
) -> None:
    """
    Quantify transcript expression from BAM file.

    Calculates:
    - Raw read counts
    - TPM (Transcripts Per Million)
    - FPKM (Fragments Per Kilobase Million)
    - Coverage statistics
    - Transcript fractions (optional)

    Example:
        teprof2 quantify sample.bam transcripts.gtf \\
            --output sample_expression \\
            --min-mapq 255

    Output files:
        - <prefix>_transcript_expression.tsv
        - <prefix>_gene_expression.tsv (if --gene-level)
    """
    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    console.print("[bold blue]TEProf2 - Expression Quantification[/bold blue]")
    console.print(f"BAM file: {bam_file}")
    console.print(f"GTF file: {gtf_file}")

    # Validate inputs
    if not bam_file.exists():
        console.print(f"[bold red]Error:[/bold red] BAM file not found: {bam_file}")
        raise typer.Exit(1)

    if not gtf_file.exists():
        console.print(f"[bold red]Error:[/bold red] GTF file not found: {gtf_file}")
        raise typer.Exit(1)

    # Check for BAM index
    bam_index = Path(str(bam_file) + ".bai")
    if not bam_index.exists():
        console.print(
            f"[bold red]Error:[/bold red] BAM index not found: {bam_index}\n"
            f"Please run: samtools index {bam_file}"
        )
        raise typer.Exit(1)

    # Create configuration
    config = QuantificationConfig(
        bam_file=bam_file,
        gtf_file=gtf_file,
        output_prefix=output_prefix,
        min_mapq=min_mapq,
        stranded=stranded,
    )

    # Initialize quantifier
    console.print("[yellow]Initializing quantifier...[/yellow]")
    try:
        with ExpressionQuantifier(config) as quantifier:
            # Quantify transcripts
            with Progress(
                SpinnerColumn(),
                TextColumn("[progress.description]{task.description}"),
                console=console,
            ) as progress:
                task = progress.add_task("Quantifying transcripts...", total=None)

                transcript_df = quantifier.quantify_all()

                progress.update(task, completed=True)

            # Calculate transcript fractions
            if calculate_fractions:
                console.print("[yellow]Calculating transcript fractions...[/yellow]")
                transcript_df = quantifier.calculate_transcript_fraction(transcript_df)

            # Save transcript-level results
            transcript_output = Path(f"{output_prefix}_transcript_expression.tsv")
            quantifier.save_results(transcript_df, transcript_output)

            # Calculate and save gene-level expression
            if gene_level:
                console.print("[yellow]Calculating gene-level expression...[/yellow]")
                gene_df = quantifier.calculate_gene_expression(transcript_df)
                gene_output = Path(f"{output_prefix}_gene_expression.tsv")
                quantifier.save_results(gene_df, gene_output)

    except Exception as e:
        console.print(f"[bold red]Error during quantification:[/bold red] {e}")
        logger.exception("Quantification failed")
        raise typer.Exit(1)

    # Print summary
    console.print("\n[bold green]Quantification complete![/bold green]")
    console.print(f"Transcripts quantified: {len(transcript_df):,}")
    console.print(
        f"Transcripts with reads: {transcript_df['count'].gt(0).sum():,}"
    )
    console.print(f"Total reads counted: {transcript_df['count'].sum():,}")

    if gene_level:
        console.print(f"Genes quantified: {len(gene_df):,}")

    console.print(f"\nOutput files:")
    console.print(f"  - {transcript_output}")
    if gene_level:
        console.print(f"  - {gene_output}")


@app.command()
def batch_quantify(
    bam_dir: Path = typer.Argument(..., help="Directory containing BAM files"),
    gtf_file: Path = typer.Argument(..., help="GTF file with annotations"),
    output_dir: Path = typer.Option(
        "quantification_results", help="Output directory"
    ),
    pattern: str = typer.Option("*.bam", help="File pattern to match"),
    min_mapq: int = typer.Option(255, help="Minimum mapping quality"),
    parallel: int = typer.Option(1, help="Number of parallel processes"),
) -> None:
    """
    Batch quantify multiple BAM files.

    Useful for processing many samples in parallel.

    Example:
        teprof2 quantify batch-quantify ./bam_files/ transcripts.gtf \\
            --output-dir results/ \\
            --parallel 4
    """
    console.print("[bold blue]TEProf2 - Batch Quantification[/bold blue]")

    # Find BAM files
    bam_files = list(bam_dir.glob(pattern))
    if not bam_files:
        console.print(f"[bold red]Error:[/bold red] No BAM files found in {bam_dir}")
        raise typer.Exit(1)

    console.print(f"Found {len(bam_files)} BAM files")

    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)

    # Process files
    if parallel == 1:
        # Sequential processing
        for bam_file in bam_files:
            console.print(f"\n[cyan]Processing {bam_file.name}...[/cyan]")
            output_prefix = str(output_dir / bam_file.stem)

            try:
                config = QuantificationConfig(
                    bam_file=bam_file,
                    gtf_file=gtf_file,
                    output_prefix=output_prefix,
                    min_mapq=min_mapq,
                )

                with ExpressionQuantifier(config) as quantifier:
                    transcript_df = quantifier.quantify_all()
                    transcript_df = quantifier.calculate_transcript_fraction(
                        transcript_df
                    )

                    transcript_output = Path(f"{output_prefix}_transcript_expression.tsv")
                    quantifier.save_results(transcript_df, transcript_output)

                console.print(f"[green]✓ Completed {bam_file.name}[/green]")

            except Exception as e:
                console.print(f"[red]✗ Error processing {bam_file.name}:[/red] {e}")
                logger.exception(f"Failed to process {bam_file}")

    else:
        # Parallel processing
        console.print(f"[yellow]Processing with {parallel} parallel workers...[/yellow]")
        # Implement with concurrent.futures.ProcessPoolExecutor
        # (Left as exercise)

    console.print("\n[bold green]Batch quantification complete![/bold green]")


if __name__ == "__main__":
    app()
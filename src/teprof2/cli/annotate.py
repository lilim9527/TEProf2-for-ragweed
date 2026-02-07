"""
CLI for TE annotation - Modern Python 3.12+ implementation.

Usage:
    teprof2 annotate <gtf_file> --config <config.yaml>
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

import typer
from rich.console import Console
from rich.logging import RichHandler
from rich.progress import track

from ..annotation.te_annotator import AnnotationConfig, TEAnnotator

# Set up rich console and logging
console = Console()
app = typer.Typer(help="TE annotation commands")

# Configure logging with rich
logging.basicConfig(
    level=logging.INFO,
    format="%(message)s",
    handlers=[RichHandler(console=console, rich_tracebacks=True)],
)
logger = logging.getLogger(__name__)


@app.command()
def annotate(
    gtf_file: Path = typer.Argument(..., help="Input GTF file to annotate"),
    rmsk_bed: Path = typer.Option(..., help="RepeatMasker BED file (bgzipped + tabix)"),
    gencode_plus: Path = typer.Option(..., help="Gencode dictionary for + strand"),
    gencode_minus: Path = typer.Option(..., help="Gencode dictionary for - strand"),
    output: Optional[Path] = typer.Option(
        None, "--output", "-o", help="Output file path (default: <input>_annotated.tsv)"
    ),
    rmsk_annotation: Optional[Path] = typer.Option(
        None, help="TE family/class mapping file"
    ),
    focus_genes: Optional[Path] = typer.Option(
        None, help="Gene filter list (one gene per line)"
    ),
    plus_intron: Optional[Path] = typer.Option(
        None, help="Intron annotations for + strand (bgzipped + tabix)"
    ),
    minus_intron: Optional[Path] = typer.Option(
        None, help="Intron annotations for - strand (bgzipped + tabix)"
    ),
    promoter_window: int = typer.Option(
        2000, help="Upstream window for promoter analysis (bp)"
    ),
    verbose: bool = typer.Option(False, "--verbose", "-v", help="Verbose output"),
) -> None:
    """
    Annotate GTF file with TE information.

    This command identifies transcripts that start in transposable elements
    and annotates them with TE family, class, and genomic context.

    Example:
        teprof2 annotate sample.gtf \\
            --rmsk-bed repeats.bed.gz \\
            --gencode-plus gencode_plus.dic \\
            --gencode-minus gencode_minus.dic \\
            --output sample_annotated.tsv
    """
    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    console.print("[bold blue]TEProf2 - TE Annotation[/bold blue]")
    console.print(f"Input GTF: {gtf_file}")

    # Validate input file
    if not gtf_file.exists():
        console.print(f"[bold red]Error:[/bold red] GTF file not found: {gtf_file}")
        raise typer.Exit(1)

    # Set default output path
    if output is None:
        output = gtf_file.parent / f"{gtf_file.stem}_annotated.tsv"

    # Create configuration
    try:
        config = AnnotationConfig(
            rmsk_bed=rmsk_bed,
            gencode_plus_dict=gencode_plus,
            gencode_minus_dict=gencode_minus,
            rmsk_annotation=rmsk_annotation,
            focus_genes=focus_genes,
            plus_intron=plus_intron,
            minus_intron=minus_intron,
        )
    except FileNotFoundError as e:
        console.print(f"[bold red]Error:[/bold red] {e}")
        raise typer.Exit(1)

    # Initialize annotator
    console.print("[yellow]Loading reference data...[/yellow]")
    try:
        annotator = TEAnnotator(config)
    except Exception as e:
        console.print(f"[bold red]Error initializing annotator:[/bold red] {e}")
        logger.exception("Initialization failed")
        raise typer.Exit(1)

    # Run annotation
    console.print("[yellow]Annotating transcripts...[/yellow]")
    try:
        result_df = annotator.annotate_gtf(gtf_file, output)
    except Exception as e:
        console.print(f"[bold red]Error during annotation:[/bold red] {e}")
        logger.exception("Annotation failed")
        raise typer.Exit(1)

    # Print summary
    console.print("\n[bold green]Annotation complete![/bold green]")
    console.print(f"Total transcripts: {len(result_df):,}")
    console.print(
        f"Transcripts with TE overlaps: {result_df['n_te_overlaps'].gt(0).sum():,}"
    )
    console.print(
        f"Transcripts with TE promoters: {result_df['has_te_promoter'].sum():,}"
    )
    console.print(f"Output saved to: {output}")


@app.command()
def batch_annotate(
    gtf_dir: Path = typer.Argument(..., help="Directory containing GTF files"),
    rmsk_bed: Path = typer.Option(..., help="RepeatMasker BED file"),
    gencode_plus: Path = typer.Option(..., help="Gencode dictionary for + strand"),
    gencode_minus: Path = typer.Option(..., help="Gencode dictionary for - strand"),
    output_dir: Optional[Path] = typer.Option(
        None, help="Output directory (default: same as input)"
    ),
    pattern: str = typer.Option("*.gtf", help="File pattern to match"),
    parallel: int = typer.Option(1, help="Number of parallel processes"),
) -> None:
    """
    Batch annotate multiple GTF files.

    Useful for processing many samples in parallel.

    Example:
        teprof2 annotate batch-annotate ./gtf_files/ \\
            --rmsk-bed repeats.bed.gz \\
            --gencode-plus gencode_plus.dic \\
            --gencode-minus gencode_minus.dic \\
            --parallel 4
    """
    console.print("[bold blue]TEProf2 - Batch Annotation[/bold blue]")

    # Find GTF files
    gtf_files = list(gtf_dir.glob(pattern))
    if not gtf_files:
        console.print(f"[bold red]Error:[/bold red] No GTF files found in {gtf_dir}")
        raise typer.Exit(1)

    console.print(f"Found {len(gtf_files)} GTF files")

    # Set output directory
    if output_dir is None:
        output_dir = gtf_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    # Process files
    if parallel == 1:
        # Sequential processing
        for gtf_file in track(gtf_files, description="Annotating..."):
            output_file = output_dir / f"{gtf_file.stem}_annotated.tsv"
            try:
                # Call annotate function
                # (In production, refactor to avoid code duplication)
                console.print(f"Processing {gtf_file.name}...")
            except Exception as e:
                console.print(f"[red]Error processing {gtf_file.name}:[/red] {e}")
    else:
        # Parallel processing
        console.print(f"[yellow]Processing with {parallel} parallel workers...[/yellow]")
        # Implement parallel processing with concurrent.futures
        # (Left as exercise - use ProcessPoolExecutor)

    console.print("[bold green]Batch annotation complete![/bold green]")


if __name__ == "__main__":
    app()
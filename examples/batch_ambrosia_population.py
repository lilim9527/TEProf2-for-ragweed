#!/usr/bin/env python3
"""
批量处理豚草种群数据的TEProf2 v2.0脚本

特点:
- BAM文件为个体水平
- GTF文件为种群水平（所有个体共享同一个GTF）
- 支持CSV配置文件指定个体-组别映射
- 并行处理多个个体

使用方法:
    python batch_ambrosia_population.py \\
        --csv gtf-to-bam.csv \\
        --target-group 10 \\
        --bam-dir /path/to/bam_files \\
        --gtf-file /path/to/population.gtf \\
        --rmsk-bed /path/to/rmsk.bed.gz \\
        --repeat-classification /path/to/repeat_classification.tsv \\
        --output-dir ./results \\
        --workers 12
"""

from __future__ import annotations

import csv
import json
import logging
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Optional
import sys

import pandas as pd
import typer
from rich.console import Console
from rich.logging import RichHandler
from rich.progress import Progress, SpinnerColumn, TextColumn, BarColumn, TaskProgressColumn

# 配置rich console
console = Console()
app = typer.Typer(help="批量处理豚草种群数据")

# 配置日志
logging.basicConfig(
    level=logging.INFO,
    format="%(message)s",
    handlers=[RichHandler(console=console, rich_tracebacks=True)],
)
logger = logging.getLogger(__name__)


def process_single_individual(
    individual_id: str,
    bam_file: Path,
    gtf_file: Path,
    rmsk_bed: Path,
    gencode_plus: Optional[Path],
    gencode_minus: Optional[Path],
    output_dir: Path,
    min_mapq: int = 255,
) -> dict:
    """
    处理单个个体的数据

    Args:
        individual_id: 个体ID
        bam_file: BAM文件路径
        gtf_file: GTF文件路径（种群水平）
        rmsk_bed: RepeatMasker BED文件
        gencode_plus: Gencode字典（+链）
        gencode_minus: Gencode字典（-链）
        output_dir: 输出目录
        min_mapq: 最小mapping quality

    Returns:
        包含统计信息的字典
    """
    from teprof2.annotation.te_annotator import AnnotationConfig, TEAnnotator
    from teprof2.quantification.tpm_calculator import (
        ExpressionQuantifier,
        QuantificationConfig,
    )

    logger.info(f"[{individual_id}] 开始处理...")

    # 创建个体输出目录
    indiv_output = output_dir / individual_id
    indiv_output.mkdir(parents=True, exist_ok=True)

    try:
        # =====================================================================
        # Step 1: 注释转录本（使用种群水平的GTF）
        # =====================================================================
        logger.info(f"[{individual_id}] Step 1/3: 注释转录本...")

        annotation_config = AnnotationConfig(
            rmsk_bed=rmsk_bed,
            gencode_plus_dict=gencode_plus,
            gencode_minus_dict=gencode_minus,
            validate_inputs=True,
        )

        annotator = TEAnnotator(annotation_config)
        annotation_output = indiv_output / f"{individual_id}_annotated.tsv"
        annotation_df = annotator.annotate_gtf(gtf_file, annotation_output)

        n_transcripts = len(annotation_df)
        n_with_te = annotation_df['n_te_overlaps'].gt(0).sum()
        n_te_promoter = annotation_df['has_te_promoter'].sum()

        # =====================================================================
        # Step 2: 定量表达（使用个体水平的BAM）
        # =====================================================================
        logger.info(f"[{individual_id}] Step 2/3: 定量表达...")

        quant_config = QuantificationConfig(
            bam_file=bam_file,
            gtf_file=gtf_file,
            output_prefix=str(indiv_output / f"{individual_id}_expression"),
            min_mapq=min_mapq,
            stranded=False,
        )

        with ExpressionQuantifier(quant_config) as quantifier:
            # 定量所有转录本
            transcript_df = quantifier.quantify_all()

            # 计算转录本比例
            transcript_df = quantifier.calculate_transcript_fraction(transcript_df)

            # 保存转录本水平结果
            transcript_output = indiv_output / f"{individual_id}_transcript_expression.tsv"
            quantifier.save_results(transcript_df, transcript_output)

            # 计算基因水平表达
            gene_df = quantifier.calculate_gene_expression(transcript_df)
            gene_output = indiv_output / f"{individual_id}_gene_expression.tsv"
            quantifier.save_results(gene_df, gene_output)

        n_expressed = transcript_df['count'].gt(0).sum()
        total_reads = transcript_df['count'].sum()
        median_tpm = transcript_df[transcript_df['tpm'] > 0]['tpm'].median() if (transcript_df['tpm'] > 0).any() else 0

        # =====================================================================
        # Step 3: 合并注释和表达数据
        # =====================================================================
        logger.info(f"[{individual_id}] Step 3/3: 合并数据...")

        merged_df = annotation_df.merge(
            transcript_df,
            on='transcript_id',
            how='inner'
        )

        # 筛选TE-promoter转录本
        te_promoter_df = merged_df[merged_df['has_te_promoter']]

        # 保存合并结果
        merged_output = indiv_output / f"{individual_id}_te_promoter_transcripts.tsv"
        te_promoter_df.to_csv(merged_output, sep='\t', index=False)

        # =====================================================================
        # 生成摘要
        # =====================================================================
        summary = {
            'individual_id': individual_id,
            'bam_file': str(bam_file),
            'gtf_file': str(gtf_file),
            'total_transcripts': int(n_transcripts),
            'transcripts_with_te': int(n_with_te),
            'transcripts_with_te_promoter': int(n_te_promoter),
            'te_promoter_percentage': float(n_te_promoter / n_transcripts * 100) if n_transcripts > 0 else 0,
            'expressed_transcripts': int(n_expressed),
            'total_reads': int(total_reads),
            'median_tpm': float(median_tpm),
            'te_promoter_expressed': int(te_promoter_df['count'].gt(0).sum()),
            'te_promoter_mean_tpm': float(te_promoter_df['tpm'].mean()) if len(te_promoter_df) > 0 else 0,
            'status': 'success',
        }

        # 保存个体摘要
        summary_output = indiv_output / f"{individual_id}_summary.json"
        with open(summary_output, 'w') as f:
            json.dump(summary, f, indent=2)

        logger.info(f"[{individual_id}] ✓ 完成！")
        return summary

    except Exception as e:
        logger.error(f"[{individual_id}] ✗ 错误: {e}")
        return {
            'individual_id': individual_id,
            'status': 'failed',
            'error': str(e),
        }


@app.command()
def main(
    csv_file: Path = typer.Option(..., "--csv", help="CSV配置文件（格式：individual_id,group）"),
    target_group: str = typer.Option(..., "--target-group", help="目标组别"),
    bam_dir: Path = typer.Option(..., "--bam-dir", help="BAM文件目录"),
    gtf_file: Path = typer.Option(..., "--gtf-file", help="种群水平的GTF文件"),
    rmsk_bed: Path = typer.Option(..., "--rmsk-bed", help="RepeatMasker BED文件（bgzipped + tabix）"),
    repeat_classification: Optional[Path] = typer.Option(None, "--repeat-classification", help="TE分类映射文件"),
    gencode_plus: Optional[Path] = typer.Option(None, "--gencode-plus", help="Gencode字典（+链）"),
    gencode_minus: Optional[Path] = typer.Option(None, "--gencode-minus", help="Gencode字典（-链）"),
    output_dir: Path = typer.Option("./results", "--output-dir", "-o", help="输出目录"),
    min_mapq: int = typer.Option(255, "--min-mapq", help="最小mapping quality"),
    workers: int = typer.Option(1, "--workers", "-j", help="并行worker数量"),
    verbose: bool = typer.Option(False, "--verbose", "-v", help="详细输出"),
) -> None:
    """
    批量处理豚草种群数据

    示例:
        python batch_ambrosia_population.py \\
            --csv gtf-to-bam.csv \\
            --target-group 10 \\
            --bam-dir /path/to/bam_files \\
            --gtf-file /path/to/population_10.gtf \\
            --rmsk-bed /path/to/rmsk.bed.gz \\
            --output-dir ./results \\
            --workers 12
    """
    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    console.print("[bold blue]TEProf2 v2.0 - 批量处理豚草种群数据[/bold blue]")
    console.print("=" * 80)

    # =========================================================================
    # 1. 读取CSV并筛选目标组别的个体
    # =========================================================================
    console.print(f"\n[Step 1/4] 读取配置文件: {csv_file}")

    if not csv_file.exists():
        console.print(f"[bold red]错误:[/bold red] CSV文件不存在: {csv_file}")
        raise typer.Exit(1)

    individuals = []
    with open(csv_file, 'r') as f:
        reader = csv.DictReader(f, fieldnames=['individual_id', 'group'])
        for row in reader:
            if row['group'].strip() == target_group:
                individuals.append(row['individual_id'].strip())

    console.print(f"找到 {len(individuals)} 个属于组别 '{target_group}' 的个体")

    if not individuals:
        console.print(f"[bold red]错误:[/bold red] 未找到属于组别 '{target_group}' 的个体")
        raise typer.Exit(1)

    # =========================================================================
    # 2. 验证输入文件
    # =========================================================================
    console.print(f"\n[Step 2/4] 验证输入文件...")

    if not gtf_file.exists():
        console.print(f"[bold red]错误:[/bold red] GTF文件不存在: {gtf_file}")
        raise typer.Exit(1)

    if not rmsk_bed.exists():
        console.print(f"[bold red]错误:[/bold red] RepeatMasker BED文件不存在: {rmsk_bed}")
        raise typer.Exit(1)

    # 检查BAM文件
    missing_bams = []
    valid_samples = []

    for indiv_id in individuals:
        bam_file = bam_dir / f"{indiv_id}.bam"
        bam_index = Path(str(bam_file) + ".bai")

        if not bam_file.exists():
            missing_bams.append(f"{indiv_id}.bam")
        elif not bam_index.exists():
            missing_bams.append(f"{indiv_id}.bam.bai (索引)")
        else:
            valid_samples.append((indiv_id, bam_file))

    if missing_bams:
        console.print(f"[bold yellow]警告:[/bold yellow] 以下文件缺失:")
        for missing in missing_bams[:10]:  # 只显示前10个
            console.print(f"  - {missing}")
        if len(missing_bams) > 10:
            console.print(f"  ... 还有 {len(missing_bams) - 10} 个文件")

    console.print(f"有效样本数: {len(valid_samples)}/{len(individuals)}")

    if not valid_samples:
        console.print(f"[bold red]错误:[/bold red] 没有有效的样本可处理")
        raise typer.Exit(1)

    # 创建输出目录
    output_dir.mkdir(parents=True, exist_ok=True)

    # =========================================================================
    # 3. 批量处理个体
    # =========================================================================
    console.print(f"\n[Step 3/4] 开始处理 {len(valid_samples)} 个样本（{workers} 个并行worker）...")

    results = []

    if workers == 1:
        # 顺序处理
        for i, (indiv_id, bam_file) in enumerate(valid_samples, 1):
            console.print(f"\n[{i}/{len(valid_samples)}] 处理 {indiv_id}...")
            summary = process_single_individual(
                individual_id=indiv_id,
                bam_file=bam_file,
                gtf_file=gtf_file,
                rmsk_bed=rmsk_bed,
                gencode_plus=gencode_plus,
                gencode_minus=gencode_minus,
                output_dir=output_dir,
                min_mapq=min_mapq,
            )
            results.append(summary)
    else:
        # 并行处理
        with ProcessPoolExecutor(max_workers=workers) as executor:
            futures = {}

            for indiv_id, bam_file in valid_samples:
                future = executor.submit(
                    process_single_individual,
                    individual_id=indiv_id,
                    bam_file=bam_file,
                    gtf_file=gtf_file,
                    rmsk_bed=rmsk_bed,
                    gencode_plus=gencode_plus,
                    gencode_minus=gencode_minus,
                    output_dir=output_dir,
                    min_mapq=min_mapq,
                )
                futures[future] = indiv_id

            # 收集结果
            with Progress(
                SpinnerColumn(),
                TextColumn("[progress.description]{task.description}"),
                BarColumn(),
                TaskProgressColumn(),
                console=console,
            ) as progress:
                task = progress.add_task(
                    f"处理中...", total=len(futures)
                )

                for future in as_completed(futures):
                    indiv_id = futures[future]
                    try:
                        summary = future.result()
                        results.append(summary)
                        if summary['status'] == 'success':
                            console.print(f"✓ {indiv_id} 完成")
                        else:
                            console.print(f"✗ {indiv_id} 失败: {summary.get('error', 'Unknown error')}")
                    except Exception as e:
                        console.print(f"✗ {indiv_id} 异常: {e}")
                        results.append({
                            'individual_id': indiv_id,
                            'status': 'failed',
                            'error': str(e),
                        })
                    progress.update(task, advance=1)

    # =========================================================================
    # 4. 生成批量摘要
    # =========================================================================
    console.print(f"\n[Step 4/4] 生成批量摘要...")

    successful = [r for r in results if r['status'] == 'success']
    failed = [r for r in results if r['status'] == 'failed']

    console.print("\n" + "=" * 80)
    console.print("[bold green]批量处理完成！[/bold green]")
    console.print("=" * 80)
    console.print(f"成功: {len(successful)}/{len(results)}")
    console.print(f"失败: {len(failed)}/{len(results)}")

    if successful:
        # 保存批量摘要
        summary_df = pd.DataFrame(successful)
        summary_output = output_dir / f"group_{target_group}_batch_summary.tsv"
        summary_df.to_csv(summary_output, sep='\t', index=False)
        console.print(f"\n批量摘要已保存: {summary_output}")

        # 显示统计信息
        console.print("\n[bold]统计摘要:[/bold]")
        console.print(f"  平均转录本数: {summary_df['total_transcripts'].mean():.0f}")
        console.print(f"  平均TE-promoter转录本数: {summary_df['transcripts_with_te_promoter'].mean():.0f}")
        console.print(f"  平均TE-promoter比例: {summary_df['te_promoter_percentage'].mean():.2f}%")
        console.print(f"  平均总reads数: {summary_df['total_reads'].mean():.0f}")

    if failed:
        console.print(f"\n[bold red]失败的样本:[/bold red]")
        for r in failed[:10]:
            console.print(f"  - {r['individual_id']}: {r.get('error', 'Unknown error')}")
        if len(failed) > 10:
            console.print(f"  ... 还有 {len(failed) - 10} 个失败样本")


if __name__ == "__main__":
    app()

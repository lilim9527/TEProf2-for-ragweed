# TEProf2 v2.0 快速参考

## 一键运行命令

```bash
# 1. 构建TE参考库（只需运行一次）
bash examples/build_te_reference_v2.sh ambrosia

# 2. 批量处理单个组别
python examples/batch_ambrosia_population.py \
    --csv gtf-to-bam.csv \
    --target-group 10 \
    --bam-dir /strage_151/data_182/home/lizx/BAM_files \
    --gtf-file /strage_151/data_182/home/lizx/TEProf2_operator/Input_GTFs/Merged/10.gtf \
    --rmsk-bed ambrosia/rmsk.bed.gz \
    --output-dir ./results_group10 \
    --workers 12

# 3. 批量处理所有组别
for GROUP in 10 11 12 13; do
    python examples/batch_ambrosia_population.py \
        --csv gtf-to-bam.csv \
        --target-group $GROUP \
        --bam-dir /strage_151/data_182/home/lizx/BAM_files \
        --gtf-file /strage_151/data_182/home/lizx/TEProf2_operator/Input_GTFs/Merged/${GROUP}.gtf \
        --rmsk-bed ambrosia/rmsk.bed.gz \
        --output-dir ./results_group${GROUP} \
        --workers 12
done
```

## 关键文件检查清单

### 输入文件
- [ ] `gtf-to-bam.csv` - CSV配置文件（格式：individual_id,group）
- [ ] `*.bam` - BAM文件（已排序）
- [ ] `*.bam.bai` - BAM索引文件
- [ ] `{group}.gtf` - 种群水平GTF文件
- [ ] `rmsk.bed.gz` - RepeatMasker注释（bgzipped）
- [ ] `rmsk.bed.gz.tbi` - tabix索引

### 输出文件（每个个体）
- `{individual}_annotated.tsv` - TE注释
- `{individual}_transcript_expression.tsv` - 转录本表达
- `{individual}_gene_expression.tsv` - 基因表达
- `{individual}_te_promoter_transcripts.tsv` - TE-promoter转录本（重要！）
- `{individual}_summary.json` - 个体摘要

### 批量输出
- `group_{N}_batch_summary.tsv` - 组别汇总统计

## 常用参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--workers` | 1 | 并行worker数（建议：12-40） |
| `--min-mapq` | 255 | STAR=255, HISAT2=60 |
| `--verbose` | False | 详细输出 |

## 故障排查

```bash
# 检查BAM文件和索引
ls -lh /path/to/bam_files/*.bam
ls -lh /path/to/bam_files/*.bam.bai

# 为BAM文件建立索引
samtools index sample.bam

# 检查tabix索引
ls -lh ambrosia/rmsk.bed.gz*

# 重建tabix索引
tabix -f -p bed ambrosia/rmsk.bed.gz

# 测试单个样本（不并行）
python examples/batch_ambrosia_population.py \
    --csv gtf-to-bam.csv \
    --target-group 10 \
    --bam-dir /path/to/bam_files \
    --gtf-file /path/to/10.gtf \
    --rmsk-bed ambrosia/rmsk.bed.gz \
    --output-dir ./test_output \
    --workers 1 \
    --verbose
```

## 与旧版对比

| 操作 | 旧版 | 新版 |
|------|------|------|
| 注释 | `rmskhg38_annotate_gtf_final.py` | 自动集成 |
| 定量 | `commandsmax_speed.py` + parallel | 自动集成 |
| 并行 | GNU parallel | Python ProcessPoolExecutor |
| 配置 | 多个bash变量 | 单个CSV文件 |
| 输出 | 分散的文本文件 | 结构化TSV + JSON |

## 性能建议

- **小规模测试**（<10个样本）: `--workers 1-4`
- **中等规模**（10-50个样本）: `--workers 8-12`
- **大规模**（>50个样本）: `--workers 20-40`

**内存估算**: 每个worker约需2-4GB内存

## 更多信息

详细文档: [examples/POPULATION_WORKFLOW_GUIDE.md](POPULATION_WORKFLOW_GUIDE.md)

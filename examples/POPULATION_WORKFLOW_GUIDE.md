# TEProf2 v2.0 豚草种群数据处理指南

本指南说明如何使用TEProf2 v2.0处理豚草（Ambrosia artemisiifolia）种群数据。

## 数据特点

- **BAM文件**: 个体水平（每个个体一个BAM文件）
- **GTF文件**: 种群水平（同一组别的所有个体共享一个GTF文件）
- **分组**: 通过CSV文件指定个体与组别的对应关系

## 工作流程

### 1. 构建TE参考库

首先需要为豚草基因组构建TE参考库。

```bash
# 运行建库脚本
bash examples/build_te_reference_v2.sh ambrosia

# 这将生成以下文件:
# - ambrosia/rmsk.bed.gz          # RepeatMasker注释（bgzipped + tabix）
# - ambrosia/rmsk.bed.gz.tbi      # tabix索引
# - ambrosia/repeat_classification.tsv  # TE分类映射
# - ambrosia/ambrosia_TElib.fa    # TE库文件
```

**重要**: 确保 `rmsk.bed.gz` 和 `rmsk.bed.gz.tbi` 都存在，TEProf2需要tabix索引来快速查询。

### 2. 准备输入文件

#### 2.1 CSV配置文件

创建一个CSV文件（例如 `gtf-to-bam.csv`），格式如下：

```csv
individual_id,group
sample001,10
sample002,10
sample003,10
sample004,11
sample005,11
```

- 第一列: 个体ID（对应BAM文件名，不含.bam后缀）
- 第二列: 组别编号

#### 2.2 BAM文件

确保所有BAM文件：
- 已排序
- 已建立索引（.bam.bai文件）
- 文件名格式: `{individual_id}.bam`

```bash
# 检查BAM文件
ls /path/to/bam_files/*.bam
ls /path/to/bam_files/*.bam.bai

# 如果缺少索引，运行:
samtools index sample001.bam
```

#### 2.3 GTF文件

每个组别一个GTF文件，例如：
- `population_10.gtf` - 组别10的种群水平GTF
- `population_11.gtf` - 组别11的种群水平GTF

### 3. 运行批量处理

#### 3.1 基本用法

```bash
python examples/batch_ambrosia_population.py \
    --csv gtf-to-bam.csv \
    --target-group 10 \
    --bam-dir /path/to/bam_files \
    --gtf-file /path/to/population_10.gtf \
    --rmsk-bed ambrosia/rmsk.bed.gz \
    --repeat-classification ambrosia/repeat_classification.tsv \
    --output-dir ./results_group10 \
    --workers 12
```

#### 3.2 参数说明

| 参数 | 说明 | 必需 |
|------|------|------|
| `--csv` | CSV配置文件路径 | 是 |
| `--target-group` | 要处理的目标组别 | 是 |
| `--bam-dir` | BAM文件目录 | 是 |
| `--gtf-file` | 种群水平GTF文件 | 是 |
| `--rmsk-bed` | RepeatMasker BED文件（bgzipped + tabix） | 是 |
| `--repeat-classification` | TE分类映射文件 | 否 |
| `--gencode-plus` | Gencode字典（+链） | 否 |
| `--gencode-minus` | Gencode字典（-链） | 否 |
| `--output-dir` | 输出目录 | 否（默认: ./results） |
| `--min-mapq` | 最小mapping quality | 否（默认: 255） |
| `--workers` | 并行worker数量 | 否（默认: 1） |
| `--verbose` | 详细输出 | 否 |

#### 3.3 并行处理

使用 `--workers` 参数控制并行度：

```bash
# 使用12个并行worker
python examples/batch_ambrosia_population.py \
    --csv gtf-to-bam.csv \
    --target-group 10 \
    --bam-dir /path/to/bam_files \
    --gtf-file /path/to/population_10.gtf \
    --rmsk-bed ambrosia/rmsk.bed.gz \
    --output-dir ./results_group10 \
    --workers 12
```

**建议**:
- 对于服务器: `--workers 12-40`（根据CPU核心数）
- 对于本地测试: `--workers 1-4`

### 4. 输出文件

对于每个个体，会生成以下文件：

```
results_group10/
├── sample001/
│   ├── sample001_annotated.tsv                    # 注释结果
│   ├── sample001_transcript_expression.tsv        # 转录本表达
│   ├── sample001_gene_expression.tsv              # 基因表达
│   ├── sample001_te_promoter_transcripts.tsv      # TE-promoter转录本
│   └── sample001_summary.json                     # 个体摘要
├── sample002/
│   └── ...
└── group_10_batch_summary.tsv                     # 批量摘要
```

#### 4.1 主要输出文件说明

**1. `{individual}_annotated.tsv`**
- 转录本的TE注释信息
- 包含: transcript_id, has_te_promoter, n_te_overlaps等

**2. `{individual}_transcript_expression.tsv`**
- 转录本水平的表达定量
- 包含: transcript_id, count, TPM, FPKM, coverage等

**3. `{individual}_te_promoter_transcripts.tsv`**
- 合并了注释和表达信息的TE-promoter转录本
- 这是最重要的结果文件

**4. `group_{N}_batch_summary.tsv`**
- 所有个体的汇总统计
- 包含: 转录本数、TE-promoter比例、表达量等

### 5. 处理多个组别

如果有多个组别需要处理，可以使用循环：

```bash
#!/bin/bash

# 定义组别列表
GROUPS="10 11 12 13"

# 循环处理每个组别
for GROUP in $GROUPS; do
    echo "处理组别 $GROUP..."

    python examples/batch_ambrosia_population.py \
        --csv gtf-to-bam.csv \
        --target-group $GROUP \
        --bam-dir /path/to/bam_files \
        --gtf-file /path/to/population_${GROUP}.gtf \
        --rmsk-bed ambrosia/rmsk.bed.gz \
        --output-dir ./results_group${GROUP} \
        --workers 12

    echo "组别 $GROUP 完成！"
done

echo "所有组别处理完成！"
```

### 6. 常见问题

#### Q1: 提示"BAM index not found"

**解决方案**:
```bash
# 为所有BAM文件建立索引
cd /path/to/bam_files
for bam in *.bam; do
    samtools index $bam
done
```

#### Q2: 提示"tabix index not found"

**解决方案**:
```bash
# 重新建立tabix索引
cd ambrosia/
tabix -f -p bed rmsk.bed.gz
```

#### Q3: 内存不足

**解决方案**:
- 减少 `--workers` 数量
- 使用更大内存的节点
- 分批处理个体

#### Q4: 某些个体处理失败

**解决方案**:
- 检查失败个体的BAM文件是否损坏
- 查看个体输出目录中的日志
- 使用 `--verbose` 参数获取详细错误信息

### 7. 性能优化

#### 7.1 并行处理

```bash
# 最优worker数 = CPU核心数 / 2
# 例如40核服务器，使用20个worker
--workers 20
```

#### 7.2 MAPQ设置

```bash
# STAR比对: 使用255（唯一比对）
--min-mapq 255

# HISAT2比对: 使用60
--min-mapq 60
```

### 8. 与旧版TEProf2的主要区别

| 特性 | 旧版TEProf2 | 新版TEProf2 v2.0 |
|------|-------------|------------------|
| 语言 | Bash + Python 2.7 | Python 3.12+ |
| 架构 | 脚本集合 | 模块化包 |
| 并行 | GNU parallel | ProcessPoolExecutor |
| 配置 | 命令行参数 | Config类 + 命令行 |
| 输出 | 文本文件 | TSV + JSON |
| 依赖管理 | 手动 | pyproject.toml |
| 类型检查 | 无 | 完整类型提示 |

### 9. 下一步分析

处理完成后，可以进行：

1. **差异表达分析**: 比较不同组别的TE-promoter转录本表达
2. **TE家族分析**: 统计不同TE家族的promoter活性
3. **可视化**: 绘制TE-promoter转录本的表达热图
4. **功能富集**: 分析TE-promoter基因的功能

## 技术支持

如有问题，请：
1. 查看日志文件
2. 使用 `--verbose` 参数获取详细信息
3. 检查输入文件格式
4. 提交issue到GitHub仓库

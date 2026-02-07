#!/bin/bash
set -e

# =================================================================================
#    Pipeline: 为TEProf2 v2.0构建普通豚草的De Novo TE参考库
#
#   功能:
#   1. 使用内置的基因组路径，结合混合策略 (RepeatModeler + DeepTE) 构建TE库
#   2. 自动运行RepeatMasker并生成TEProf2 v2.0所需的全套参考文件
#
#   TEProf2 v2.0 需要的参考文件:
#   - rmsk.bed.gz: RepeatMasker注释的BED文件（bgzipped + tabix索引）
#   - repeat_classification.tsv: TE分类映射文件（可选）
#   - gencode_plus.dic / gencode_minus.dic: 基因编码字典（如果有基因注释）
#
# =================================================================================

echo "--- 开始为TEProf2 v2.0构建TE参考文件 ---"

# --- CONFIGURATION ---
# =================================================================================
# --- 输入与输出 ---
GENOME_FASTA="/strage_151/data_182/home/lizx/data/Ambrosia_genome_dir/Ambrosia_artemisiifolia_genome.fa"
PREFIX=$1 # 输出文件的前缀
THREADS=40

# --- 工具绝对路径 ---
PATH_REPEATMODELER="/strage_151/data_182/home/lizx/miniconda3/envs/repeatmodeler_env/bin"
PATH_DEEPTE_SCRIPT="/strage_151/data_182/home/lizx/DeepTE/DeepTE-master/DeepTE.py"
PATH_DEEPTE_PYTHON="/strage_151/data_182/home/lizx/miniconda3/envs/deepte_env/bin/python"
PATH_DEEPTE_MODEL_DIR="/strage_151/data_182/home/lizx/DeepTE/lib/Plants_model"
PATH_BASE_PLANT_LIB="/strage_151/data_182/home/lizx/DeepTE/Viridiplantae_v4.0.fasta"

# --- 检查输入参数 ---
if [ -z "$PREFIX" ]; then
  echo "错误: 请提供一个输出文件前缀作为参数！"
  echo "用法: $0 <prefix>"
  exit 1
fi
# =================================================================================

# --- 阶段一：构建TE共有序列库 ---
# =================================================================================
echo "--- 阶段一：开始构建TE共有序列库 ---"

# 1. 准备工作目录
mkdir -p "$PREFIX" && cd "$PREFIX"

# 2. 复制基因组文件
echo "--> 正在复制基因组文件..."
cp "$GENOME_FASTA" ./
LOCAL_GENOME=$(basename "$GENOME_FASTA")

# 3. 运行RepeatModeler
echo "--> 正在运行 RepeatModeler..."
${PATH_REPEATMODELER}/BuildDatabase -name "$PREFIX" "$LOCAL_GENOME"
${PATH_REPEATMODELER}/RepeatModeler -database "$PREFIX" -threads "$THREADS" > repeatmodeler-stdout.log 2>&1

# 4. 分离已知与未知TE
echo "--> 正在分离已知与未知TE..."
seqkit grep -vnrp '#Unknown' "${PREFIX}-families.fa" > "${PREFIX}_knownTE.fa"
seqkit grep -nrp '#Unknown' "${PREFIX}-families.fa" > "${PREFIX}_unknownTE.fa"

# 5. 使用DeepTE分类未知TE
echo "--> 正在使用 DeepTE 分类未知TE..."
mkdir -p working_dir output_dir
${PATH_DEEPTE_PYTHON} "$PATH_DEEPTE_SCRIPT" -d working_dir -o output_dir -i "${PREFIX}_unknownTE.fa" -sp P -m_dir "$PATH_DEEPTE_MODEL_DIR" > deepte-stdout.log 2>&1

# 6. 后处理并合并成最终库
echo "--> 正在合并所有TE序列以生成最终库..."
awk '/^>/ {sub(" ", "#", $1); sub("\t", "/", $1); print $1; next} {print}' output_dir/opt_DeepTE.fasta > "${PREFIX}_unknownTE_classified.fa"
cat "$PATH_BASE_PLANT_LIB" "${PREFIX}_knownTE.fa" "${PREFIX}_unknownTE_classified.fa" > "${PREFIX}_TElib.fa"

echo "--- 阶段一完成！最终TE库文件为: ${PREFIX}_TElib.fa ---"

# --- 阶段二：生成TEProf2 v2.0所需的参考文件 ---
# =================================================================================
echo "--- 阶段二：开始生成TEProf2 v2.0所需的参考文件 ---"

# 7. 使用新生成的库运行RepeatMasker
echo "--> 正在使用新库运行 RepeatMasker..."
RepeatMasker -lib "${PREFIX}_TElib.fa" -gff -dir . -pa "$THREADS" "$LOCAL_GENOME"

# 8. 生成TEProf2参考文件1: 注释位置的BED文件（bgzipped + tabix）
echo "--> 正在生成注释位置的BED文件 (rmsk.bed.gz)..."
# 从RepeatMasker的.out文件生成BED格式
# 格式: chrom start end name score strand
awk 'NR>3 {print $5"\t"$6-1"\t"$7"\t"$10"\t"$1"\t"$9}' "${LOCAL_GENOME}.out" | \
    sort -k1,1 -k2,2n > rmsk.bed

# 压缩并建立索引
bgzip -f rmsk.bed
tabix -f -p bed rmsk.bed.gz
echo "--> rmsk.bed.gz 和 rmsk.bed.gz.tbi 已生成"

# 9. 生成TEProf2参考文件2: TE分类映射文件
echo "--> 正在生成TE分类映射文件 (repeat_classification.tsv)..."
# 从TE库的FASTA头部提取分类信息
# 格式: repeat_name  repeat_class  repeat_family
grep ">" "${PREFIX}_TElib.fa" | \
    sed 's/>//' | \
    awk -F'#' '{
        name=$1
        classification=$2
        # 分割class/family
        split(classification, parts, "/")
        class=parts[1]
        family=(parts[2] != "" ? parts[2] : parts[1])
        print name"\t"class"\t"family
    }' | \
    sort | uniq > repeat_classification.tsv

echo "--> repeat_classification.tsv 已生成"

# 10. 生成使用说明
echo ""
echo "--- Pipeline全部完成！ ---"
echo ""
echo "生成的TEProf2 v2.0参考文件:"
echo "  1. RepeatMasker注释: $(pwd)/rmsk.bed.gz"
echo "  2. TE分类映射:       $(pwd)/repeat_classification.tsv"
echo "  3. TE库文件:         $(pwd)/${PREFIX}_TElib.fa"
echo ""
echo "使用方法:"
echo "  python batch_ambrosia_population.py \\"
echo "      --csv gtf-to-bam.csv \\"
echo "      --target-group 10 \\"
echo "      --bam-dir /path/to/bam_files \\"
echo "      --gtf-file /path/to/population_10.gtf \\"
echo "      --rmsk-bed $(pwd)/rmsk.bed.gz \\"
echo "      --repeat-classification $(pwd)/repeat_classification.tsv \\"
echo "      --output-dir ./results \\"
echo "      --workers 12"
echo ""
echo "注意事项:"
echo "  - rmsk.bed.gz 必须有对应的 .tbi 索引文件"
echo "  - 如果有基因注释，还需要准备 gencode_plus.dic 和 gencode_minus.dic"
echo "  - BAM文件必须已排序并建立索引（.bam.bai）"
echo ""
# =================================================================================

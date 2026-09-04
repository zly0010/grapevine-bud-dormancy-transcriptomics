# Historical execution records confirm this parser was run during early legacy processing.
# Only local path literals were sanitized for this public copy. Historical executed
# source bytes were not separately hashed. The GSE337039 processed branch is historical;
# final GSE337039 validation uses the documented raw-read branch.

"""
统一解析所有数据集，生成标准格式的表达矩阵
"""
import os
import csv
import gzip

RAW_DIR = r"data/legacy/data/raw"
PROCESSED_DIR = r"data/legacy/data/processed"

def merge_gse124820():
    """合并 GSE124820 的 individual count files"""
    print("\n处理 GSE124820 (合并 individual counts)...")
    counts_dir = os.path.join(RAW_DIR, "GSE124820", "counts")

    # 读取所有 count 文件
    all_counts = {}
    sample_names = []

    for f in sorted(os.listdir(counts_dir)):
        if f.endswith('.counts.txt.gz'):
            sample_name = f.replace('_Aligned.counts.txt.gz', '')
            sample_names.append(sample_name)

            filepath = os.path.join(counts_dir, f)
            with gzip.open(filepath, 'rt') as fh:
                for line in fh:
                    parts = line.strip().split('\t')
                    if len(parts) == 2:
                        gene_id = parts[0]
                        count = parts[1]
                        if gene_id not in all_counts:
                            all_counts[gene_id] = {}
                        all_counts[gene_id][sample_name] = count

    # 写入合并后的矩阵
    output_dir = os.path.join(PROCESSED_DIR, "GSE124820")
    os.makedirs(output_dir, exist_ok=True)

    output_file = os.path.join(output_dir, "counts_matrix.txt")
    with open(output_file, 'w') as f:
        # 写表头
        f.write("gene_id\t" + "\t".join(sample_names) + "\n")
        # 写数据
        for gene_id in sorted(all_counts.keys()):
            values = [all_counts[gene_id].get(s, "0") for s in sample_names]
            f.write(f"{gene_id}\t" + "\t".join(values) + "\n")

    print(f"  基因数: {len(all_counts)}")
    print(f"  样本数: {len(sample_names)}")
    print(f"  保存: {output_file}")
    return True

def parse_csv_counts(gse_id, filename, delimiter=','):
    """解析 CSV 格式的 count 文件"""
    print(f"\n处理 {gse_id}...")
    filepath = os.path.join(RAW_DIR, gse_id, filename)

    output_dir = os.path.join(PROCESSED_DIR, gse_id)
    os.makedirs(output_dir, exist_ok=True)

    with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
        # 读取表头
        header = f.readline().strip().split(delimiter)
        sample_names = [h.strip('"') for h in header[1:]]

        # 读取数据
        genes = []
        counts = []
        for line in f:
            parts = line.strip().split(delimiter)
            if len(parts) > 1:
                gene_id = parts[0].strip('"')
                values = [v.strip('"') for v in parts[1:]]
                genes.append(gene_id)
                counts.append(values)

    # 写入标准格式
    output_file = os.path.join(output_dir, "counts_matrix.txt")
    with open(output_file, 'w') as f:
        f.write("gene_id\t" + "\t".join(sample_names) + "\n")
        for i, gene_id in enumerate(genes):
            f.write(f"{gene_id}\t" + "\t".join(counts[i]) + "\n")

    print(f"  基因数: {len(genes)}")
    print(f"  样本数: {len(sample_names)}")
    print(f"  保存: {output_file}")
    return True

def parse_tab_counts(gse_id, filename):
    """解析 Tab 分隔的 count 文件"""
    print(f"\n处理 {gse_id}...")
    filepath = os.path.join(RAW_DIR, gse_id, filename)

    output_dir = os.path.join(PROCESSED_DIR, gse_id)
    os.makedirs(output_dir, exist_ok=True)

    with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
        # 读取表头
        header = f.readline().strip().split('\t')
        sample_names = header[1:]  # 第一列是 gene_id

        # 读取数据
        genes = []
        counts = []
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) > 1:
                gene_id = parts[0]
                values = parts[1:]
                genes.append(gene_id)
                counts.append(values)

    # 写入标准格式
    output_file = os.path.join(output_dir, "counts_matrix.txt")
    with open(output_file, 'w') as f:
        f.write("gene_id\t" + "\t".join(sample_names) + "\n")
        for i, gene_id in enumerate(genes):
            f.write(f"{gene_id}\t" + "\t".join(counts[i]) + "\n")

    print(f"  基因数: {len(genes)}")
    print(f"  样本数: {len(sample_names)}")
    print(f"  保存: {output_file}")
    return True

def parse_gse184114():
    """解析 GSE184114 (两个 count 矩阵合并)"""
    print(f"\n处理 GSE184114 (合并 Acclimation + Deacclimation)...")

    output_dir = os.path.join(PROCESSED_DIR, "GSE184114")
    os.makedirs(output_dir, exist_ok=True)

    # 读取 Acclimation
    accl_file = os.path.join(RAW_DIR, "GSE184114", "GSE184114_Acclimation_raw_gene_count_matrix.txt")
    deaccl_file = os.path.join(RAW_DIR, "GSE184114", "GSE184114_Deacclimation_raw_gene_count_matrix.txt")

    all_genes = {}
    all_samples = []

    for filepath, prefix in [(accl_file, "Accl_"), (deaccl_file, "Deaccl_")]:
        with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
            header = f.readline().strip().split('\t')
            samples = [f"{prefix}{h}" for h in header]
            all_samples.extend(samples)

            for line in f:
                parts = line.strip().split('\t')
                if len(parts) > 1:
                    gene_id = parts[0]
                    values = parts[1:]
                    if gene_id not in all_genes:
                        all_genes[gene_id] = {}
                    for i, val in enumerate(values):
                        all_genes[gene_id][samples[i]] = val

    # 写入
    output_file = os.path.join(output_dir, "counts_matrix.txt")
    with open(output_file, 'w') as f:
        f.write("gene_id\t" + "\t".join(all_samples) + "\n")
        for gene_id in sorted(all_genes.keys()):
            values = [all_genes[gene_id].get(s, "0") for s in all_samples]
            f.write(f"{gene_id}\t" + "\t".join(values) + "\n")

    print(f"  基因数: {len(all_genes)}")
    print(f"  样本数: {len(all_samples)}")
    print(f"  保存: {output_file}")
    return True

if __name__ == "__main__":
    os.makedirs(PROCESSED_DIR, exist_ok=True)

    # 处理各数据集
    merge_gse124820()
    parse_csv_counts("GSE273240", "GSE273240_tetralone_ABA_gene_count.csv")
    parse_gse184114()
    parse_csv_counts("GSE337039", "GSE337039_geo_processed_read_counts.csv")
    parse_tab_counts("GSE277812", "GSE277812_raw_counts.txt")

    print(f"\n{'='*60}")
    print("所有数据集解析完成！")
    print(f"{'='*60}")

    # 统计
    print("\n处理结果:")
    for gse_id in ["GSE124820", "GSE273240", "GSE184114", "GSE337039", "GSE277812"]:
        output_dir = os.path.join(PROCESSED_DIR, gse_id)
        if os.path.exists(output_dir):
            files = os.listdir(output_dir)
            print(f"  {gse_id}: {', '.join(files)}")

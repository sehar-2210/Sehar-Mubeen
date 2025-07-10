# analysis_SNP.py
import pandas as pd
from collections import defaultdict
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr
# File paths
gff_file = "Saccharomyces_cerevisiae.R64-1-1.61.gff3"
vcf_file = "saccharomyces_cerevisiae.vcf"
expression_file = "merged_gene_expression.csv"
# Step 1: Extract exon coordinates from GFF3
def extract_exons(gff_file):
    exons = defaultdict(list)
    with open(gff_file) as file:
        for line in file:
            if line.startswith("#"):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9 or parts[2].lower() != "exon":
                continue
            chrom = parts[0]
            start = int(parts[3])
            end = int(parts[4])
            attributes = parts[8]
            gene_id = None
            for item in attributes.split(";"):
                if "Parent" in item or "GeneID" in item:
                    gene_id = item.split("=")[-1]
                    break
            if gene_id:
                exons[gene_id].append((chrom, start, end))
    return exons
# Step 2: Extract SNPs from VCF
def extract_snps(vcf_file):
    snps = defaultdict(list)
    with open(vcf_file) as file:
        for line in file:
            if line.startswith("#"):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 2:
                continue
            chrom = parts[0]
            pos = int(parts[1])
            snps[chrom].append(pos)
    return snps
# Step 3: Calculate SNP density
def calculate_snp_density(exons, snps):
    results = []
    for gene_id, exon_list in exons.items():
        total_exon_length = 0
        total_snps = 0
        for chrom, start, end in exon_list:
            exon_length = end - start + 1
            total_exon_length += exon_length
            snp_count = sum(1 for pos in snps.get(chrom, []) if start <= pos <= end)
            total_snps += snp_count
        if total_exon_length > 0:
            density = total_snps / total_exon_length
            results.append((gene_id, density))
    return pd.DataFrame(results, columns=["GeneID", "SNP_Density"])
# Step 4: Merge SNP density with gene expression
def merge_expression(snp_df, expression_file):
    exp_df = pd.read_csv(expression_file)
    # Rename and fix columns based on your file
    exp_df.rename(columns={"gene_id": "GeneID"}, inplace=True)
    exp_df["Expression"] = exp_df["fpkm_Sample1"]  
    exp_df = exp_df[["GeneID", "Expression"]]
    merged = pd.merge(snp_df, exp_df, on="GeneID")
    return merged
# Step 5: Visualizations + Pearson correlation
def create_visualizations(df):
    #  Skip if not enough data
    if df.shape[0] < 2:
        print("\n Not enough data for correlation or plots. Skipping visualizations.")
        df.to_csv("debug_empty_merged_df.csv", index=False)
        return
    # 1. Scatter Plot
    plt.figure(figsize=(8, 6))
    sns.scatterplot(data=df, x="Expression", y="SNP_Density", color='purple')
    plt.title("SNP Density vs Gene Expression")
    plt.xlabel("Gene Expression")
    plt.ylabel("SNP Density")
    plt.tight_layout()
    plt.savefig("scatter_plot.png")
    plt.close()
    # 2. Bar Plot (Top 20 Genes by SNP Density)
    top20 = df.sort_values(by="SNP_Density", ascending=False).head(20)
    plt.figure(figsize=(12, 6))
    sns.barplot(data=top20, x="GeneID", y="SNP_Density", palette="viridis")
    plt.xticks(rotation=90)
    plt.title("Top 20 Genes by SNP Density")
    plt.xlabel("Gene ID")
    plt.ylabel("SNP Density")
    plt.tight_layout()
    plt.savefig("bar_plot_top20_snp_density.png")
    plt.close()
    # 3. Heatmap of Correlation
    plt.figure(figsize=(6, 4))
    corr_matrix = df[["Expression", "SNP_Density"]].corr()
    sns.heatmap(corr_matrix, annot=True, cmap="coolwarm", fmt=".2f")
    plt.title("Correlation Matrix")
    plt.tight_layout()
    plt.savefig("correlation_heatmap.png")
    plt.close()
    # 4. Pearson Correlation
    corr, pval = pearsonr(df["Expression"], df["SNP_Density"])
    print(f"\n Pearson correlation:")
    print(f"   r = {corr:.3f}")
    print(f"   p-value = {pval:.3e}")
# Main Execution
def main():
    print(" Step 1: Parsing GFF3 file...")
    exons = extract_exons(gff_file)
    print(" Step 2: Parsing VCF file...")
    snps = extract_snps(vcf_file)
    print(" Step 3: Calculating SNP densities...")
    snp_df = calculate_snp_density(exons, snps)
    print(" Step 4: Merging with expression data...")
    merged_df = merge_expression(snp_df, expression_file)
    merged_df.to_csv("merged_results.csv", index=False)
    print(" Step 5: Creating visualizations...")
    create_visualizations(merged_df)
    print("\n Analysis complete! Results and plots saved.")
if __name__ == "__main__":
    main()

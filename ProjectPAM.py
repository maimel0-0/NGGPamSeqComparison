import re
from pathlib import Path

import pandas as pd
from Bio.SeqIO.FastaIO import SimpleFastaParser

import matplotlib.pyplot as plt
import seaborn as sns

from scipy.stats import ttest_ind

#parse the fna files
def parse_fasta(fasta_file_path, genome_type: str):
    ngg_pattern = re.compile(r'(?=([ACGT]GG))')
    results = []
    with open(fasta_file_path, "r") as file:
        for header, sequence in SimpleFastaParser(file):
            genome_id = header.split()[0]
            genome_str = sequence.upper()

#ngg count
            total_triplets = len(genome_str) - 2
            if total_triplets <= 0:
                continue

            ngg_count = len(ngg_pattern.findall(genome_str))
            ngg_percentage = (ngg_count / total_triplets) * 100

            results.append({
                "Genome_ID": genome_id,
                "Type": genome_type,
                "Genome_Length": len(genome_str),
                "NGG_Count": ngg_count,
                "NGG_Percentage": round(ngg_percentage, 4)
            })

    return results


def main():

    all_results = []

    # Analyze bacteria
    for file in Path("./data/bacteria").glob("*.fna"):
        all_results.extend(parse_fasta(file, "Bacteria"))

    # Analyze phages
    for file in Path("./data/phages").glob("*.fna"):
        all_results.extend(parse_fasta(file, "Phage"))

    df = pd.DataFrame(all_results)

    print(df)

#get the averages
    average_df = (
        df.groupby("Type")["NGG_Percentage"]
        .mean()
        .reset_index()
    )

    print(average_df)

# plot the data
    plt.figure(figsize=(6, 5))

    sns.boxplot(
        data=df,
        x="Type",
        y="NGG_Percentage"
    )

    plt.title("Distribution of NGG Percentages")
    plt.xlabel("Genome Type")
    plt.ylabel("NGG Percentage")

    plt.show()

    plt.figure(figsize=(6, 5))

    sns.boxplot(
        data=df,
        x="Type",
        y="NGG_Percentage"
    )

    sns.stripplot(
        data=df,
        x="Type",
        y="NGG_Percentage",
        color="black"
    )

    plt.title("NGG Percentages by Genome Type")

    plt.show()

    plt.savefig("NGG_Boxplot.png", dpi=300)
    plt.show()


    bacteria = df[df["Type"] == "Bacteria"]["NGG_Percentage"]
    phages = df[df["Type"] == "Phage"]["NGG_Percentage"]

    # perform t-test
    t_stat, p_value = ttest_ind(bacteria, phages)

    print("\nT-test Results")
    print("T-statistic:", t_stat)
    print("P-value:", p_value)

    if p_value < 0.05:
        print("The difference is statistically significant.")
    else:
        print("The difference is not statistically significant.")







main()



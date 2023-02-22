import pandas as pd
import re
import logging
import os
import subprocess
from io import StringIO

# """
# Identify interesting regions for strain-specific primers from VCF file
# """

# Specify the path to the VCF file
vcf_file = "Alignment/Protease1_Bvel10075.test.vcf"

# retrieve INFO column key:value pairs from metadata
f = open(vcf_file, "r")
line = f.readline().rstrip()
new_col_names = []
while "##" in line:
    if "##INFO=" in line:
        new_col_names.append(re.search("ID=[A-Z]+", line).group().lstrip("ID="))
    line = f.readline()

# Define the column names
col_names = line.split("\t")[0:9]
f.close()

# Read in the VCF file as a dataframe
df = pd.read_csv(vcf_file, sep='\t', skiprows = 29, names=col_names, usecols=range(0,9))

# Split the INFO column into separate fields
def splitinfo(str):
    str = str.split(";")
    info = {var:value for var, value in [pair.split("=") for pair in str if "=" in pair]}
    return info

df = pd.concat([df, df["INFO"].apply(splitinfo).apply(pd.Series)], axis=1)
df = df.drop("INFO", axis=1)

# Filter the variant calls on quality metrics
df[["DP", "IMF", "QUAL", "VDB", "MQ", "MQSB"]] = df[["DP", "IMF", "QUAL", "VDB", "MQ", "MQSB"]].apply(pd.to_numeric)
df = df[(df["DP"] > 1000) & (df["QUAL"] > 30) & (df["MQ"] > 20) & (df["MQSB"] > 0.5) & (df["VDB"] > 0.8) & ((df["IMF"] > 0.80) | pd.isnull(df["IMF"]))]
df = df.sort_values(["DP", "IMF", "QUAL"], ascending=[False, False, False])
df.to_csv("Sorted_VCF.tsv", sep="\t", index=False)

# Retrieve interesting primer positions: more than 2 SNPs in <25 bp (=> multiple SNPs on primer)
pos_primer_pos = df.sort_values(["POS"])
pos_primer_pos["POS_prev1"] = pos_primer_pos["POS"].shift(1)
pos_primer_pos["POS_prev2"] = pos_primer_pos["POS"].shift(2)
pos_primer_pos["POS_prev3"] = pos_primer_pos["POS"].shift(3)
pos_primer_pos = pos_primer_pos[pos_primer_pos["POS"] - pos_primer_pos["POS_prev1"] < 25 ]

pos_primer_pos.to_csv("PossiblePrimers.tsv", sep="\t", index=False)

# """
# Extract Primer targets from reference genome
# """

# retrieve interesting primer loci from reference fasta; write to query fasta file for later BLASTing
with open("tmp.fa", "w") as temp, open("/scratch/brbloemen/SHIME/GMMs/ProteaseIGMMs/VCF/Reference/Bvel10075.fasta", "r") as ref:
    # parse genome and put into single string
    ref_name = ref.readline().strip(">").rstrip()
    genome = "".join(line.rstrip("\n") for line in ref)
    gen_length = len(genome)

    positions = [*map(int, pos_primer_pos["POS"].tolist())]
    extract_pos = []

    start, end = max(positions[0]-100, 0), min(positions[0] + 200, gen_length)
    for index, pos in enumerate(positions):
        
        # first position or the new position is past the end of the previous region -> start new region
        if index == 0 or pos > end:
            start, end = max(pos - 100, 0), min(pos + 200, gen_length)
            extract_pos.append([start, end])

        # the new positions is within range of the previous region -> make region longer
        else:
            end = min(pos + 200, gen_length)
            extract_pos[-1][1] = end
    # extract query sequences (shift to 0-based, vcf positions were 1-based)
    query_seqs = [genome[region[0]-1:region[1]-1] for region in extract_pos]

    for index, region in enumerate(extract_pos):
        query_name = ref_name + " " + str(region[0]-1) + ":" + str(region[1]-1)
        temp.write(f">{query_name}\n")
        temp.write(f"{query_seqs[index]}\n")


# """
# BLAST against variant genomes
# """

# logging.info(f"Blasting upstream and downstream parts of ARG reads: {2*num_reads} sequences")
# Set the input FASTA file
fasta_file = 'tmp.fa'

# Set the BLAST command and arguments
blast_cmd = 'blastn'
blast_db = '/scratch/brbloemen/SHIME/GMMs/ProteaseIGMMs/Assemblies/BvelGMM_proteaseI_isolates/BvelGMM_ProteaseI_BLASTdb.fasta'
blast_threads = 30
blast_args = '-query %s -db %s -outfmt 0 -num_threads %s -max_target_seqs 50 -task megablast' % (fasta_file, blast_db, blast_threads)

# Build the full command
command = '%s %s' % (blast_cmd, blast_args)

# Run the BLAST search
with open("blast_output.txt", "w") as outputf:
    subprocess.run(["bash", "-c", f"module load blast; {command}"], stdout=outputf)
# subprocess.PIPE

# # Read the TSV file from the stdout into a DataFrame
# BLAST_table = pd.read_csv(StringIO(result.stdout.decode('utf-8')), sep='\t', 
# header=None, names=["qseqid","sseqid","pident","length","sacc","staxids","sscinames","sblastname","stitle"])
# logging.info(f"Done BLASTing")

# # remove temporary fasta
# os.remove("tmp.fa")

import pysam
import argparse
import os
import sys
import re
import subprocess
import pandas as pd
import datetime
from io import StringIO
import plotly.express as px
import logging
logging.basicConfig(level = logging.DEBUG,format='%(asctime)s %(message)s',
                    datefmt='%x %X')

parser = argparse.ArgumentParser(description='Retrieve ARG reads, blast upstream and downstream part of the read')
parser.add_argument("input_file", help="SAM/BAM formatted file with alignments to resistance genes")
parser.add_argument("-o","--output", type=str, help="Output filename", default=None)
parser.add_argument("--tid_arg", type=float, help="Minimum Percentage Identity to AMR gene template", default=97.0)
parser.add_argument("-db","--db", type=str, help="database to blast against", default="nt")
parser.add_argument("-t","--threads", type=int, help="Number of threads for BLAST", default=5)
args = parser.parse_args()
if args.output is None:
    args.output = re.search(r'([^/]+)$',args.input_file).group()
    args.output = re.sub('.bam', '', args.output)

logging.info(f"Parsing alignment to ARG database: {args.input_file}")
# Open the BAM file in read mode
with pysam.AlignmentFile(args.input_file, "rb", threads=args.threads) as bam_file:
    refseq_lens = {ref["SN"]:ref["LN"] for ref in bam_file.header.to_dict()["SQ"]}

    num_reads=0
    with open("tmp.fa", "w") as temp:
        # Iterate over all reads in the BAM file
        for read in bam_file:
            # Check if the read is mapped
            if not read.is_unmapped:
                # Get the reference name and position of the alignment
                ref_len = refseq_lens[read.reference_name]
                query_start = read.query_alignment_start
                query_end = read.query_alignment_end

                # rename query sequences to contain information about ARG, read identity and timing
                query_name = [read.query_name.split(sep=" ")[i] for i in (0,5)]
                query_name.append(read.reference_name)
                query_name = ";".join(query_name)

                # Get the sequence of the read
                seq = read.query_sequence

                # Get template identity (percentage):
                tid = 100*read.get_cigar_stats()[0][7]/ref_len

                # Check if the read fully covers the template sequence
                if tid > args.tid_arg:
                    num_reads += 1
                    # Perform a BLAST search on the upstream and downstream sequences
                    temp.write(f">{query_name};upstr\n")
                    temp.write(f"{seq[:query_start]}\n")
                    temp.write(f">{query_name};dwnstr\n")
                    temp.write(f"{seq[query_end:]}\n")

"""
BLAST part
"""
logging.info(f"Blasting upstream and downstream parts of ARG reads: {2*num_reads} sequences")
# Set the input FASTA file
fasta_file = 'tmp.fa'

# Set the BLAST command and arguments
blast_cmd = 'blastn'
blast_args = '-query %s -db %s -outfmt "6 qseqid sseqid pident length sacc staxids sscinames sblastname stitle" -num_threads %s -max_target_seqs 1 -task megablast' % (fasta_file, args.db, args.threads)

# Build the full command
command = '%s %s' % (blast_cmd, blast_args)

# Run the BLAST search
result = subprocess.run(["bash", "-c", f"module load blast; {command}"], stdout=subprocess.PIPE)

# Read the TSV file from the stdout into a DataFrame
BLAST_table = pd.read_csv(StringIO(result.stdout.decode('utf-8')), sep='\t', 
header=None, names=["qseqid","sseqid","pident","length","sacc","staxids","sscinames","sblastname","stitle"])
logging.info(f"Done BLASTing")

# remove temporary fasta
os.remove("tmp.fa")

# Split the queryID column
BLAST_table[["ReadID", "Time", "ARG", "Rel_pos"]] = BLAST_table["qseqid"].str.split(";", expand=True)
BLAST_table.drop("qseqid", axis=1)

# convert Time string to datetime value
def convert_to_datetime(string):
    datetime_string = string.split('=')[1]
    return datetime.datetime.strptime(datetime_string, '%Y-%m-%dT%H:%M:%SZ')

BLAST_table["Time"] = BLAST_table["Time"].apply(convert_to_datetime)

# extract species name
def extract_species_name(string):
    pattern = re.compile(r'[A-Z][a-z]+ [a-z]+')

    species = re.sub(r'\[|\]|\(|\)', '', string)
    species = re.sub('_', ' ', species)
    species = re.sub(r'Synthetic|scf', '', species)
    species = re.sub('E.coli', 'Escherichia coli', species)
    species = re.sub("veillonella", "Veillonella", species)
    species = re.sub("Clostridium difficil+e", "Clostridioides difficile", species)
    if not re.search(pattern, species):
        species = "Unknown"
    else:
        species = re.search(pattern, species).group()

    return species
BLAST_table["Species"] = BLAST_table["stitle"].apply(extract_species_name)

# print BLAST output to file
BLAST_table.to_csv(f"{args.output}_BLAST.tsv", sep='\t', index=True)

# Retain only the best hit per ReadID and Rel_pos
BLAST_table['Matching_bp'] = BLAST_table['pident'] * BLAST_table['length']
grouped = BLAST_table.groupby(['ARG', 'Species', 'ReadID', 'Rel_pos'])

max_index = grouped['Matching_bp'].idxmax()
ARGlinks = BLAST_table.loc[max_index]

# Count number of resistance gene links, optionally filter for up and downstream BLAST hits to be the same

    # # Group the rows by "species" and "ID"
    # grouped = ARGlinks.groupby(['ARG', 'Species', 'ReadID'])

    # # Define a function that returns True if the group contains both "up" and "down" in the "direction" column, and False otherwise
    # def has_both(group):
    #     return set(group['Rel_pos']) == {'upstr', 'dwnstr'}

    # # Use the apply() method to apply the function to each group
    # filtered_df = grouped.apply(has_both)

    # # Select the rows where the function returns True
    # ARGlinks = ARGlinks[filtered_df]

ARGlinks = ARGlinks.groupby(['ARG', 'Species']).sum(['Matching_bp'])

# Write the DataFrame to a TSV file
# ARGlinks.to_csv(sys.stdout, sep='\t', index=True)

# Reset the index and create columns for the 'Species' and 'ARG' columns

ARGlinks = ARGlinks.reset_index()
# ARGlinks['Rgene_species_links'] = ARGlinks['Rgene_species_links']/2
ARGlinks.to_csv(f"{args.output}_ARGlinks.tsv", sep='\t')
# Create the scatterplot
fig = px.scatter(ARGlinks, x="Species", y="ARG", size="Matching_bp", size_max=50)

# Save the figure to a PNG file
fig.write_html(f"{args.output}_ARGlink.html")
fig.write_image(f"{args.output}_ARGlink.png")
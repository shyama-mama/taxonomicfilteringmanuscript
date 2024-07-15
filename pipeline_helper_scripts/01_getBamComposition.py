import pysam
import re
import pandas as pd
import argparse
import textwrap
import sys
import concurrent.futures


TAXID = 'taxid'
SEQID = 'seqid'
FASTQHEADER = "fastqheader"
ORIGINAL_SEQID = "original_seqid"
ORIGINAL_START = "original_position_start"
ORIGINAL_END = "original_position_end"
ORIGINAL_SEQ_LEN = "original_seq_length"
MAPPED_SEQID = "mapped_seqid"
MAPPED_START = "mapped_position_start"
MAPPED_END = "mapped_position_end"
MAPPED_LEN = "mapped_seq_length"
MAPPING_QUAL = "mapping_quality"
IS_UNMAPPED = "is_unmapped"
GARGAMMEL_TAG = "gargammel_tag"
ORIGINAL_TAXID = "original_taxid"

# Make new dataframe
bam_df = pd.DataFrame(columns=[
    FASTQHEADER,
    ORIGINAL_SEQID,
    ORIGINAL_START,
    ORIGINAL_END,
    ORIGINAL_SEQ_LEN,
    MAPPED_SEQID,
    MAPPED_START,
    MAPPED_END,
    MAPPED_LEN,
    MAPPING_QUAL,
    IS_UNMAPPED,
    GARGAMMEL_TAG
])



# Extract information from bam file

def process_read(alignment):
    read_name = alignment.query_name
    read_name_list = read_name.split(':')
    print(read_name_list)
    original_seqid = re.sub(r"^[MFR]_", "", read_name_list[0])
    original_position_start = read_name_list[2]
    original_position_end = read_name_list[3]
    original_seq_length = re.findall(r"\d+", read_name_list[5])[0] if re.findall(r"\d+", read_name_list[5]) else None
    gargammel_tag = re.findall(r"[cbe]", read_name_list[5])[0] if re.findall(r"[cbe]", read_name_list[5]) else None
    mapped_seqid = alignment.reference_name
    mapping_quality = alignment.mapping_quality
    mapped_position_start = alignment.reference_start
    mapped_seq_length = alignment.query_length
    mapped_position_end = alignment.reference_end
    is_unmapped = alignment.is_unmapped
    new_row = {
        FASTQHEADER: read_name,
        ORIGINAL_SEQID: original_seqid,
        ORIGINAL_START: original_position_start,
        ORIGINAL_END: original_position_end,
        ORIGINAL_SEQ_LEN: original_seq_length,
        MAPPED_SEQID: mapped_seqid,
        MAPPED_START: mapped_position_start,
        MAPPED_END: mapped_position_end,
        MAPPED_LEN: mapped_seq_length,
        MAPPING_QUAL: mapping_quality,
        IS_UNMAPPED: is_unmapped,
        GARGAMMEL_TAG: gargammel_tag
    }
    return new_row

def process_read_emperical(alignment):
    read_name = alignment.query_name
    mapped_seqid = alignment.reference_name
    mapping_quality = alignment.mapping_quality
    mapped_position_start = alignment.reference_start
    mapped_seq_length = alignment.query_length
    mapped_position_end = alignment.reference_end
    is_unmapped = alignment.is_unmapped
    new_row = {
        FASTQHEADER: read_name,
        ORIGINAL_SEQID: None,
        ORIGINAL_START: None,
        ORIGINAL_END: None,
        ORIGINAL_SEQ_LEN: None,
        MAPPED_SEQID: mapped_seqid,
        MAPPED_START: mapped_position_start,
        MAPPED_END: mapped_position_end,
        MAPPED_LEN: mapped_seq_length,
        MAPPING_QUAL: mapping_quality,
        IS_UNMAPPED: is_unmapped,
        GARGAMMEL_TAG: None
    }
    return new_row

def make_bam_df(bam_file, bam_df):
    
    # Open the BAM file
    bam = pysam.AlignmentFile(bam_file, "rb")

    # Initialize a counter for the reads
    read_count = 0
    curr = 0


    # Iterate over the alignments
    for alignment in bam:
        read_name = alignment.query_name
        read_name_list = read_name.split(':')
        original_seqid = re.sub(r"^[MFR]_", "", read_name_list[0])
        original_position_start = read_name_list[2]
        original_position_end = read_name_list[3]
        original_seq_length = re.findall(r"\d+", read_name_list[5])[0] if re.findall(r"\d+", read_name_list[5]) else None
        gargammel_tag = re.findall(r"[cbe]", read_name_list[5])[0] if re.findall(r"[cbe]", read_name_list[5]) else None
        mapped_seqid = alignment.reference_name
        mapping_quality = alignment.mapping_quality
        mapped_position_start = alignment.reference_start
        mapped_seq_length = alignment.query_length
        mapped_position_end = alignment.reference_end
        is_unmapped = alignment.is_unmapped
        new_row = {
            FASTQHEADER: read_name,
            ORIGINAL_SEQID: original_seqid,
            ORIGINAL_START: original_position_start,
            ORIGINAL_END: original_position_end,
            ORIGINAL_SEQ_LEN: original_seq_length,
            MAPPED_SEQID: mapped_seqid,
            MAPPED_START: mapped_position_start,
            MAPPED_END: mapped_position_end,
            MAPPED_LEN: mapped_seq_length,
            MAPPING_QUAL: mapping_quality,
            IS_UNMAPPED: is_unmapped,
            GARGAMMEL_TAG: gargammel_tag
        }
        bam_df = bam_df.append(new_row, ignore_index=True)
        read_count += 1

        print(curr)

        diff = pct_done - curr
        if diff > 1:
            curr = pct_done
            print(pct_done)

    # Close the BAM file
    bam.close()

    # Return the read count
    return read_count, bam_df



parser = argparse.ArgumentParser(prog='get_bam_composition.py',
   usage='python %(prog)s [-h] -b <input bam> -t <taxid to seqid map> -o <output filename>',
   formatter_class=argparse.RawDescriptionHelpFormatter,
   description=textwrap.dedent('''\
   author:
     Shyamsundar Ravishankar (shyamsundar.ravishankar@adelaide.edu.au)

   description:
     %(prog)s Script to extract bam information 
   '''))
parser.add_argument("-b", "--input_bam", dest="input_bam", default=None, help="BAM file", required=True)
parser.add_argument("-o", "--output_file", dest="output_file", default="out.csv", help="Path to output file")
parser.add_argument("-e", "--emperical", dest="emperical", action="store_true", help="If this is an emperical sample.", required=False)
args = parser.parse_args()

output_file = args.output_file
emperical = args.emperical

# Get the input args 
# taxid_seqid_map = args.taxid_seqid_map
# taxid_seqid_df = pd.read_csv(taxid_seqid_map, header=None, sep='\t', dtype=str)
# taxid_seqid_df.columns = [TAXID, SEQID]

bam_file = pysam.AlignmentFile(args.input_bam, 'rb') 

# read_count, bam_df = make_bam_df(bam_file, bam_df)
# Set the number of threads
num_threads = 16


# Create a thread pool
with concurrent.futures.ThreadPoolExecutor(max_workers=num_threads) as executor:
    # Submit read processing tasks to the thread pool
    if(emperical):
        read_tasks = [executor.submit(process_read_emperical, alignment) for alignment in bam_file]
    else:
        read_tasks = [executor.submit(process_read, alignment) for alignment in bam_file]

    # Get the results as the tasks complete
    results = [task.result() for task in concurrent.futures.as_completed(read_tasks)]

# Create a DataFrame from the results
bam_df = pd.DataFrame(results)

# Close the BAM file
bam_file.close()

# print(read_count)

# bam_df = pd.merge(bam_df, taxid_seqid_df, left_on=ORIGINAL_SEQID, right_on=SEQID, how='left')
# bam_df = bam_df.drop(SEQID, axis=1)
# bam_df = bam_df.rename({TAXID: ORIGINAL_TAXID}, axis=1)


# print(len(bam_df))

bam_df.to_csv(output_file, index=False)

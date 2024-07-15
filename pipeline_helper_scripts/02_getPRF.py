import pandas as pd 
import argparse
import textwrap
from pathlib import Path

ORIGINAL_TAXID = 'original_taxid'
IS_UNMAPPED = 'is_unmapped'
GARGAMMEL_TAG = 'gargammel_tag'
MAPPED_START_POS = 'mapped_position_start'
MAPPING_QUALITY = 'mapping_quality'
MAPQS = [0, 10, 20, 30]
ENDO = ["e"]
ENDO_CONT = ["e", "c"]

HUMAN = "9606"

BACT_FILT = 0
CONT_FILT = 0
ENDO_FILT = 0




# these are mapped reads ('mapping_position_start != -1') 
# that have the same taxid as the reference

# Translation
# df[IS_UNMAPPED] means it is unmapped
# ~df[IS_UNMAPPED] means it is not unmapped or it is mapped 

def calculate_true_positive(df, tag_list, mode=0, value=0):
    if mode == 0:
        return len(df[(~df[IS_UNMAPPED]) & (df[GARGAMMEL_TAG].isin(tag_list))].index)
    elif mode == 1:
        return len(df[(~df[IS_UNMAPPED]) & (df[MAPPING_QUALITY] > value) & (df[GARGAMMEL_TAG].isin(tag_list))].index)

# These are all reads that aligned to our reference and are not from reference taxa
# For this we remove reads with 'mapping_position_start = -1' (unmapped reads)
def calculate_false_positive(df, tag_list, mode=0, value=0):
    if mode == 0:
        return len(df[(~df[IS_UNMAPPED]) & (~df[GARGAMMEL_TAG].isin(tag_list))].index)
    elif mode == 1:
        return len(df[(~df[IS_UNMAPPED]) & (df[MAPPING_QUALITY] > value) & (~df[GARGAMMEL_TAG].isin(tag_list))].index)

# number of reads that are unmapped but are actually from taxid 
def calculate_false_negative(df, tag_list, mode=0, value=0, false_negative_kraken=0):
    if mode == 0:
        return len(df[(df[IS_UNMAPPED]) & (df[GARGAMMEL_TAG].isin(tag_list))].index)+false_negative_kraken
    elif mode == 1:
        len1 = len(df[(~df[IS_UNMAPPED]) & (df[MAPPING_QUALITY] <= value) & (df[GARGAMMEL_TAG].isin(tag_list))].index)
        len2 = len(df[(df[IS_UNMAPPED]) & (df[GARGAMMEL_TAG].isin(tag_list))].index)
        sum = len1 + len2 + false_negative_kraken
        return sum

# Unmapped reads that are not from reference taxa
def calculate_true_negative(df, tag_list, mode=0, value=0, true_negative_kraken=0):
    if mode == 0:
        return len(df[(df[IS_UNMAPPED]) & (~df[GARGAMMEL_TAG].isin(tag_list))].index) + true_negative_kraken
    elif mode == 1:
        len1 = len(df[(~df[IS_UNMAPPED]) & (df[MAPPING_QUALITY] <= value) & (~df[GARGAMMEL_TAG].isin(tag_list))].index)
        len2 = len(df[(df[IS_UNMAPPED]) & (~df[GARGAMMEL_TAG].isin(tag_list))].index)
        sum = len1 + len2 + true_negative_kraken
        return sum


def calculate_precision(true_positive, false_positive):
    precision = true_positive / (true_positive + false_positive)
    return precision

def calculate_recall(true_positive, false_negative):
    recall = true_positive / (true_positive + false_negative)
    return recall

def calculate_f_measure(precision, recall):
    f_measure = (2 * precision * recall)/(precision + recall)
    return f_measure

def display_all_values(df, tag_list, mode=0, value=0, true_negative_kraken=0, false_negative_kraken=0):
    true_positive = calculate_true_positive(df, tag_list, mode=mode, value=value)
    true_negative = calculate_true_negative(df, tag_list, mode=mode, value=value, true_negative_kraken=true_negative_kraken)
    false_positive = calculate_false_positive(df, tag_list, mode=mode, value=value)
    false_negative = calculate_false_negative(df, tag_list, mode=mode, value=value, false_negative_kraken=false_negative_kraken)

    total = true_negative + true_positive + false_negative + false_positive

    precision = calculate_precision(true_positive, false_positive)
    recall = calculate_recall(true_positive, false_negative)

    f_measure = calculate_f_measure(precision, recall)
    
    return [total, true_positive, true_negative, false_positive, false_negative, precision, recall, f_measure]

# Arguements
parser = argparse.ArgumentParser(prog='kraken_visualise.py',
   usage='python %(prog)s [-h] -i <input file> -e <taxid of endo sequences> -o <output file>',
   formatter_class=argparse.RawDescriptionHelpFormatter,
   description=textwrap.dedent('''\
   author:
     Shyamsundar Ravishankar (shyamsundar.ravishankar@adelaide.edu.au)

   description:
     %(prog)s Script to evaluate calculate precision, recall and f-score from parsed bam file
   '''))
parser.add_argument("-i", "--input_file", dest="input_file", default=None, help="kraken2 output file", required=True)
parser.add_argument("-o", "--output_file", dest="output_file", default="out.csv", help="Path to output file")
parser.add_argument("-s", "--species", dest="species", default=None, help="TaxID of endogenous species", required=True)
parser.add_argument("-t", "--stats_file", dest="stats_file", default=None, help="Path to stats file")

# parser.add_argument("-c", "--include_contamination", dest="include_contamination", default=False, help="If contamination is same as endo sepcies")
args = parser.parse_args()

input_file = args.input_file
name = Path(input_file).stem.replace('.bam_composition', '')


output_file = args.output_file
species = args.species

stats_file = args.stats_file

if stats_file is not None:
    filtered_stats_df = pd.read_csv(stats_file, dtype=int)
    BACT_FILT = filtered_stats_df['bact'].iloc[0]
    CONT_FILT = filtered_stats_df['cont'].iloc[0]
    ENDO_FILT = filtered_stats_df['endo'].iloc[0]

# If  
tag_list = ENDO
true_negative_kraken = BACT_FILT + CONT_FILT
false_negative_kraken = ENDO_FILT

if str(species) == HUMAN:
    true_negative_kraken = BACT_FILT
    false_negative_kraken = ENDO_FILT + CONT_FILT
    tag_list = ENDO_CONT


df = pd.read_csv(input_file, dtype=str)
df[MAPPING_QUALITY] = pd.to_numeric(df[MAPPING_QUALITY])
mapping = {'True': True, 'False': False}
df[IS_UNMAPPED] = df[IS_UNMAPPED].map(mapping)

value_list = display_all_values(df, tag_list, true_negative_kraken=true_negative_kraken, false_negative_kraken=false_negative_kraken)
as_string = ','.join(str(item) for item in value_list)
print(name, "bwa aln", "mapped or unmapped", as_string, sep=',')

for mapq in MAPQS:
    value_list = display_all_values(df, tag_list, mode=1, value=mapq, true_negative_kraken=true_negative_kraken, false_negative_kraken=false_negative_kraken)
    as_string = ','.join(str(item) for item in value_list)
    comment = "MapQ > " + str(mapq)
    print(name, "bwa aln", comment, as_string, sep=',')









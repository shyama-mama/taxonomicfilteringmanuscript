import argparse
import textwrap
import sys
import pandas as pd

# CLASSIFIED, FASTQ_HEADER, PREDICTED_TAXID, SEQ_LEN, LIST
PREDICTED_TAXID = 'predicted_taxid'
CLASSIFIED = 'classified'
FASTQ_HEADER = 'fastqheader'
SEQ_LEN = 'seqlen'
LIST = 'list'
ZERO = '0'
PRIMARY_PCT_KMERS_USED = 'primary_pct_kmers_used'


HUMAN = '9606'
BISON = '43346'
DOG = '9615'

# Classification
# Homo sapiens (species) - 9606
# Homo (genus) - 9605
# Homininae (sub-family) - 207598
# Hominidae (family) - 9604 - Up to here in kraken2
# Simiiformes (infraorder) - 314293
# Primates (order) - 9443
# Mammalia (class) - 40674
# Chordata (phylum) - 7711
# Animalia (kingdom) - 33208
human_order = ['9606', '9605', '207598', '9604', '314295', '9526', '314293', '9443', '9598', '9595', '9601']

# Classification
# Bison bison bison (sub-species) - 43346
# Bison bison (species) - 9901
# Bison (genus) - 9900
# Bovinae (sub-family) - 27592
# Bovidae (family) - 9895 - Family
# Ruminantia (suborder) - 9845
# Artiodactyla (order) - 91561
# Mammalia (class) - 40674
# Chordata (phylum) - 7711
# Animalia (kingdom) - 33208
bison_order = ['43346', '9900', '27592', '9895']

# Classification 
# Canis lupus familiaris (sub-species) - 9615
# Canis lupus (species) - 9612
# Canis (genus) - 9611
# Canidae (family) - 9608
# dog_order = ['9615', '9612', '9611', '9608', '0']
dog_order = ['9615', '9612', '9611', '9608', '379584', '33554', '314145', '1437010', '9347', '32525', '40674', '32524', '7711', '33208', '286419']


# Make sense of and calculate some values from kraken2 k-mer breakdown 
def decode_kraken2_list(kraken2_breakdown, seqlen, taxid, classification):
    # Check if classified
    classified = True
    whacky = False
    if classification == "U":
        classified = False

    # Get kraken's primary tax prediction
    primary_taxid = str(taxid)
    breakdown_list = kraken2_breakdown.split(" ")
    seqlen = int(seqlen)

    cumulativeDict = {}
    total_k_mers_used = 0

    # for each element in breakdown list get the taxid and number of kmers
    # Store the cumulative total of kmers belonging to each tax id. 
    for element in breakdown_list:
        elem_list = element.split(':')
        key = str(elem_list[0])
        k_mers = int(elem_list[1])
        total_k_mers_used = total_k_mers_used + k_mers
        if key in cumulativeDict:
            cumulativeDict[key] = cumulativeDict[key] + k_mers
        else:
            cumulativeDict[key] = k_mers
    
    # In case we have a werid read
    if seqlen == 0 or total_k_mers_used == 0:
        # return 'NA', 0, 0, 0, 0, 0, 0, 0, True
        return 0.0

    if classified:
        # Remove primary taxID 
        # k_mers_used_for_predicted_taxid = cumulativeDict.pop(primary_taxid, None)
        # Changing from None to 0 (default value) because of this particular line
        # C	M_NC_000014.9:-:60570660:60570779:119e-2	131567	119	0:21 1185412:2 0:49 4182:2 0:11
        # Kraken has predicted it as 131567 but there is not 131567 in the breakdown list. 
        # Taxid is 131567 which is cellular organisms (no rank), biota (synonym).
        # 1185412: Defluviitaleaceae bacterium Ra1766G1 (species).
        # The other is seseme seeds 
        # 4182: Sesamum indicum (species), Sesamum indicum L. (authority), Sesamum orientale (synonym), Sesamum orientale L. (authority), beniseed (common name), gingelly (common name), hu ma (common name), koba (common name), sesame (genbank common name).
        
        # Adding this to catch out more of these whacky examples. Not sure how to deal with these for now just classifying them as whacky. 
        if not primary_taxid in cumulativeDict:
            whacky = True
        
        k_mers_used_for_predicted_taxid = cumulativeDict.pop(primary_taxid, 0)    
    else:
        k_mers_used_for_predicted_taxid = 0

    # remove 0 classification as it is just unclassified 
    k_mers_with_zero = cumulativeDict.pop(ZERO, 0) 
    

    # check if empty dict and find the second closest classification
    if cumulativeDict:
        k_mers_used_for_second_taxid = int(sorted(set(cumulativeDict.values()), reverse=True)[0])
        second_taxid = str(sorted(cumulativeDict, key=cumulativeDict.get, reverse=True)[0])
    else:
        k_mers_used_for_second_taxid = 0
        second_taxid = str('NA')

    # Calculating some values
    zero_pct_kmers_used = float(k_mers_with_zero/total_k_mers_used)
    zero_pct_total_seqlen = float(k_mers_with_zero/seqlen)

    primary_pct_kmers_used = float(k_mers_used_for_predicted_taxid/total_k_mers_used)
    primary_pct_total_seqlen = float(k_mers_used_for_predicted_taxid/seqlen)

    secondary_pct_kmers_used = float(k_mers_used_for_second_taxid/total_k_mers_used)
    secondary_pct_total_seqlen = float(k_mers_used_for_second_taxid/seqlen)

    pct_seq_used = float(total_k_mers_used/seqlen)
    
    # return second_taxid, zero_pct_kmers_used, zero_pct_total_seqlen, primary_pct_kmers_used, primary_pct_total_seqlen, secondary_pct_kmers_used, secondary_pct_total_seqlen, pct_seq_used, whacky
    return float(primary_pct_kmers_used)


# Arguements
parser = argparse.ArgumentParser(prog='kraken_extract.py',
   usage='python %(prog)s [-h] -i <input file> -s <taxid of endo sequences>',
   formatter_class=argparse.RawDescriptionHelpFormatter,
   description=textwrap.dedent('''\
   author:
     Shyamsundar Ravishankar (shyamsundar.ravishankar@adelaide.edu.au)

   description:
     %(prog)s Script that gets fastq headers of species of interest from kraken output, also keeps track of what reads were removed
   '''))
parser.add_argument("-i", "--input_file", dest="input_file", default=None, help="kraken2 output file", required=True)
parser.add_argument("-o", "--output_file", dest="output_file", default="out.csv", help="Path to output file")
parser.add_argument("-s", "--species", dest="species", default=None, help="Species to positively select", required=True)
parser.add_argument("-t", "--stats_file", dest="stats_file", default="filtered_stats.csv", help="Summary of what was filtered out")
parser.add_argument("-f", "--filter_value", dest="filter_value", default=0.0, help="pct_kmers_used threshold to remove contaminants not in species")
parser.add_argument("-e", "--emperical", dest="emperical", action="store_true", help="If this is an emperical sample.", required=False)
args = parser.parse_args()

output_file = args.output_file
stats_file = args.stats_file
filter_value = float(args.filter_value)
emperical = args.emperical

input_file = args.input_file
input_df = pd.read_csv(input_file, header=None, sep='\t', dtype=str)
input_df.columns = [CLASSIFIED, FASTQ_HEADER, PREDICTED_TAXID, SEQ_LEN, LIST]

print(len(input_df.index))

# calculate PRIMARY_PCT_KMERS_USED
#input_df[PRIMARY_PCT_KMERS_USED] = input_df.apply(lambda row : decode_kraken2_list(row[LIST], row[SEQ_LEN], row[PREDICTED_TAXID], row[CLASSIFIED]), axis = 1, result_type='expand')

total_size = len(input_df.index)

# Drop columns we don't need
input_df.drop(CLASSIFIED, inplace=True, axis=1)
input_df.drop(SEQ_LEN, inplace=True, axis=1)
input_df.drop(LIST, inplace=True, axis=1)


species = str(args.species)

order = ['0']

if (species == HUMAN):
  order = human_order
elif (species == BISON):
  order = bison_order
elif (species == DOG):
  order = dog_order


# Extract the desired rows from df that is in our order 
# fastq_headers = input_df[input_df[PREDICTED_TAXID].isin(order)][FASTQ_HEADER]

#df_indexes = input_df.index[(~input_df[PREDICTED_TAXID].isin(order)) & (input_df[PRIMARY_PCT_KMERS_USED] >= filter_value)]
df_indexes = input_df.index[(~input_df[PREDICTED_TAXID].isin(order))]

# get rows from df not in df1
fastq_headers = input_df.drop(df_indexes)

# Filtered rows 
filtered_rows = input_df.loc[df_indexes]

size = len(fastq_headers)
fil_size = len(filtered_rows)
total = size + fil_size
print("Found " + str(size) + " reads belonging to " + str(species) + " filtered " + str(fil_size) + " reads.")

if (size == 0):
  sys.exit(1)


sf = open(stats_file, "w")
sf.write("bact,cont,endo\n")
if(emperical):
  sf.write("0,0,0\n")
else:
  filtered_cont = filtered_rows[FASTQ_HEADER].str.split(":").str[5].str.count('c').sum()
  filtered_bact = filtered_rows[FASTQ_HEADER].str.split(":").str[5].str.count('b').sum()
  filtered_endo = filtered_rows[FASTQ_HEADER].str.split(":").str[5].str.count('e').sum()
  sf.write(str(filtered_bact) + "," + str(filtered_cont) + "," + str(filtered_endo)  + "\n")
sf.close()

f = open(output_file, "w")
for header in fastq_headers[FASTQ_HEADER]:
  f.write(header.strip() + '\n')
f.close()

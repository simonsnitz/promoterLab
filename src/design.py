
from pprint import pprint
import json

import streamlit as st
from promoter_calculator import promoter_calculator



# Output file name
out_file_name = "mepr_slips.json"


@st.cache_data
def create_library(operator, base_promoter):
# Create promoter library
    operator_length = len(operator)
    init_pos = 67       # this is the position in base_promoter where the start of the -35 region is found
    insert_pos = 67 - operator_length
    end_pos = len(base_promoter)
    total_promoters = end_pos - operator_length +1
    promoters = []
    for pos in range(insert_pos, total_promoters):
        start_prom = base_promoter[0:pos]
        end_prom = base_promoter[(pos+operator_length):end_pos]
        promoter_seq = start_prom.lower() + operator + end_prom.lower()
        promoters.append({"seq":promoter_seq})
    
    return promoters

@st.cache_data
def calculate_overlap(base_promoter, promoters):
    overlap_pos = [i for i in range(67,101)]
    scoring_counter = []
    for i in range(0, len(base_promoter)):
        if i in overlap_pos:
            score = 1
        else:
            score = 0
        scoring_counter.append(score)


    for i in range(0,len(promoters)):
        prom = promoters[i]["seq"]
        score = 0
        for j in range(0, len(prom)):
            if prom[j].isupper():

                score += scoring_counter[j]

        promoters[i]["score"] = score

    return promoters


@st.cache_data
def calculate_tx_rates(promoters):

    for i in range(0, len(promoters)):
        seq = promoters[i]["seq"]
        results = promoter_calculator(sequence=seq, threads=4)
        # Find the best promoter
        max_promoter = 0
        # prom_deets = 0
        for promoter in results:
            if promoter.strand == "+":
                if promoter.Tx_rate > max_promoter:
                    max_promoter = promoter.Tx_rate
                    prom_deets = promoter
        
        promoters[i]["tx_rate"] = max_promoter

    return promoters





if __name__ == "__main__":

    # Define inputs
    base_promoter = "aaaaatggcgcccatcggcgccatttttttatggccatgtattaaaatatatttttcaaaagtatcgTTGACGgcgtatctcttgctttcTATAATgctatcgatcgatcgtctaattgagctgtcaccggatgtgctt"
    # RamR operator
    # operator = "ATAATGAGTGCTTACTCACTCATAAT"
    # MepR operator
    operator = "ATTTGGTTAGACATCTAACGAAAT"

    promoters = create_library(operator, base_promoter)
    promoters = calculate_overlap(base_promoter, promoters)
    #print(promoters)


# truncate promoter list for development
# promoters = promoters[0:5]


'''

# Create an overlap scoring function
pos_35_init = [67, 68, 69]
pos_35_term = [70, 71, 72]
pos_10ext = [88, 89]
pos_10_init = [90, 91, 92]
pos_10_term = [93, 94, 95]
pos_disc = [96, 97, 98]

scoring_counter = []
for i in range(0, len(base_promoter)):
    if i in pos_35_init:
        score = 7.5
    elif i in pos_35_term:
        score = 5.5
    elif i in pos_10ext:
        score = 3
    elif i in pos_10_init:
        score = 8.5
    elif i in pos_10_term:
        score = 8
    elif i in pos_disc:
        score = 1.5
    else:
        score = 0
    scoring_counter.append(score)


for i in range(0,len(promoters)):
    prom = promoters[i]["seq"]
    score = 0
    for j in range(0, len(prom)):
        if prom[j].isupper():
            score += scoring_counter[j]

    promoters[i]["score"] = score


# Calculate transcription rates
for i in range(0, len(promoters)):
    seq = promoters[i]["seq"]
    results = promoter_calculator(seq)
    # Find the best promoter
    max_promoter = 0
    prom_deets = 0
    for promoter in results:
        if promoter.strand == "+":
            if promoter.Tx_rate > max_promoter:
                max_promoter = promoter.Tx_rate
                prom_deets = promoter
    
    promoters[i]["tx_rate"] = max_promoter

with open(out_file_name, "w+") as f:
    json.dump(promoters, f)

# pprint(promoters)

'''
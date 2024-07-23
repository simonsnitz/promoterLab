from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from tempfile import NamedTemporaryFile
import streamlit.components.v1 as components


import math
import base64
import streamlit as st




def create_genbank(part_name, promoter, op_positions, overlap_score):

    part_primer_F = "GCACTGCTCACTGCATCGATACATG"
    part_primter_R = "GGATGAGTCAGCATCGAGCACAGAC"
    upstream_spacer = "aggactacgtggcgagattagatgcagctacaagacacgctggcgtgaacaagaggtgagctcaattcacattggttgatacgcgcaatcttggattgcgcgaa"
    downstream_spacer = "tacatggcgtggtctaggatctgacactcgcacgatagactggactgcgatcctgtgctgagacatcaatgatcactccaaggttatcttagtcgatc"
    upstream_bsaI_cut = "ggtctcaGCCA"
    downstream_bsaI_cut = "ATTGtgagacc"

    full_promoter = promoter.promoter_sequence
    promoter_bsaI_primer_len = len(full_promoter) + 22 + 50
    spacer_append_len = math.ceil((300 - promoter_bsaI_primer_len +1)/2) *2

    upstream_spacer = upstream_spacer[0:int(spacer_append_len/2)]
    downstream_spacer = downstream_spacer[0:int(spacer_append_len/2)]

    UP_start = 25 + len(upstream_spacer) + 11
    UP_end = UP_start + promoter.UP_position[1] - promoter.UP_position[0]

    operator_len = op_positions[1] - op_positions[0]
    operator_start = UP_start + (op_positions[0] - promoter.UP_position[0])
    operator_end = operator_start + operator_len    

    disc_len = promoter.disc_position[1] - promoter.disc_position[0]
    disc_start = UP_start + (promoter.disc_position[0] - promoter.UP_position[0])
    disc_end = disc_start + disc_len

    promSpacer_len = promoter.spacer_position[1] - promoter.spacer_position[0]
    promSpacer_start = UP_start + (promoter.spacer_position[0] - promoter.UP_position[0])
    promSpacer_end = promSpacer_start + promSpacer_len

    hex10_len = promoter.hex10_position[1] - promoter.hex10_position[0]
    hex10_start = UP_start + (promoter.hex10_position[0] - promoter.UP_position[0])
    hex10_end = hex10_start + hex10_len

    hex35_len = promoter.hex35_position[1] - promoter.hex35_position[0]
    hex35_start = UP_start + (promoter.hex35_position[0] - promoter.UP_position[0])
    hex35_end = hex35_start + hex35_len

    annotations = [
    {"type": "misc_feature",
    "label": "Part-Primer_F",
    "color": "#11ff00",
    "start": 0,
    "end": 25,
    "strand": 1},

    {"type": "misc_feature",
    "label": "spacer",
    "color": "#bababa",
    "start": 25,
    "end": 25 + len(upstream_spacer),
    "strand": 1},

    {"type": "misc_feature",
    "label": "BsaI",
    "color": "#ff70e7",
    "start": 25 + len(upstream_spacer),
    "end": 25 + len(upstream_spacer) + 6,
    "strand": 1},

    {"type": "misc_feature",
    "label": "Cut",
    "color": "#24edff",
    "start": 25 + len(upstream_spacer) + 7,
    "end": 25 + len(upstream_spacer) + 11,
    "strand": 1},

    {"type": "misc_feature",
    "label": str(part_name)+ " Promoter",
    "color": "#ffae00",
    "start": 25 + len(upstream_spacer) + 11,
    "end": 25 + len(upstream_spacer) + 11 + len(full_promoter),
    "strand": 1},

    {"type": "misc_feature",
    "label": "UP element",
    "color": "#ffff00",
    "start": UP_start,
    "end": UP_end,
    "strand": 1},

    {"type": "misc_feature",
    "label": "Operator",
    "color": "#636363",
    "start": operator_start,
    "end": operator_end,
    "strand": 1},

    {"type": "misc_feature",
    "label": "spacer",
    "color": "#ffff00",
    "start": promSpacer_start,
    "end": promSpacer_end,
    "strand": 1},

    {"type": "misc_feature",
    "label": "disc",
    "color": "#ffff00",
    "start": disc_start,
    "end": disc_end,
    "strand": 1},

    {"type": "misc_feature",
    "label": "ITR",
    "color": "#ffff00",
    "start": disc_end,
    "end": disc_end + 20,
    "strand": 1},

    {"type": "misc_feature",
    "label": "-10",
    "color": "#ff1e00",
    "start": hex10_start,
    "end": hex10_end,
    "strand": 1},

    {"type": "misc_feature",
    "label": "-35",
    "color": "#ff1e00",
    "start": hex35_start,
    "end": hex35_end,
    "strand": 1},

    {"type": "misc_feature",
    "label": "Cut",
    "color": "#24edff",
    "start": 25 + len(upstream_spacer) + 11 + len(full_promoter),
    "end": 25 + len(upstream_spacer) + 11 + len(full_promoter) + 4,
    "strand": -1},

    {"type": "misc_feature",
    "label": "BsaI",
    "color": "#ff70e7",
    "start": 25 + len(upstream_spacer) + 11 + len(full_promoter) + 5,
    "end": 25 + len(upstream_spacer) + 11 + len(full_promoter) + 11,
    "strand": -1},

    {"type": "misc_feature",
    "label": "spacer",
    "color": "#bababa",
    "start": 25 + len(upstream_spacer) + 11 + len(full_promoter) + 11,
    "end": 25 + len(upstream_spacer) + 11 + len(full_promoter) + 11 + len(downstream_spacer),
    "strand": -1},

    {"type": "misc_feature",
    "label": "Part-Primer_R",
    "color": "#11ff00",
    "start": 25 + len(upstream_spacer) + 11 + len(full_promoter) + 11 + len(downstream_spacer),
    "end": 25 + len(upstream_spacer) + 11 + len(full_promoter) + 11 + len(downstream_spacer) + 25,
    "strand": -1},
    ]

    # Create the plasmid sequence record

    seq = part_primer_F + upstream_spacer + upstream_bsaI_cut + promoter.promoter_sequence + downstream_bsaI_cut + downstream_spacer + part_primter_R
    
    plasmid_sequence = Seq(seq)
    record = SeqRecord(plasmid_sequence,
                    id=str(part_name),
                    name=part_name,
                    description='Designed by promoterLab_V1.\
                                Transcription rate: '+str(round(promoter.Tx_rate,1))+". \
                                Overlap score: "+str(overlap_score),
                    annotations={"molecule_type": "DNA"})

    # Set topology as circular
    record.annotations["topology"] = "linear"



    for annotation in annotations:
        record.features.append(\
            SeqFeature(FeatureLocation(\
                start = annotation["start"], \
                end = annotation["end"], \
                strand = annotation["strand"]), \
                type = annotation["type"], \
                qualifiers={ \
                    "ApEinfo_fwdcolor":[annotation["color"]],\
                    "label": annotation["label"]}))
        


    # Converts gbk into straight text
    def get_gbk(record):
        outfileloc=NamedTemporaryFile()
        with open(outfileloc.name, "w") as handle:
            SeqIO.write(record, handle, "genbank")
        with open(outfileloc.name) as handle:
            record=handle.read()
        outfileloc.close()

        return record

    record = get_gbk(record)
    
    return record





def download_button(object_to_download, download_filename):
    # Generates a link to download the given object_to_download.
    try:
        # some strings <-> bytes conversions necessary here
        b64 = base64.b64encode(object_to_download.encode()).decode()

    except AttributeError as e:
        b64 = base64.b64encode(object_to_download).decode()

    dl_link = f"""
    <html>
    <head>
    <title>Start Auto Download file</title>
    <script src="http://code.jquery.com/jquery-3.2.1.min.js"></script>
    <script>
    $('<a href="data:text/csv;base64,{b64}" download="{download_filename}">')[0].click()
    </script>
    </head>
    </html>
    """
    return dl_link


def download_promoter_df(part_name, promoter, op_positions, overlap_score):
    

    data = create_genbank(part_name, promoter, op_positions, overlap_score)
    '''
    Having trouble downloading the GenBank file when not locally hosted, so we'll just display it in the browser.
    '''
    components.html(
        download_button(data, str(part_name)+"_promoter-part.gb"),
        height=0,
    )


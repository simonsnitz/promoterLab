from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from tempfile import NamedTemporaryFile
import streamlit.components.v1 as components

import requests
from dnachisel import *
import base64
import streamlit as st




    # Input protein accession ID, output sequence in fasta format
def accID2seq(accID: str):
    URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi/?db=protein&id="+accID+"&rettype=fasta"
    response = requests.get(URL)
    if response.ok:
        fasta = response.text.split("\n")
        fasta = [i for i in fasta if len(i) != 0]
        fasta = "".join(i for i in fasta if i[0] != ">")
        return fasta
    else:
        print("FATAL: Bad eFetch request for accID2seqeuence"+ str(response.status_code))
        return None
    



def codon_opt(protein_seq:str):

    if protein_seq != None:
        # Create a random DNA seq given the protein seq. Append a stop codon.
        protein_dna_seq = reverse_translate(protein_seq+"*")

        # DEFINE THE OPTIMIZATION PROBLEM
        problem = DnaOptimizationProblem(
            sequence=protein_dna_seq,
            constraints=[
                AvoidPattern("BsaI_site"),
                EnforceGCContent(mini=0.35, maxi=0.65, window=50),
                EnforceTranslation(location=(0, len(protein_dna_seq)))
            ],
            objectives=[CodonOptimize(species='e_coli', location=(0, len(protein_dna_seq)))]
        )

        # SOLVE THE CONSTRAINTS, OPTIMIZE WITH RESPECT TO THE OBJECTIVE

        problem.resolve_constraints()
        problem.optimize()

        # GET THE FINAL SEQUENCE (AS STRING OR ANNOTATED BIOPYTHON RECORDS)

        final_sequence = problem.sequence  # string
        return final_sequence








def create_genbank(part_name, CDS, protein_name):

    upstream_CDS = "GCACTGCTCACTGCATCGATACATGacggtctcaT"
    downstream_CDS = "taATCCtgagacctgGGATGAGTCAGCATCGAGCACAGAC"

    annotations = [
    {"type": "misc_feature",
    "label": "Part-Primer_F",
    "color": "#00ffea",
    "start": 0,
    "end": 25,
    "strand": 1},

    {"type": "misc_feature",
    "label": "BsaI",
    "color": "#ff70e7",
    "start": 27,
    "end": 33,
    "strand": 1},

    {"type": "misc_feature",
    "label": str(protein_name),
    "color": "#cfa8ff",
    "start": 35,
    "end": 35+len(CDS),
    "strand": 1},

    {"type": "CDS",
    "label": str(protein_name),
    "color": "#cfa8ff",
    "start": 35,
    "end": 35+len(CDS),
    "strand": 1,
    "translation": CDS},

    {"type": "misc_feature",
    "label": "Cut",
    "color": "#70ff9d",
    "start": 35+len(CDS)+2,
    "end": 35+len(CDS)+6,
    "strand": -1},

    {"type": "misc_feature",
    "label": "BsaI",
    "color": "#ff70e7",
    "start": 35+len(CDS)+7,
    "end": 35+len(CDS)+13,
    "strand": -1},

    {"type": "misc_feature",
    "label": "Part-Primer-R",
    "color": "#00ffea",
    "start": 35+len(CDS)+15,
    "end": 35+len(CDS)+40,
    "strand": -1},
    ]

    # Create the plasmid sequence record

    seq = upstream_CDS + CDS + downstream_CDS
    
    plasmid_sequence = Seq(seq)
    record = SeqRecord(plasmid_sequence,
                    id=str(part_name),
                    name=str(protein_name)+"_CDS",
                    description='This is a genetic assembly part designed by promoterLab',
                    annotations={"molecule_type": "DNA"})

    # Set topology as circular
    record.annotations["topology"] = "linear"



    for annotation in annotations:
        if "translation" in annotation.keys():
            record.features.append(\
                SeqFeature(FeatureLocation(\
                    start = annotation["start"], \
                    end = annotation["end"], \
                    strand = annotation["strand"]), \
                    type = annotation["type"], \
                    qualifiers={ \
                        "ApEinfo_fwdcolor":[annotation["color"]],\
                        "translation":[annotation["translation"]],\
                        "label": annotation["label"]}))
        else:
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


def download_cds_df(part_name, CDS, protein_name):
    

    data = create_genbank(part_name, CDS, protein_name)
    '''
    Having trouble downloading the GenBank file when not locally hosted, so we'll just display it in the browser.
    '''
    components.html(
        download_button(data, str(protein_name)+"_CDS-part.gb"),
        height=0,
    )




if __name__ == "__main__":
    codon_opt("MTTIRWRRMSIHSERITLADSPLHWAHTLNGSMRTHFEVQRLERGRGAYLARSRFGAGELYSAIAPSQVLRHFNDQRNANEAEHSYLIQIRSGALGVASGGRKVILANGDCSIVDSRQDFTLSSNSSTQGVVIRFPVSWLGAWVSNPEDLIARRVDAEIGWGRALSASVSNLDPLRIDDLGSNVNSIAEHVAMLISLASSAVSSEDGGVALRKMREVKRVLEQSFADANLEPESVSSQLGISKRYLHYVFAACGTTFGRELLEIRLGKAYRMLCATSGSGAVLKVAMSSGFSDSSHFSKKFKERYGVSPVSLVRQA")

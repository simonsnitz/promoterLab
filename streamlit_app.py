import streamlit as st
import pandas as pd
import json

from src.design import create_library
from src.design import calculate_overlap
from src.design import calculate_tx_rates
from src.genbank.CDS_utils import accID2seq, codon_opt, download_cds_df
from src.genbank.promoter_utils import download_promoter_df

from promoter_calculator import promoter_calculator


st.set_page_config(page_title="promoterCAD", layout='wide', initial_sidebar_state='auto')


hide_streamlit_style = '''
<style>
#MainMenu {visibility: hidden;}
footer {visibility: hidden;}
</style>
'''
st.markdown(hide_streamlit_style, unsafe_allow_html=True)

# Removes border around forms
css = r'''
    <style>
        [data-testid="stForm"] {border: 0px}
    </style>
'''
st.markdown(css, unsafe_allow_html=True)

# Removes the full-screen button for various elements
style_fullscreen_button_css = """
    button[title="View fullscreen"] {
        display: none;
    }
    button[title="View fullscreen"]:hover {
        display: none;
        }
    """
st.markdown(
    "<style>"
    + style_fullscreen_button_css
    + "</styles>",
    unsafe_allow_html=True,
)

# Initialize state variables
if "data" not in st.session_state:
        st.session_state.data = False

if 'SUBMITTED' not in st.session_state:
    st.session_state.SUBMITTED =  False


def _connect_form_cb(connect_status):
    st.session_state.SUBMITTED = connect_status
    st.session_state.data = False




# HEADER: Title and basic input
head = st.container()
head1, head2, head3 = head.columns((1,1,1))

head2.image("images/promoterLab_logo.png", use_column_width=True)
head2.markdown("<h3 style='text-align: center; color: black;'>Design ligand-inducible promoters for prokaryotes</h3>", unsafe_allow_html=True)

selection_container = st.container()
sel1, sel2, sel3 = selection_container.columns((1,2,1))

promoter = "aaaaatggcgcccatcggcgccatttttttatggccatgtattaaaatatatttttcaaaagtatcgTTGACGgcgtatctcttgctttcTATAATgctatcgatcgatcgtctaattgagctgtcaccggatgtgctt"


protein = sel2.text_input(label="Protein ID", value="WP_013083972.1")
base_promoter = sel2.text_input(label="Base promoter", value=promoter)
operator = sel2.text_input(label="Operator", value="ATTTGGTTAGACATCTAACGAAAT")
promoter_name = sel2.text_input(label="Promoter_name", value="Pmttr")

operation = sel2.radio("operation", ["Create CDS assembly part", "Predict Tx rate", "Design promoters"])



def predict_tx_rate(seq):
    results = promoter_calculator(sequence=seq, threads=4)
    # Find the best promoter
    max_promoter = 0
    prom_deets = 0
    for promoter in results:
        if promoter.strand == "+":
            if promoter.Tx_rate > max_promoter:
                max_promoter = promoter.Tx_rate
                prom_deets = promoter

    return prom_deets




# FORM
with st.form(key='promoterLab'):

    # SUBMIT BUTTON
    submit = st.container()
    submit_spacer_1, submit_button, submit_spacer_2 = submit.columns([5,1,5])
    submitted = submit_button.form_submit_button("Submit", use_container_width=True, on_click=_connect_form_cb, args=(True,))


    # RUN PROMOTERLAB
    if st.session_state.SUBMITTED:

        if operation == "Predict Tx rate":

            prom_deets = predict_tx_rate(base_promoter)
            max_promoter = round(prom_deets.Tx_rate, 2)

            st.subheader("Transcription rate: "+str(max_promoter))
            st.markdown("Promoter sequence: "+str(prom_deets.promoter_sequence))
            st.text("UP element: "+str(prom_deets.UP))
            st.text("-35: "+str(prom_deets.hex35))
            st.text("Spacer: "+str(prom_deets.spacer))
            st.text("-10: "+str(prom_deets.hex10))
            st.text("disc: "+str(prom_deets.disc))
            st.text("ITR: "+str(prom_deets.ITR))



        elif operation == "Create CDS assembly part":

            protein_fasta = accID2seq(protein)

            CDS = codon_opt(protein_fasta)

            output = st.container()
            spacer1, out, spacer2 = output.columns([1,3,1])   

            out.form_submit_button(label="Download CDS Part", type="primary", on_click=download_cds_df, args=(protein, CDS))




        elif operation == "Design promoters":


            output = st.container()
            spacer1, outputTable, spacer2 = output.columns([1,3,1])


            promoters = create_library(operator, base_promoter)
            promoters = calculate_overlap(base_promoter, promoters)
            with st.spinner("calculating transcription rates"):
                promoters = calculate_tx_rates(promoters)

            # For testing purposes:
            # with open("promoters.json", "r") as f:
            #     promoters = json.load(f)


            # select the best promoter
            strong = [p for p in promoters if p["tx_rate"] > 15000]
            highest_overlap = max([s["score"] for s in strong])
            best = [b for b in strong if b["score"] == highest_overlap][0]
            overlap_score = best["score"]

            # get position of operator
            op_positions = []
            for i in range(0,len(best["seq"])):
                if best["seq"][i].isupper() == True:
                    op_positions.append(i)
            op_positions = [min(op_positions),max(op_positions)]


            prom_deets = predict_tx_rate(best['seq'])

            output = st.container()
            spacer1, out, spacer2 = output.columns([1,3,1])  

            out.form_submit_button(label="Download Promoter Part", type="primary", on_click=download_promoter_df, args=(promoter_name, prom_deets, op_positions, overlap_score))


            # Plot

            # x = [i for i in range(0, len(promoters))]
            # promoter_seqs = [i["seq"] for i in promoters]
            # overlap_scores = [i["score"] for i in promoters]
            # tx_rates = [i["tx_rate"] for i in promoters]
            # for i in range(0, len(promoters)):
            #     promoters[i]["index"] = i

            # outputTable.dataframe(data=promoters, height=350, width=1000)

            # chart_data = pd.DataFrame(promoters, columns=["seq", "score", "tx_rate", "index"])

            # st.scatter_chart(data=chart_data, x="tx_rate", y="score", color="index", height=600)
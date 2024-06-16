import streamlit as st
import requests
import json
import plotly.figure_factory as ff
import pandas as pd

from src.design import create_library
from src.design import calculate_overlap
from src.design import calculate_tx_rates


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

# head2.image("images/Snowprint_Logo.png", use_column_width=True)
head2.image("images/promoterLab_logo.png", use_column_width=True)
head2.markdown("<h3 style='text-align: center; color: black;'>Design ligand-inducible promoters for prokaryotes</h3>", unsafe_allow_html=True)

selection_container = st.container()
sel1, sel2, sel3 = selection_container.columns((1,2,1))

promoter = "aaaaatggcgcccatcggcgccatttttttatggccatgtattaaaatatatttttcaaaagtatcgTTGACGgcgtatctcttgctttcTATAATgctatcgatcgatcgtctaattgagctgtcaccggatgtgctt"
base_promoter = sel2.text_input(label="Base promoter", value=promoter)

operator = sel2.text_input(label="Operator", value="ATTTGGTTAGACATCTAACGAAAT")


# FORM
with st.form(key='promoterLab'):

    # SUBMIT BUTTON
    submit = st.container()
    submit_spacer_1, submit_button, submit_spacer_2 = submit.columns([5,1,5])
    submitted = submit_button.form_submit_button("Submit", use_container_width=True, on_click=_connect_form_cb, args=(True,))

# RUN PROMOTERLAB
if st.session_state.SUBMITTED:

    output = st.container()
    spacer1, outputTable, spacer2 = output.columns([1,3,1])

    with st.spinner("creating promoter library"):

        promoters = create_library(operator, base_promoter)


    with st.spinner("calculating overlap scores"):

        promoters = calculate_overlap(base_promoter, promoters)


    with st.spinner("calculating transcription rates"):

        promoters = calculate_tx_rates(promoters)



    # Plot

    x = [i for i in range(0, len(promoters))]
    promoter_seqs = [i["seq"] for i in promoters]
    overlap_scores = [i["score"] for i in promoters]
    tx_rates = [i["tx_rate"] for i in promoters]
    for i in range(0, len(promoters)):
        promoters[i]["index"] = i

    outputTable.dataframe(data=promoters, height=350, width=1000)

    chart_data = pd.DataFrame(promoters, columns=["seq", "score", "tx_rate", "index"])

    st.scatter_chart(data=chart_data, x="tx_rate", y="score", color="index", height=600)
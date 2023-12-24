import bioinformatics
import streamlit as st

st.title("Bioinformatics Application - BFCAI\n")
st.header("An introduction\n")
st.markdown("**This application for educational purposes only**")
st.markdown("""**The biology field is intersected with many of `AI` applications, which is important to use to identify
the sequences of DNA, the translation into protein sequence and searching for `PATTERNS` !**
""")
st.divider()
st.header("Searching")
st.markdown("""**This section is for searching for patterns inside a sequence**\n
- Upload your `.fasta` file that contains the sequence""")
file = st.file_uploader("Fasta file", type='fasta')
with open(file, 'r') as fasta:
    text = fasta.read()
    text = text.split('\n')[1]
    
st.write(text)
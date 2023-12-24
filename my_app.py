from bioinformatics import *
import streamlit as st
cnt = 0
st.title("Bioinformatics Application - BFCAI\n")
st.header("An introduction\n")
st.markdown("**This application for educational purposes only**")
st.markdown("""**The biology field is intersected with many of `AI` applications, which is important to use to identify
the sequences of DNA, the translation into protein sequence and searching for `PATTERNS` !**
""")
st.divider()
tab1, tab2 = st.tabs(['DNA Translation', 'DNA Sequence Searching'])
with tab1:
    file = st.file_uploader("Fasta file", type=['fasta','fa'], key=0)
    text = None
    if file is not None:
        text = file.getvalue().decode(encoding='utf-8')
        text = text.split('\n')[1]

    st.subheader("This is your DNA sequence:")
    if text is not None:
        st.write(text)
        st.subheader("The translated protein:")
        st.markdown(f'**:green[{translation(text)}]**')
    else:
        st.write(":red[No DNA sequence is added yet !]")

with tab2:
    file = st.file_uploader("Fasta file", type='fasta',key=1)
    text = None
    if file is not None:
        text = file.getvalue().decode(encoding='utf-8')
        text = text.split('\n')[1]

    st.subheader("This is your DNA sequence:")
    if text is not None:
        selected = -1
        st.write(text)
        pattern = st.text_input("**Pattern required**", placeholder="Enter your pattern to search in UPPERCASE")
        with st.container():
            options = ['Naive Matching',
            'Bad Character','Good Suffix','Boyer Moore',
            'KMP Searching','K-Meer Searching']
            selected = st.selectbox("Pick a search technique",sorted(options),index=None)
            if 'prev_t' not in st.session_state:
                st.session_state.prev_t = None

        if selected == 'Bad Character':
            if 'cnt' not in st.session_state:
                    st.session_state.cnt = 0
            
            if st.session_state.prev_t != pattern:
                st.session_state.cnt = 0
            st.session_state.prev_t = pattern
            info = boyer_moore_bad_char(text,pattern)
            found = info['found in']
            show_text = st.empty()

            if len(found) > 0:
                found_pattern = text[:found[st.session_state.cnt]] + ':blue[→]' + text[found[st.session_state.cnt]:]
                show_text.write(found_pattern)
                st.write(":red[**Skipped alignments:**] " + str(info['skipped alignments']))


                if st.button('Next location'):
                    st.session_state.cnt += 1
                    if st.session_state.cnt >= len(found):
                        st.session_state.cnt = 0
                    found_pattern = text[:found[st.session_state.cnt]] + ':blue[→]' + text[found[st.session_state.cnt]:]
                    show_text.write(found_pattern)
            else:
                show_text.write("The pattern is not found")
        elif selected == 'Good Suffix':
            if 'cnt' not in st.session_state:
                    st.session_state.cnt = 0
            
            if st.session_state.prev_t != pattern:
                st.session_state.cnt = 0
            st.session_state.prev_t = pattern
            info = boyer_moore_good_suffix(text,pattern)
            found = info['found in']
            show_text = st.empty()

            if len(found) > 0:
                found_pattern = text[:found[st.session_state.cnt]] + ':blue[→]' + text[found[st.session_state.cnt]:]
                show_text.write(found_pattern)
                st.write(":red[**Skipped alignments:**] " + str(info['skipped alignments']))


                if st.button('Next location'):
                    st.session_state.cnt += 1
                    if st.session_state.cnt >= len(found):
                        st.session_state.cnt = 0
                    found_pattern = text[:found[st.session_state.cnt]] + ':blue[→]' + text[found[st.session_state.cnt]:]
                    show_text.write(found_pattern)
            else:
                show_text.write("The pattern is not found")
        elif selected == 'Boyer Moore':
            if 'cnt' not in st.session_state:
                    st.session_state.cnt = 0
            
            if st.session_state.prev_t != pattern:
                st.session_state.cnt = 0
            st.session_state.prev_t = pattern
            info = boyer_moore(text,pattern)
            found = info['found in']
            show_text = st.empty()

            if len(found) > 0:
                found_pattern = text[:found[st.session_state.cnt]] + ':blue[→]' + text[found[st.session_state.cnt]:]
                show_text.write(found_pattern)
                st.write(":red[**Skipped alignments:**] " + str(info['skipped alignments']))


                if st.button('Next location'):
                    st.session_state.cnt += 1
                    if st.session_state.cnt >= len(found):
                        st.session_state.cnt = 0
                    found_pattern = text[:found[st.session_state.cnt]] + ':blue[→]' + text[found[st.session_state.cnt]:]
                    show_text.write(found_pattern)
            else:
                show_text.write("The pattern is not found")
        elif selected == 'K-Meer Searching':
            if 'cnt' not in st.session_state:
                    st.session_state.cnt = 0
            
            if st.session_state.prev_t != pattern:
                st.session_state.cnt = 0
            st.session_state.prev_t = pattern
            info = k_meer_search(text,pattern)
            found = info['found in']
            show_text = st.empty()

            if len(found) > 0:
                found_pattern = text[:found[st.session_state.cnt]] + ':blue[→]' + text[found[st.session_state.cnt]:]
                show_text.write(found_pattern)
                st.write(":red[**Skipped alignments:**] " + str(info['skipped alignments']))


                if st.button('Next location'):
                    st.session_state.cnt += 1
                    if st.session_state.cnt >= len(found):
                        st.session_state.cnt = 0
                    found_pattern = text[:found[st.session_state.cnt]] + ':blue[→]' + text[found[st.session_state.cnt]:]
                    show_text.write(found_pattern)
            else:
                show_text.write("The pattern is not found")
        elif selected == 'KMP Searching':
            if 'cnt' not in st.session_state:
                    st.session_state.cnt = 0
            
            if st.session_state.prev_t != pattern:
                st.session_state.cnt = 0
            st.session_state.prev_t = pattern
            found = KMP(text,pattern)
            show_text = st.empty()

            if len(found) > 0:
                found_pattern = text[:found[st.session_state.cnt]] + ':blue[→]' + text[found[st.session_state.cnt]:]
                show_text.write(found_pattern)

                if st.button('Next location'):
                    st.session_state.cnt += 1
                    if st.session_state.cnt >= len(found):
                        st.session_state.cnt = 0
                    found_pattern = text[:found[st.session_state.cnt]] + ':blue[→]' + text[found[st.session_state.cnt]:]
                    show_text.write(found_pattern)
            else:
                show_text.write("The pattern is not found")
        elif selected == 'Naive Matching':
            if 'cnt' not in st.session_state:
                    st.session_state.cnt = 0
            
            if st.session_state.prev_t != pattern:
                st.session_state.cnt = 0
            st.session_state.prev_t = pattern
            info = naive_matching(text,pattern)
            found = info['found in']
            show_text = st.empty()

            if len(found) > 0:
                found_pattern = text[:found[st.session_state.cnt]] + ':blue[→]' + text[found[st.session_state.cnt]:]
                show_text.write(found_pattern)
                st.write(":red[**Skipped alignments:**] " + str(info['skipped alignments']))


                if st.button('Next location'):
                    st.session_state.cnt += 1
                    if st.session_state.cnt >= len(found):
                        st.session_state.cnt = 0
                    found_pattern = text[:found[st.session_state.cnt]] + ':blue[→]' + text[found[st.session_state.cnt]:]
                    show_text.write(found_pattern)
            else:
                show_text.write("The pattern is not found")
                
    else:
        st.write("No DNA sequence is added yet !")


    


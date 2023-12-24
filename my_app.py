from bioinformatics import *
import streamlit as st
col1, col2, col3 = st.columns(3)
with col1:
    st.image('images/bionformatics', use_column_width=True)
with col2:
    st.image('images/bionformatics2', use_column_width=True)
with col3:
    st.image('images/bionformatics3', use_column_width=True)

st.title(":blue[Bioinformatics Application - BFCAI]")
st.header("Welcome To Our Website")
st.markdown("**This application for educational purposes only**")
st.markdown("""**The biology field is intersected with many of `AI` applications, which is important to use to identify
the sequences of DNA, the translation into protein sequence and searching for `PATTERNS` !**
""")
st.divider()
st.subheader(":blue[**Tools**]")
tab1, tab2, tab3 = st.tabs(['DNA Translation', 'DNA Sequence Searching', 'Approximate Matching'])
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
    file = st.file_uploader("Fasta file", type=['fasta','fa'],key=1)
    text = None
    if file is not None:
        text = file.getvalue().decode(encoding='utf-8')
        text = text.split('\n')[1]

    st.subheader("This is your DNA sequence:")
    if text is not None:
        selected = -1
        st.write(text)
        pattern = st.text_input("**Pattern required**", placeholder="Enter your pattern to search in UPPERCASE", value=None)
        with st.container():
            options = ['Naive Matching',
            'Bad Character','Good Suffix','Boyer Moore',
            'KMP Searching','K-Meer Searching']
            selected = st.selectbox("Pick a search technique",sorted(options),index=None)
            if 'prev_t' not in st.session_state:
                st.session_state.prev_t = None

        if selected == 'Bad Character' and pattern != None:
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
        elif selected == 'Good Suffix' and pattern != None:
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
        elif selected == 'Boyer Moore' and pattern != None:
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
        elif selected == 'K-Meer Searching' and pattern != None:
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
        elif selected == 'KMP Searching' and pattern != None:
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
        elif selected == 'Naive Matching' and pattern != None:
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
        st.write(":red[No DNA sequence is added yet !]")
with tab3:
    file = st.file_uploader("Fasta file", type=['fasta','fa'],key=2)
    text = None
    if file is not None:
        text = file.getvalue().decode(encoding='utf-8')
        text = text.split('\n')[1]

    st.subheader("This is your DNA sequence:")
    if text is not None:
        selected = -1
        st.write(text)
        pattern2 = st.text_input("**Pattern required**", placeholder="Enter your pattern to search in UPPERCASE", value=None)
        if 'prev_t2' not in st.session_state:
                st.session_state.prev_t2 = None
        if 'cnt2' not in st.session_state:
                    st.session_state.cnt2 = 0
        if st.session_state.prev_t2 != pattern2:
            st.session_state.cnt2 = 0
        
        st.session_state.prev_t2 = pattern2

        if pattern2 != None:
            found, array = approximate_matching(text,pattern2)
            new_text, new_pattern, edits, percentage, ix = interpret_solution_approximate(found, array, 
                                       st.session_state.cnt2, text, pattern2)
            show_index = st.empty()
            show_index.write(f"**Found at index {ix}**")
            st.write("**`Only part of the text will be shown where the pattern is found`**")
            show_text = st.empty()
            show_text.write(new_text + '\n\n' + new_pattern)

            show_per = st.empty()
            show_edits = st.empty()
            show_per.write(f':green[**Matching Percentage %:** ] `{percentage * 100}`')
            show_edits.write(f':red[**Num. of edits performed:** ] {edits}')
            if st.button('Next location', key=3):
                st.session_state.cnt2 += 1
                if st.session_state.cnt2 >= len(found):
                    st.session_state.cnt2 = 0
                new_text, new_pattern, edits, percentage, ix = interpret_solution_approximate(found, array, 
                                       st.session_state.cnt2, text, pattern2)
                show_index.write(f"**Found at index {ix}**")
                show_text.write(new_text + '\n\n' + new_pattern)
                show_per.write(f':green[**Matching Percentage %:** ] `{percentage * 100}`')
                show_edits.write(f':red[**Num. of edits performed:** ] {edits}')
        
    else:
        st.write(":red[No DNA sequence is added yet !]")

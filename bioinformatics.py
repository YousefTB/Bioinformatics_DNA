import numpy as np
import bisect
import time
def naive_matching(t,p):
    found = []
    pattern_length = len(p)
    text_length = len(t)
    num_of_alignments = text_length - pattern_length + 1
    skipped = 0
    for i in range(num_of_alignments):
        if t[i:i+pattern_length] == p:
            found.append(i)
            
    info = {'found in':found, 'performed alignments':num_of_alignments - skipped, 'skipped alignments':skipped}
    return info

def preprocess_pattern_bad_char(p):
    rows = ['A','C','G','T']
    table = np.zeros((len(chars), len(p)), dtype=np.int32)

    for i in range(len(rows)):
        counter = -1
        for j in range(len(p)):
            if rows[i] == p[j]:
                counter = -1
                table[i,j] = counter
            else:
                counter += 1
                table[i,j] = counter
                
    return (table, rows)

def prerocess_pattern_good_suffix(p):
    table = np.full((len(p),2),-1, dtype=np.int32)
    for i in range(1,len(table)):
        remaining = p[:len(p) - i]
        matched = p[len(p) - i:]
        table[i,0] = i
        found = naive_matching(remaining,matched)['found in']
        if found != []:
            index = found[-1]
            table[i,1] = len(p) - i - index - 1
        else:
            for j in range(len(matched), 0, -1):
                ix = len(matched) - j
                if p[:len(matched) - ix] == matched[ix:]:
                    table[i,1] = len(p) - len(matched[ix:]) - 1
                    break
        if table[i,1] == -1:
            table[i,1] = len(p) - 1
    table[0] = [0,1]
    return table

def failure_link_automation(p):
    faliure_link = np.zeros((len(p) + 1, 1), dtype=np.int32)
    for i in range(1,len(p)):
        prefix = p[:i+1]
        suffix = p[1:len(prefix)]
        for j in range(len(suffix)):
            if prefix[:len(suffix[j:])] == suffix[j:]:
                faliure_link[i + 1] = len(suffix[j:])
                break
    return faliure_link

def boyer_moore_bad_char(t,p):
    start = time.time()
    found = []
    skipped = 0
    pattern_length = len(p)
    text_length = len(t)
    num_of_alignments = text_length - pattern_length + 1
    table, rows = preprocess_pattern_bad_char(p)
    
    i = 0
    while(i < num_of_alignments):
        alignments = 0
        if t[i:i+pattern_length] == p:
            found.append(i)
        else:
            for j in range(i + pattern_length - 1, i - 1, -1):
                if t[j] != p[j - i]:
                    key = rows.index(t[j])
                    alignments = table[key, j - i]
                    break

        skipped += alignments
        i += alignments + 1
        
    info = {'found in':found, 'performed alignments':num_of_alignments - skipped, 'skipped alignments':skipped}
    print("Time taken:", time.time() - start)
    return info

def boyer_moore_good_suffix(t,p):
    start = time.time()
    found = []
    skipped = 0
    pattern_length = len(p)
    text_length = len(t)
    num_of_alignments = text_length - pattern_length + 1
    table = prerocess_pattern_good_suffix(p)
    i = 0
    while(i < num_of_alignments):
        match_size = 0
        if t[i:i+pattern_length] == p:
            found.append(i)
        else:
            for j in range(i + pattern_length - 1, i - 1, -1):
                if t[j] != p[j - i]:
                    break
                else:
                    match_size += 1
            
        alignments = table[match_size,1]
        skipped += alignments
        i += alignments + 1
        
    info = {'found in':found, 'performed alignments':num_of_alignments - skipped, 'skipped alignments':skipped}
    print("Time taken:", time.time() - start)
    return info

def boyer_moore(t,p):
    start = time.time()
    found = []
    skipped = 0
    pattern_length = len(p)
    text_length = len(t)
    num_of_alignments = text_length - pattern_length + 1
    table_bad, rows = preprocess_pattern_bad_char(p)
    table_good = prerocess_pattern_good_suffix(p)
    
    i = 0
    while(i < num_of_alignments):
        alignment_bad = 0
        alignment_good = 0
        match_size = 0
        if t[i:i+pattern_length] == p:
            found.append(i)
        else:
            for j in range(i + pattern_length - 1, i - 1, -1):
                if t[j] != p[j - i]:
                    key = rows.index(t[j])
                    alignment_bad = table_bad[key, j - i]
                    break
                else:
                    match_size += 1
        alignment_good = table_good[match_size, 1]
        alignments = max(alignment_good,alignment_bad)
        skipped += alignments
        i += alignments + 1
        
    info = {'found in':found, 'performed alignments':num_of_alignments - skipped, 'skipped alignments':skipped}
    print("Time taken:", time.time() - start)
    return info


def k_meer_search(t,p):
    start = time.time()
    found = []
    length = len(p) - 1
    num_of_alignments = len(t) - len(p) + 1
    text_table = []
    
    for i in range(len(text) - length + 1):
        text_table.append((text[i:i+length], i))
        
    found_counter = 0
    text_table.sort(key=lambda x: x[0])
    keys = [key[0] for key in text_table]
    left_occ = bisect.bisect_left(keys,pattern[:length])
    right_occ = bisect.bisect_right(keys,pattern[:length])
    for i in range(left_occ,right_occ+1):
        if text[text_table[i][1]:text_table[i][1] + len(p)] == pattern:
            found_counter += 1
            found.append(text_table[i][1])
            
    info = {'found in':found, 'performed alignments':found_counter, 'skipped alignments':num_of_alignments - found_counter}
    print("Time taken:", time.time() - start)
    return info

def KMP(t,p):
    start = time.time()
    found = []
    text_length = len(t)
    pattern_length = len(p)
    text = t + ' '
    pattern = p + ' '
    table = failure_link_automation(p)
    i = 0
    j = 0
    while(i != text_length or j != 0):
        if text[i] == pattern[j]:
            i += 1
            j += 1
        elif j == 0:
            i += 1
        else:
            if j == pattern_length:
                found.append(i-j)
            j = table[j,0]
            
    print("Time taken:", time.time() - start)
    return found
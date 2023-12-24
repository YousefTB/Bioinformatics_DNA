import numpy as np
import bisect

def reverse_complement(seq):
    dic = {'G':'C', 'C':'G', 'A':'T', 'T':'A'}
    reversed_seq = list(reversed(seq))
    for i in range(len(reversed_seq)):
        reversed_seq[i] = dic[reversed_seq[i]]
        
    return reversed_seq

def translation(seq):
    dic = {"TTT" : "F", "CTT" : "L", "ATT" : "I", "GTT" : "V",
           "TTC" : "F", "CTC" : "L", "ATC" : "I", "GTC" : "V",
           "TTA" : "L", "CTA" : "L", "ATA" : "I", "GTA" : "V",
           "TTG" : "L", "CTG" : "L", "ATG" : "M", "GTG" : "V",
           "TCT" : "S", "CCT" : "P", "ACT" : "T", "GCT" : "A",
           "TCC" : "S", "CCC" : "P", "ACC" : "T", "GCC" : "A",
           "TCA" : "S", "CCA" : "P", "ACA" : "T", "GCA" : "A",
           "TCG" : "S", "CCG" : "P", "ACG" : "T", "GCG" : "A",
           "TAT" : "Y", "CAT" : "H", "AAT" : "N", "GAT" : "D",
           "TAC" : "Y", "CAC" : "H", "AAC" : "N", "GAC" : "D",
           "TAA" : "*", "CAA" : "Q", "AAA" : "K", "GAA" : "E",
           "TAG" : "*", "CAG" : "Q", "AAG" : "K", "GAG" : "E",
           "TGT" : "C", "CGT" : "R", "AGT" : "S", "GGT" : "G",
           "TGC" : "C", "CGC" : "R", "AGC" : "S", "GGC" : "G",
           "TGA" : "*", "CGA" : "R", "AGA" : "R", "GGA" : "G",
           "TGG" : "W", "CGG" : "R", "AGG" : "R", "GGG" : "G" 
           }
    
    protien = ''
    flag_start = False
    for i in range(0,len(seq)- 2, 3):
        if dic[seq[i:i+3]] == 'M':
            flag_start = True
        elif dic[seq[i:i+3]] == '*':
            flag_start = False
            
        if flag_start:
            protien = protien + dic[seq[i:i+3]]

    return protien

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
    table = np.zeros((len(rows), len(p)), dtype=np.int32)

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
    return info

def boyer_moore_good_suffix(t,p):
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
    return info

def boyer_moore(t,p):
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
    return info


def k_meer_search(text,pattern):
    found = []
    length = len(pattern) - 1
    num_of_alignments = len(text) - len(pattern) + 1
    text_table = []
    
    for i in range(len(text) - length + 1):
        text_table.append((text[i:i+length], i))
        
    found_counter = 0
    text_table.sort(key=lambda x: x[0])
    keys = [key[0] for key in text_table]
    left_occ = bisect.bisect_left(keys,pattern[:length])
    right_occ = bisect.bisect_right(keys,pattern[:length])
    for i in range(left_occ,right_occ+1):
        if text[text_table[i][1]:text_table[i][1] + len(pattern)] == pattern:
            found_counter += 1
            found.append(text_table[i][1])
            
    info = {'found in':found, 'performed alignments':found_counter, 'skipped alignments':num_of_alignments - found_counter}
    return info

def KMP(t,p):
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
            

    return found

def approximate_matching(text,pattern):
    t = '*' + text
    p = '*' + pattern
    array = np.zeros((len(p), len(t)), dtype=np.int32)
    array[:,0] = np.arange(len(p))
    for i in range(1,len(p)):
        for j in range(1,len(t)):
            delta = 1 if p[i] != t[j] else 0
            array[i,j] = min(array[i - 1, j - 1] + delta, array[i, j - 1,] + 1, array[i - 1, j] + 1)

    found = np.argwhere(array[-1] == min(array[-1])).flatten().tolist()
    return (found, array)

def interpret_solution_approximate(found_array, dp_array, index, text, pattern):
    t = '*' + text
    p = '*' + pattern
    solution = []
    j = found_array[index]
    i = len(p) - 1
    while(i != 0):
        ix1 = dp_array[i, j - 1] + 1
        ix2 = dp_array[i-1, j] + 1
        delt = 1 if p[i] != t[j] else 0
        ix3 = dp_array[i - 1, j - 1] + delt
        solution.append((i,j,dp_array[i,j]))
        if dp_array[i,j] == ix1:
            j = j -1
        elif dp_array[i,j] == ix2:
            i = i - 1
        elif dp_array[i,j] == ix3:
            i = i - 1
            j = j - 1
            
    solution = list(reversed(solution))
    new_pattern = ''
    new_text = t[1:solution[0][1] + 1]
    stop = 0
    prev_path = solution[0]
    new_pattern += p[prev_path[0]]
    new_text += t[prev_path[1]]
    for i in range(1,len(solution)):
        curr_path = solution[i]
        if curr_path[2] > prev_path[2] and curr_path[0] == prev_path[0]:
            new_pattern += '-'
            new_text += t[curr_path[1]]
            prev_path = curr_path
            continue
        elif curr_path[2] > prev_path[2] and curr_path[1] == prev_path[1]:
            new_text += '-'
            new_pattern += p[curr_path[0]]
            prev_path = curr_path
            continue
        new_pattern += p[curr_path[0]]
        new_text += t[curr_path[1]]
        prev_path = curr_path


    new_text += t[solution[-1][1] + 1:]
    edits = solution[-1][2]
    percentage = (len(pattern) - edits) / len(pattern)
    return (new_text[solution[0][1]: solution[-1][1] + 30], new_pattern, edits, percentage, solution[0][1] - 1)

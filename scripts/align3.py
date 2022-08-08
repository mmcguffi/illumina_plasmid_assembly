from collections import Counter
from itertools import groupby
from operator import itemgetter

import pandas as pd
from Bio import Align


def rotate(strg, n):
    n = n % len(strg)
    return strg[n:] + strg[:n]

def slice_plas(strg, start, size, direction):
    # if direction == -1:
    #     strg = strg.reverse_complement()
    return strg[start:(start+size)]

# def rm_not_found(hits):
#     return [(k[0],k[1],v) for k,v in hits.items() if k[0] >= 0]

def get_chunk_locs(ref, query, drct, chunk_size):
    #CHUNK_SIZE = 6
    #hits = Counter([(query.find(slice_plas(ref, i, CHUNK_SIZE,  drct)) - i,  drct) for i in range(len(ref))])
    if drct == -1:
        query = query.reverse_complement()
    hits = [(query.find(slice_plas(ref, i, chunk_size,  drct)) - i,  drct, i) for i in range(len(ref))]
    hits = list(filter(lambda _: _[0] >= 0, hits))

    return hits
    return rm_not_found(hits)

def get_adj_offest_dir(ref, query, chunk_size):
    
    fwd = get_chunk_locs(ref, query, 1, chunk_size)
    rev = get_chunk_locs(ref, query, -1, chunk_size)
    
    try: 
        best = Counter([(_[0],_[1]) for _ in fwd+rev]).most_common()[0]
        best_amt = best[0][0]
        strand = best[0][1]

        best_correct_strand = list(filter(lambda _: _[1] == strand, fwd+rev))
        offset_of_offest = min(list(filter(lambda _: _[0] ==best_amt, best_correct_strand)), key=lambda _: _[2])
        offest_adj = offset_of_offest[0] + offset_of_offest[2]
    
    except IndexError:
        offest_adj = 0
        strand = 1
    
    return (offest_adj, strand)

def get_indels_circular_canon(query, reference, offset, direction):
    
    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    aligner.match = 2
    aligner.mismatch = -5
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -4
    aligner.target_end_gap_score = 0.0
    aligner.query_end_gap_score = 0.0
    
    if direction == 1:
        pass
    elif direction == -1:
        query = query.reverse_complement()
    else:
        raise ValueError("direction must be '1' or '-1'")

    align = aligner.align(reference, rotate(query, offset))[0]

    return align

def calc_indel_df(align):

    
    align = format(align, '').split()
    inserts = [i for i,_ in enumerate(align[0]) if _ == '-']
    mismatches = [i for i,_ in enumerate(align[1]) if _ == '.']
    deletions = [i for i,_ in enumerate(align[2]) if _ == '-']

    #reduce all indices of '-' to ranges of consecutive '-'
    ranges =[]
    for Type, pos in {"insert":inserts, "deletion":deletions, "mismatch":mismatches}.items():
        for k,g in groupby(enumerate(pos), lambda x : x[0] - x[1]):
            group = (map(itemgetter(1), g))
            group = list(map(int, group))
            ranges.append((group[0], group[-1] + 1, Type))
    
    indels = pd.DataFrame(ranges, columns=['qstart','qend','Type'])
    indels = indels.sort_values(by=['qstart','qend'])
    
    return indels#, rotation, direction

def get_indels(sequencing, ref):
    df_aligns = []
    for chunk_size in [16,9,7,6,5]:
        offest_adj, strand = get_adj_offest_dir(ref, sequencing, chunk_size)
        align = get_indels_circular_canon(sequencing, ref, offest_adj, strand)
        df = calc_indel_df(align)
        if df.empty:
            import streamlit as st
            st.write(chunk_size)
            return df
        df_score = sum(abs(df['qstart'] - df['qend'])) 
        df_aligns.append((df_score, df))
    best_df = min(df_aligns, key=lambda _: _[0])[1]

    return best_df

from Bio import SeqIO
import pandas as pd
from Bio import SeqIO

import bokeh_plot as bk
import align3 as al

def gbk_to_df(sg):
    
    #gbk['seq'] = str(sg.seq)
    gbks = []
    for feature in sg.features:
        gbk = {}
        
        if feature.type == "source":
            continue
        
        gbk['Type'] = feature.type
        
        parts = feature.location.parts
        if len(parts) == 1:
                   
            gbk['qstart'] = feature.location.start
            gbk['qend'] = feature.location.end
        elif len(parts) == 2:

            s1 = feature.location.parts[0].start
            e1 = feature.location.parts[1].end
            
            s2 = feature.location.parts[1].start
            e2 = feature.location.parts[0].end
            
            if (e1-s2) > (e2-s1):
                s = s2
                e = e2
            else:
                s = s1
                e = e1
                
            gbk['qstart'] = s
            gbk['qend'] = e
            
        else:
            continue
    
 
        gbk['qlen'] = len(sg.seq) #gbks['qend'] - gbks['qstart']

        # if (gbk['qend'] - gbk['qstart']) == gbk['qlen']:
        #     continue
        
        gbk['sframe'] = feature.location.strand
        
        for k,v in feature.qualifiers.items():
            gbk[k] = ', '.join([str(elem) for elem in v])
        
        gbks.append(gbk)

    gbks = pd.DataFrame(gbks)
    return gbks


gbk  = SeqIO.read(snakemake.input[0], "genbank")
gbk2 = SeqIO.read(snakemake.input[1], "genbank")
# gbk = SeqIO.read("references/MCG1_S281.gbk", "genbank")
# gbk2  = SeqIO.read("output/annotated_plasmids/MCG1_S281_pLann.gbk", "genbank")
inSeq = str(gbk.seq)

gbk_df = gbk_to_df(gbk)

try:
    gbk_df['score'] = abs(gbk_df['qend'] - gbk_df['qstart']) * (gbk_df['identity'].astype(float)/100) * (gbk_df['match_length'].astype(float)/100)
except KeyError:
    gbk_df['score'] = abs(gbk_df['qend'] - gbk_df['qstart'])
gbk_df['pi_permatch'] = gbk_df['score']
gbk_df['db'] = "mutations"

# x = pd.read_csv("/Users/mmcguffi/projects/pLannotate/tests/test_data/pXampl3.csv")


try:
    gbk_df['fragment'] = gbk_df['fragment'].map({"True":True, "False":False})
except KeyError:
    gbk_df['fragment'] = False
gbk_df['fragment'] = gbk_df['fragment'].fillna(False)
    
indels = al.get_indels(gbk2.seq, gbk.seq)
indels['fragment'] = False
indels = indels.apply(pd.to_numeric, errors='ignore')

for index in indels.index:
    start = indels.loc[index]['qstart']
    end = indels.loc[index]['qend']
    Type = indels.loc[index]['Type']
    seq = "n" * (indels.loc[index]['qend'] - indels.loc[index]['qstart']) #indels.loc[index]['qseq']
    
    if Type == 'insert':
        left = inSeq[ : start]
        right = inSeq[start : ]
        inSeq = left + seq + right
        
    elif Type == 'deletion':
        pass
        # left = inSeq[ : start]
        # right = inSeq[end : ]
        # inSeq = left + right
    elif Type == 'mismatch':
        pass
    else:
        raise ValueError("Type not recognized")
    
    
# if indels.empty:
#     st.success("Identical sequences.")
# elif sum(abs(indels['qstart'] - indels['qend'])) >= .95 * len(inSeq):
#     st.error("Likely incorrect reference sequence.")
# else:
#     st.warning("Mutations detected.")

gbk_df = gbk_df.append(indels)
gbk_df['qlen'] = len(inSeq)
gbk_df = gbk_df.reset_index(drop=True)

bk.save_plot(gbk_df, snakemake.output[0])

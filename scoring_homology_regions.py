#!/opt/anaconda3/bin/python
import pandas as pd

with open("covid19_to_other_6_sars_consensus.csv","r+") as f:
    data = [x.replace(" ","").replace("[]","").replace("\n","") for x in f.readlines()]
f.close()

data= [x.split(",") for x in data]
data = data[0]
data = data[1:]

data=[float(x) for x in data]

# read in positions of reference sequence (or else won't align correctly!)
with open("NC_045512_alignment_sequence.txt","r+") as f:
    ref_seq =[x.replace("\n","") for x in f.readlines()]
f.close()
ref_seq = ("".join(ref_seq))


# remap sequence positions
# remove positions associated with "-" regions
data_new = []
# ref_seq = list(ref_seq)
posn_nums = []
i=0
j=0
while i<len(ref_seq):
    if(ref_seq[i]=="-"):
        data_new.append(0)
        posn_nums.append(-1)
    else:
        data_new.append(data[i])
        (posn_nums.append(j+1))
        j+=1
    i+=1

data=data_new

def compute_ranges(start,stop,data,gene,cutoff_pcntl,posn_nums):
    req_region_num = 12

    seq_len = 16

    start = posn_nums.index(start)
    stop = posn_nums.index(stop)

    # create sliding window of seq_len and identify regions where average % identity is high
    import statistics
    avg_idenity_ls = []
    index_ls = []
    '''
    # homology % cutoff
    hm = 80
    ct = 7
    while(len(index_ls) < 3):
        index_ls = []
        i = start
        while i<stop:
            # positions 1-7 homology > hm %
            count = len([i for i in (data[i:i+7]) if i >= hm])
            if(count>=ct):
                # positions 1-16 homology > hm-20 % (less stringent for rest of sequence)
                count = len([i for i in (data[i+7:i+seq_len]) if i >= hm-20])
                if(count>=10):
#                 print(posn_nums[i])
#                 print(data[i:i+7])
                    index_ls.append(posn_nums[i])#+start])
            i+=1
        hm-=5
        ct-=0.5
    '''
    hm_score = []
    index_ls = []
    i = start
    while i<stop:
        #hm_score.append(int(len([i for i in (data[i:i+seq_len]) if i >= 70])))
        #hm_score.append(int(sum(data[i:i+seq_len]))) # sum homology for that region
        hm_score.append(int(statistics.mean(data[i:i+seq_len]))) # sum homology for that region
        index_ls.append(posn_nums[i])#+start])
        i+=1






    output_data_str = ("gene:"+gene)




    return(output_data_str,index_ls,hm_score)



# protein,start position,end position
gene_data =[['pp1a',266,13483],
            ['pp1ab',266,21555,],
            ["spike",21563,25384],
            ["3a",25393,26220],
            ["envelope",26245,26472],
            ["matrix",26523,27191],
            ["7a",27394,27759],
            ["8b",27894,28259],
            ["nucleocapsid",28274,29533]
           ]

output_data = []
for g in gene_data:
    output_data.append(compute_ranges(g[1],g[2],data,g[0],99,posn_nums))


# output to file
with open("homology_scores.csv","w+") as f:
    i = 0
    while i<len(output_data): # loop through each gene
        j = 0
        f.write(output_data[i][0])
        f.write("\n")
        while j<len(output_data[i][1]):
            f.write(str(output_data[i][1][j]).replace("[","").replace("]",""))
            f.write(",")
            f.write(str(output_data[i][2][j]).replace("[","").replace("]",""))
            f.write("\n")
            j+=1
        i+=1
f.close()

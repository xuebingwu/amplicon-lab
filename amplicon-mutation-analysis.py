'''
A python3 script for analyzing error-prone PCR amplicon sequencing data

Input:
    - A fastq.gz file from amplicon-sequencing of error-prone PCR products.
    - Most reads should be of identical length. Shorter or longer reads with indels will be discarded.
    - assume mutation rate is <50% per position
    
Output:
    - plots for mutation rate at each position
    
for file in *_R1_*.fastq.gz
do
	python ~/xuebing/UTR-mutagenesis/amplicon-mutation-analysis.py $file $file
done

'''

import copy
import gzip
import numpy as np

from matplotlib import pyplot as plt
#%matplotlib inline
import matplotlib.gridspec as gridspec



fastq_gz_file='/media/rna/sequencing_data/chio-lab/dms-pilot-20230512/10128-CR31_R1_001.fastq.gz'
#fastq_gz_file='/home/xw2629/chio-lab/DMSseqT4_R1_001.fastq.gz'
fastq_gz_file='/media/rna/sequencing_data/peiguo_shi/DMS-2023-0410-silvestrol/00_fastq/LM2-Silvestrol-DMS-2_R1_001.fastq.gz'
sample_name='test'

import sys

fastq_gz_file=sys.argv[1]
sample_name = sys.argv[2]

# the script automatically determines the length and reference sequences from a fastq.gz file
# reads may have different lengths

#######################
# step 1: process a fraction of the reads to determine the normal length
#######################

num_seq_to_read = 100

# a dict: length -> counts
seq_length = {} 

# a counter for how many sequences have been read
num_seq = 0
with gzip.open(fastq_gz_file, 'rb') as f:
    line_nb = 0 # line number
    for line in f:
        line_nb = line_nb + 1
        if(line_nb % 4 == 2): # in fastq format, only 1 of every 4 lines is sequence
            seq = line.decode().strip()
            l = str(len(seq)) # the length of the sequence
            if not l in seq_length:
                seq_length[l] = 0
            seq_length[l] = seq_length[l] + 1
            num_seq = num_seq + 1
            if num_seq == num_seq_to_read:
                break
                
# find the most frequent read length (mode)
sorted_d = sorted(seq_length.items(), key=lambda x: x[1],reverse=True)
mode_len = int(sorted_d[0][0])
#print(sorted_d)
print(str(sorted_d[0][1]/num_seq_to_read*100)+'% of the first '+str(num_seq_to_read)+' reads have a length of '+str(mode_len))
#print(mode_len,sorted_d[0][1],num_seq_to_read)



#######################
# step 2: count nucleotide frequency at each position
#######################

counts = {}
counts['A']=[0]*mode_len # the number of reads with 'A' at each position
counts['C']=[0]*mode_len
counts['G']=[0]*mode_len
counts['T']=[0]*mode_len
counts['N']=[0]*mode_len

# unique sequenes and read counts
uniq_seq = {}

n_seq = 0
with gzip.open(fastq_gz_file, 'rb') as f:
    line_nb = 0
    for line in f:
        line_nb = line_nb + 1
        if(line_nb % 4 == 2):
            seq = line.decode().strip()
            if len(seq) == mode_len:
                n_seq = n_seq + 1
                if not seq in uniq_seq:
                    uniq_seq[seq] = 0
                uniq_seq[seq] = uniq_seq[seq] + 1
                for i in range(mode_len):
                    counts[seq[i]][i] = counts[seq[i]][i] + 1

print(str(n_seq)+" of "+str(round(line_nb/4))+" ("+str(round(n_seq*4/line_nb*10000)/100)+"%)"+" reads have a length of "+str(mode_len))
print(str(len(uniq_seq))+ ' unique sequences')
sorted_d = sorted(uniq_seq.items(), key=lambda x: x[1],reverse=True)

print("top 3 most abundant sequences and their counts:")
print(sorted_d[:3])


#######################
# step 3: count nucleotide frequency at each position and calculate mutation rate
#######################

consensus_seq =['']*mode_len     # the consensus sequene, likely the reference/wt 
mutation_rate = copy.deepcopy(counts) # rate of mutation to each nucleotide at each position
avg_mutation_rate = [0]*mode_len # average mutation rate (total) at each position

for base in mutation_rate.keys():
    for i in range(mode_len):
        mutation_rate[base][i] = mutation_rate[base][i] / n_seq
        if mutation_rate[base][i] > 0.5: # assuming mutation rate is <50% at each position
            mutation_rate[base][i] = 0
            consensus_seq[i] = base
        else:
            avg_mutation_rate[i] = avg_mutation_rate[i] + mutation_rate[base][i]

            
fig = plt.figure(tight_layout=True,figsize=(14, 14), dpi=80)
gs = gridspec.GridSpec(4, 2)

ax = fig.add_subplot(gs[0, :])

# plot mutation rate at each position
#plt.figure(figsize=(14, 4), dpi=80)
for base in ['T','A','C','G']:
    ax.step(range(mode_len),mutation_rate[base],label=base)
    ax.set_ylim([0,0.1])
ax.legend()
ax.set_title('Fig. 1. Rate of mutation to each nucleotide at each position')
ax.set_ylabel('Mutation rate to each nucleotide')
ax.set_xlabel('Position along the sequence')
#plt.savefig(sample_name+'-mutation_to_each_nucleotide.pdf')

ax = fig.add_subplot(gs[1, :])
ax.step(range(mode_len),avg_mutation_rate)
ax.set_ylim([0,0.1])
ax.set_ylabel('Total mutation rate')
ax.set_xlabel('Position')
ax.set_title('Fig. 2. Total mutation rate at each position')

#plt.savefig(sample_name+'-mutation_to_all_nucleotide.pdf')

# count mutations per read and mutation type
mutation_per_read = []
base2num = {'A':0,'C':1,'G':2,'T':3}
mutation_matrix = np.zeros((4,4))
mutation_dict = {}
with gzip.open(fastq_gz_file, 'rb') as f:
    line_nb = 0
    for line in f:
        line_nb = line_nb + 1
        if(line_nb % 4 == 2):
            seq = line.decode().strip()
            if len(seq) == mode_len:
                mutation_per_read.append(0)
                for i in range(mode_len):
                    if seq[i] != consensus_seq[i]:
                        mutation_per_read[-1] = mutation_per_read[-1] + 1
                        
                        if seq[i] != 'N':
                            m = consensus_seq[i]+'>'+seq[i]
                            if not (m in mutation_dict):
                                mutation_dict[m] = 0
                            mutation_dict[m] = mutation_dict[m] + 1
                            mutation_matrix[base2num[consensus_seq[i]]][base2num[seq[i]]] = mutation_matrix[base2num[consensus_seq[i]]][base2num[seq[i]]] + 1
                
# plot a histogram for the number of mutations per read
ax = fig.add_subplot(gs[2, 0])
n, bins, patches = ax.hist(mutation_per_read, max(mutation_per_read), density=True, facecolor='w', alpha=0.99,color='w')
ax.set_xlim(0,9)
dens = n[:8]
dens[-1] = sum(n[7:])
nm = ['0','1','2','3','4','5','6','>6']
ax.bar(nm,dens)
ax.set_xlim(-0.5,8)
ax.set_ylabel('Fraction of reads')
ax.set_xlabel('Number of mutations per read')
ax.set_title('Fig. 3. Number of mutations per read. Read length = '+str(mode_len))
#plt.savefig(sample_name+'-mutation_per_read.pdf')


# plot mutation matrix
mutation_matrix = mutation_matrix/sum(sum(mutation_matrix))
label = ['A','C','G','T']
ax = fig.add_subplot(gs[2, 1])
im = ax.imshow(mutation_matrix,cmap='Greens') # winter
# Show all ticks and label them with the respective list entries
ax.set_xticks(np.arange(len(label)))
ax.set_yticks(np.arange(len(label)))
ax.set_xticklabels(label)
ax.set_yticklabels(label)
ax.set_xlabel('To')
ax.set_ylabel('From')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)
# Loop over data dimensions and create text annotations.
for i in range(len(label)):
    for j in range(len(label)):
        text = ax.text(j, i, round(mutation_matrix[i,j],ndigits=2),
                       ha="center", va="center", color="w")

ax.set_title("Fig. 4. Fraction of all mutations")
#fig.tight_layout()
#plt.show()
#fig.savefig(sample_name+'-matation_matrix.pdf')

# barplot of mutation types
ax = fig.add_subplot(gs[3, 0])
ax.bar(mutation_dict.keys(),mutation_dict.values())
ax.set_xlabel("Mutation type")
ax.set_xlabel("Frequency")
ax.set_title("Fig. 5. Frequency of mutation types")
plt.savefig(sample_name+'.pdf')

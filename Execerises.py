from Bio import pairwise2
#Question 1
#importing the sequence for exercise_1
# the file was downloaded to my laptop from the following NCBI source
# https://www.ncbi.nlm.nih.gov/nuccore/NM_002520.7?report=fasta&log$=seqview
sequence_ex1 = open('/Users/duniaalmelhm/Downloads/sequence.fasta', 'r')
lines_seq_ex1 = sequence_ex1.readlines()[1:20] #reading all lines in the fasta file
#len(lines_seq_ex1) gives 21 lines
string_seq_exc1 = "".join(lines_seq_ex1) #unlisting to make the seq as a string

#### part 1
length_nt = len(string_seq_exc1)
print(length_nt) #1339
length_codons = length_nt/3
print(length_codons) #446 codons and 1 nucleotide

#### part 2
GC_content = 0
for i in string_seq_exc1:
    if i in 'GC':
        GC_content += 1
percent_GC_content = GC_content/int(length_nt)*100
print(percent_GC_content) #40.70% approximately

#### part 3
motif = "CTTAGTAGCTGTGGAGGAA"
if motif in string_seq_exc1:
    print("yes")
print(string_seq_exc1.count(motif)) # 1 time
print(string_seq_exc1.find(motif)) # at start position/index 450

#Question 2
sequence_ex2 = open('/Users/duniaalmelhm/Downloads/exercise_2_sequence.fasta', 'r')
lines_seq_ex2 = sequence_ex2.readlines()[1:14] #number of lines in this file is 15
string_seq_exc2 = "".join(lines_seq_ex2)
alignments = pairwise2.align.globalxx(string_seq_exc1, string_seq_exc2, score_only=True)
print(alignments) #886.0nt is the alignment score
comparing_score_seq1 = int(alignments)/length_nt*100
print(comparing_score_seq1) #66.17%
comparing_score_seq2 = int(alignments)/len(string_seq_exc2)*100
print(comparing_score_seq2) #98.34%
print(len(string_seq_exc2)) #901nt
#we notice more than 50% of the sequence in question 1 is similar to the sequence given in question 2 when taking into
#consideration the full length of sequence question 1. The similarity of these two sequences can be determined by the
#alignement score which is 886.0nt. A higher alignment score indicates more similarity as a known fact.








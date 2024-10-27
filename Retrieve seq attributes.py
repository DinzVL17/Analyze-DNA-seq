'''
Task 4: Load the downloaded FASTA file as a sequence record object and print attributes of the record
Author: Dinuri Vishara
Date:4/1/2023
'''

#import packages
from Bio import SeqIO
import re
#load the FASTA file
for seq_record in SeqIO.parse("ATdreb2a.fasta", "fasta"):
    #print the attributes of the record
    print(seq_record.id)
    print(seq_record.description)
    print(repr(seq_record.seq))
    print(len(seq_record))

# run the web-based nucleotide blast program on the ATDREB2A sequence
from Bio.Blast import NCBIWWW
from Bio import SeqIO
record = SeqIO.read("ATdreb2a.fasta", format="fasta")
result_handle = NCBIWWW.qblast("blastn", "nt", record.seq)

# Save the blast output as “dreb2a_blast.xml”
with open ('dreb2a_blast.xml','w') as file1:
    file1.write(result_handle.read())

# select the hits that are below the E-value threshold of 0.05
E_VALUE_THRESH = 0.05
# open the xml file
outputfile = open("dreb2a_blast.xml")
from Bio.Blast import NCBIXML
blast_records = NCBIXML.read(outputfile)

# Print the attributes of each blast hit
count = 0
for alignment in blast_records.alignments:
    for hsp in alignment.hsps:
        if hsp.expect < E_VALUE_THRESH:
            print("****Alignment****")
            print("sequence:", alignment.title)
            print("length:", alignment.length)
            print("e value:", hsp.expect)
            print("Score:", hsp.score)
            print("subject sequence:",hsp.sbjct[0:75] + "...")
            print("hit sequence length:", len(hsp.sbjct))

            # define a search string to search for the ABRE element
            pattern= "[TC]ACGT[GT]C"
            # check each BLAST hit for the element
            mo = re.finditer(pattern, hsp.sbjct)

            for item in mo:
                # print the detected sequence fragment and location
                print("Element:" , item.group())
                print("Sequence location:", item.span())
                # get the number of BLAST hits with ABRE element
                count+=1

print("Total count:", count)


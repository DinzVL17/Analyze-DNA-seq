'''
Task2: Writing custom Python methods to analyze DNA sequences
Author: Dinuri Vishara lokuwalpola
Date : 30/11/2022
'''

# calculate the AT content of a DNA sequence
def countAT(inputseq):
    ATcount = 0

    for base in inputseq:
        if base == "A" or base =="T":
            ATcount +=1

    return(ATcount)


#split multiple FASTA sequences in a single text file and return a dictionary
def splitSeq(file):
    seq={}
    # Open the FASTA file
    with open (file,'r') as f1:
        for line in f1:
            # strip lines
            lines = line.strip()
            if lines != '\n':
                # store FSTA headers in header variable
                if ">" in lines:
                    seq1 = ""
                    header = lines
                # store FASTA sequence in seq1 variable
                else:
                    seq1 += lines
                    # create the dictionary
                    seq[header] = seq1
        return(seq)


# check the type of a given sequence is DNA, mRNA or amino acid sequence
def checkSeqType(seq):

    # create a list of amino acids
    aa = [ "R", "N", "D","B","E", "Q","Z","H","I","L","K","M","F","P","S","W","Y","V"]
    for base in seq:
        # check whether characters in sequence are equal to amino acids in list
        if base in aa:
            return ("amino acid sequence")
        elif base not in aa and "U" in seq:
            return ("mRNA sequence")
        elif base not in aa and "U" not in seq:
            return ("DNA sequence")


#split the FASTA file and create a cictionary
diction = splitSeq("OSDREB_sequences.FASTA")
print(diction)

seqInfo = {}
# check the sequence type of the sequences
for value in diction.values():
    seqType = checkSeqType(value)
    #if the sequnce is a DNA sequence calculate the AT content
    if seqType == "DNA sequence":
        ATcount =countAT(value)
        # get the percentage of the ATcontent
        ATpercentage = (ATcount / len(value) *100)
        # create a dictionary with sequence as the key and sequence type plus AT content as the value
        seqInfo [value] = seqType, ATpercentage

print(seqInfo)





















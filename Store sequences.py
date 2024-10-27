'''
Task3: Extend the sequence class and write subclasses for DNA, mRNA, and amino acid sequences
Date : 7/12/2022
Author: Dinuri Vishara
'''

#Sequence class
class Sequence:
    Sequence_count = 0

    def __init__(self, Gene_name, Gene_ID, Species_name, Subspecies_name,sequence):
        self.Gene_ID = Gene_ID
        self.Gene_name = Gene_name
        self.Sequence_type = Sequence.get_Seq_Type(Sequence,sequence)
        self.Sequence_length = len(sequence)
        self.Species_name = Species_name
        self.Subspecies_name = Subspecies_name
        Sequence.Sequence_count += 1

    @staticmethod
    def fasta_Split(file):
        seqDiction = {}
        seq1 =""
        header =""
        with open(file, 'r') as file1:
            for line in file1:
                lines = line.strip()
                if lines != '\n':
                    if ">" in lines:
                        seq1 = ""
                        headerlist = []
                        lines = lines.split("-")
                        header = lines[0]
                        headerlist = lines
                    else:
                        seq1 += lines
                        seqDiction[header] = headerlist + [seq1]
        return(seqDiction)


    def get_Seq_Type(self, seq):
        aa = ["R", "N", "D", "B", "E", "Q", "Z", "H", "I", "L", "K", "M", "F", "P", "S", "W", "Y", "V"]
        for base in seq:
            if base in aa:
                return ("amino acid sequence")
            elif base not in aa and "U" in seq:
                # seqType = "mRNA sequence"
                return ("mRNA sequence")
            elif base not in aa and "U" not in seq:
                # seqType = "DNA sequence"
                return ("DNA sequence")


    def get_Character_Count(self, seq):
        characterCount = {}
        for base in seq:
            if base in characterCount:
                characterCount[base] += 1
            else:
                characterCount[base] = 1
        return (characterCount)

# DNA subclass
class DNAseq(Sequence):

    # constructer method
    def __init__(self, Gene_name,Gene_ID,  Species_name, Subspecies_name,sequence):
        super().__init__( Gene_name, Gene_ID,  Species_name, Subspecies_name, sequence)

        self.AT_content = DNAseq.get_ATcontent(DNAseq,sequence)
        self.transcribed_seq = DNAseq.transcribe_Sequence(DNAseq, sequence)

    # method to transcribe a DNA sequence into its mRNA sequence
    def transcribe_Sequence(self, sequence):
        self.transcribed_seq = sequence.replace("T","U")
        return (self.transcribed_seq)

    # method to get the AT content of a given sequence
    def get_ATcontent(self, sequence):
        ATcount = 0
        for base in sequence:
            if base == "A" or base == "T":
                ATcount += 1
                self.AT_content = str(ATcount/len(sequence)*100)
        return (self.AT_content)


# mRNA class
class MRNAseq(Sequence):
    __Amino_acid_codons = {}

    # constructer method
    def __init__(self, Gene_name, Gene_ID, Species_name, Subspecies_name, sequence):
        super(MRNAseq, self).__init__(Gene_name, Gene_ID,  Species_name, Subspecies_name, sequence)

        self.AT_content = MRNAseq.get_ATcontent(MRNAseq, sequence)
        self.translated_seq = MRNAseq.translate_Sequence(MRNAseq, sequence)

    # method to get AT content of a sequence
    def get_ATcontent(self, sequence):
        ATcount = 0
        sequence1 = sequence.replace("U", "T")
        for base in sequence1:
            if base == "A" or base == "T":
                ATcount += 1
                self.AT_content = str(ATcount/len(sequence)*100)
        return (self.AT_content)

    # method to store codon-amino acid pairs from a text file
    @classmethod
    def upload_Codons(cls, file):
        with open(file, 'r') as file1:
            for line in file1:
                if "#" not in line:
                    newline = line.strip().split()
                    codon = newline[0]
                    AminoA = newline[2]
                    cls.__Amino_acid_codons[codon] = AminoA
            return (cls.__Amino_acid_codons)

    # method to translate a given mRNA sequence into its amino acid sequence
    def translate_Sequence(self, sequence):
        arr = []
        for i in range(0, len(sequence), 3):
            c = (sequence[i: i+3])
            arr.append(c)
        num = 0
        aa = ""
        while num < len(arr):
            c = arr[num]
            aa += str(self.__Amino_acid_codons.get(c))
            self.translated_seq = aa
            if c in ["UAA", "UAG", "UGA"]:
                break
            num += 1
        return (self.translated_seq)

# protein subclass
class Proteinseq(Sequence):

    # constructer method
    def __init__(self, Gene_name, Gene_ID, Species_name, Subspecies_name,
                 Uniprot_ID, Reviewed_status, sequence):
        super(Proteinseq, self).__init__(Gene_name, Gene_ID, Species_name, Subspecies_name, sequence)

        self.uniprotID = Uniprot_ID
        self.reviewedS = Reviewed_status
        self.hydrophobicity = Proteinseq.get_Hydrophobicity(Proteinseq,sequence)

    # get the percentage of the total hydrophobic amino acid residues in a sequence
    def get_Hydrophobicity(self, sequence):
        aa = 0
        percentage = 0
        hydrophobicAcids = ["A", "I", "L", "M", "F", "W", "Y", "V"]
        for i in sequence:
            if i in hydrophobicAcids:
                aa += 1
                percentage = (aa / len(sequence) * 100)
        return (percentage)


if __name__ == "__main__":
    # create dictionary using FASTA split method
    diction = Sequence.fasta_Split("OSDREB_sequences.FASTA")

    # get values of the dictionary and store in seqInfo list
    for values in diction:
        seqInfo = diction[values]

        # check the type of the sequence in seqInfo list
        seqType = Sequence.get_Seq_Type(Sequence,seqInfo[-1])
        # If the sequence type is DNA, create object of DNAseq subclass
        if seqType == "DNA sequence":
            DNA = DNAseq(*seqInfo)
            # get details for the OSDREB1A DNA sequence
            if DNA.Gene_name == ">DREB1A_CDS":
                print("Gene ID of DREB1A sequence:", DNA.Gene_ID)
                print("Sequence length of DREB1A sequence:", DNA.Sequence_length)
                print("Sequence type of DREB1A sequence:", DNA.Sequence_type)
                print("AT content of DREB1A sequence:", DNA.get_ATcontent(seqInfo[-1]))

            # transcribe the OSDREB2B coding sequence
            if DNA.Gene_name == ">DREB2B_CDS":
                DREB2B_transcribed = DNA.transcribe_Sequence(seqInfo[-1])

                # create a new object for the resulting mRNA sequence.
                mRNA = MRNAseq(seqInfo[0],seqInfo[1],seqInfo[2],seqInfo[3], DREB2B_transcribed)
                if mRNA.Gene_name == ">DREB2B_CDS":
                    print("Length of DREB2B sequence: ", mRNA.Sequence_length)
                    print("Type of DREB2B sequence: ",mRNA.Sequence_type)
                    print("AT content of DREB2B sequence: ", mRNA.get_ATcontent(seqInfo[-1]))
                    print("mRNA sequence of DREB2B: ", DREB2B_transcribed)

                    # translate the OSDREB2B mRNA sequence
                    mRNA.upload_Codons("codon_table.txt")
                    DREB2B_translated = mRNA.translate_Sequence(DREB2B_transcribed)
                    print("Translated sequence of DREB2B: ", DREB2B_translated)
                    print("Length of translated DREB2B sequence: ", len(DREB2B_translated))

        # If the sequence type is protein, create object of Proteinseq subclass
        elif seqType == "amino acid sequence":
            protein = Proteinseq(*seqInfo)
            if protein.Gene_name == ">DREB2A_P":
                print("UniprotID of DREB2A sequence: ", protein.uniprotID)
                print("Reviewed  status of DREB2A sequence: ", protein.reviewedS)
                print("Type of DREB2A sequence: ", protein.Sequence_type)
                print("Amino acid composition of DREB2B: ",protein.get_Character_Count(seqInfo[-1]) )
                print("Hydrophobicity of DREB2A: ", protein.hydrophobicity)

    # get the number of sequences created
    print("Sequence count: ",Sequence.Sequence_count)
























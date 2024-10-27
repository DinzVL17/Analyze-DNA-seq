# creating a Sequence class to store any biological  sequence
class Sequence:
    #define class variables
    Sequence_count = 0

    #creating constructor method
    def __init__(self, Gene_ID, Gene_name, Species_name, Subspecies_name,sequence):
        self.Gene_ID = Gene_ID
        self.Gene_name = Gene_name
        self.Sequence_type = Sequence.get_Seq_Type(Sequence,sequence)
        self.Sequence_length = len(sequence)
        self.Species_name = Species_name
        self.Subspecies_name = Subspecies_name
        Sequence.Sequence_count += 1

    #create a static method to split multiple FASTA sequences in a single text file and return a dictionary
    @staticmethod
    def fasta_Split(file):
        seqDiction = {}
        headerlist = []
        seq1 =""
        header =""
        #open the file
        with open(file, 'r') as file1:
            for line in file1:
                lines = line.strip()
                if lines != '\n':
                    if ">" in lines:
                        seq1 = ""
                        #create a list containing hyphen-separated fields in the header
                        headerlist = lines.split("-")
                        #store the gene name in header variable
                        header = headerlist[0]
                    else:
                        #store the sequence in seq1 variable
                        seq1 += lines
                        #create the dictionary
                        seqDiction[header] = headerlist + [seq1]

        return(seqDiction)

    #create an instance method to check the sequence type
    def get_Seq_Type(self,seq):
        #create a list of amino acids
        aa = ["R", "N", "D", "B", "E", "Q", "Z", "H", "I", "L", "K", "M", "F", "P", "S", "W", "Y", "V"]
        for base in seq:
            #check whether characters in sequence are equal to amino acids in list
            if base in aa:
                return ("amino acid sequence")
            elif base not in aa and "U" in seq:
                return ("mRNA sequence")
            elif base not in aa and "U" not in seq:
                return ("DNA sequence")

    #create an instance method to get the count of each character in sequence
    def get_Character_Count(self,seq):
        characterCount = {}
        for base in seq:
            if base in characterCount:
                characterCount[base] += 1
            else:
                characterCount[base] = 1
        return (characterCount)


if __name__ == "__main__":
    # create the dictionary
    output = (Sequence.fasta_Split("OSDREB_CDS.FASTA"))
    print(output)
    order = [1, 0, 2, 3,4]

    for values in output:
        # get the values of the dictionary and store in seqInfo list
        seqInfo = output[values]
        # reorder the list
        seqInfo = [seqInfo[i] for i in order]

        # create objects using list of parameters (seqInfo)
        object = Sequence(*seqInfo)
        # get the details(Gene ID, sequence length, and sequence type) of DREB1A DNA sequence
        if object.Gene_name == ">DREB1A_CDS":
            print(object.Gene_ID)
            print(object.Sequence_length)
            print(object.Sequence_type)

            # get the count of the four bases of the DREB1A coding sequence
            print(object.get_Character_Count(seqInfo[4]))
















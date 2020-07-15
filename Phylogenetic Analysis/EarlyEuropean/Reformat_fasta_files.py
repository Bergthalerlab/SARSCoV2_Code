#def main():

import re

# open the sequnces file
sequences = open("./data/gisaid_sequences.fasta", "r")
editedsequences = open("./data/gisaid_sequences_edited.fasta", "w")


# reformat the fasta names so that they fit to the metadata entry names
regex1 = "hCoV-19/"
regex2 = "\|.*\|\d+.*"
regex3 = " "

for sequence in sequences:
    newline = re.sub(regex1, "", sequence)
    newline = re.sub(regex2, "", newline)
    newline = re.sub(regex3, "", newline)
    editedsequences.write(newline)

editedsequences.close()
editedsequences = open("./data/gisaid_sequences_edited.fasta", "r")
sortedsequences = open("./data/gisaid_sequences_edited_sorted.fasta", "w")

# Sort duplicate sequences and CeMM sequences out
sequence_names = []
current_sequence = "Sequence Names"

for sequence in editedsequences:
    #if fasta name
    if sequence[0] == ">":
        # and add the fasta name to the list to block it in the future to avoid duplicates
        sequence_names.append(current_sequence)
        current_sequence = sequence
        
        #remove if CeMM sequence
        match = re.search("CeMM", current_sequence)
        if match:
            sequence_names.append(current_sequence)
        
        #check if it exists already in the list
        if str(current_sequence) not in sequence_names:
            sortedsequences.write(sequence)
    
    #if fasta sequence
    if sequence[0] != ">":
        #check if it exists already in the list
        if str(current_sequence) not in sequence_names:
            # if not add the sequence to the list
            sortedsequences.write(sequence)

#for sequence in editedsequences:
#    if sequence[0] == ">":
#        if str(current_sequence) not in sequence_names:
#            print(current_sequence)
#            sequence_names.append(current_sequence.rstrip())
#        current_sequence = sequence
#    if str(current_sequence) not in sequence_names:
#        sortedsequences.write(sequence)

editedsequences.close()
sequences.close()
sortedsequences.close()


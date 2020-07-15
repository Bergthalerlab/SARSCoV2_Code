#def main():

import re
import random

# open the sequnces file
filtered = open("./results/filtered.fasta", "r")
filteredsequencenames = open("./results/filterednames.txt", "w")
subsamplednames = open("./results/filteredsubsamplednames.txt", "w")
subsampled = open("./data/filteredsubsampled.fasta", "w")


# create a list of all filtered sequence names
namelist = []

for name in filtered:
    if name[0] == ">":
        namelist.append(name)
        filteredsequencenames.write(name)

filtered.close()
print(len(namelist), " sequences in the filtered file")

#select randomly 8000 sequences and all Austrian sequences

numsequences = []

for sequencename in namelist:
    match = re.search("Austria", sequencename)
    if match:
        subsamplednames.write(sequencename)
        numsequences.append(sequencename)

for sequencename in namelist:
    match = re.search("BavPat", sequencename)
    if match:
        subsamplednames.write(sequencename)
        numsequences.append(sequencename)

print(len(numsequences), " Austrian sequences + BavPat + Wuhan-Hu-1")

chooselist = namelist

while len(numsequences) < 8000:
    randomnumber = random.randint(0, len(chooselist)-1)
    
    # check if already existis in list and if not add it
    if chooselist[randomnumber] not in numsequences:
        subsamplednames.write(chooselist[randomnumber])
        numsequences.append(chooselist.pop(randomnumber))

print(len(numsequences), "sequences selected in total")


#copy sequences of the names into subsample file
filtered = open("./results/filtered.fasta", "r")
subsampledsequences = []

for entryline in filtered:
    
    if entryline[0] == ">":
        currentsequence = entryline
        if currentsequence in numsequences:
            subsampledsequences.append(currentsequence)
    
    if currentsequence in numsequences:
        subsampled.write(entryline)

subsampled.close()

# to be sure count the Austrian sequences in subsampling file
subsampled = open("./data/filteredsubsampled.fasta", "r")
numaustrians = []

for sequences in subsampled:
    match = re.search("Austria", sequences)
    if match:
        numaustrians.append(sequences)

print(len(subsampledsequences), "sequences in total in filtered subsampled file in /data")
print(len(numaustrians), "Austrian sequences in final subsampling file")


filtered.close()
filteredsequencenames.close()
subsamplednames.close()



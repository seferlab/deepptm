import os
import sys
import xlrd
import pandas as pd
import random
import numpy as np

def preparefasta(modifications,species_list, outfolder, left):
    """ Prepares Fasta file
    """
    for modification in modifications:
        fname = "{0}.txt".format(modification)
        for species in species_list:
            positives = set()
            negatives = set()
            with open(fname, "r") as infile:
                for line in infile:
                    splitted = line.split("\t")
                    assert(splitted[3] == modification)
                    if splitted[5] != species:
                        continue
                    
                    location = int(splitted[2])
                    sequence = splitted[6]
                    assert(sequence[location-1] == "K")
                
                    positions = []
                    for ind,char in enumerate(sequence):
                        if char == "K":
                            positions.append(ind+1)
                    
                    sequence = "O" * left + sequence + "O" * left
                    localcount = 0
                    for position in positions:
                        part = sequence[position-1:position+2*left]
                        if "O" in part:
                            continue
                        
                        assert(sequence[left+position-1] == "K")
                        if position == location:
                            localcount += 1
                            positives.add(part)
                        else:
                            negatives.add(part)
                    assert(localcount <= 1)
                    
            print(species, modification)
            print("Positives:", len(positives))
            negatives = negatives - positives # Clean overlaps
            print("Negatives:", len(negatives))

            # CSV Output
            outfpath = "{0}/{1}_{2}.csv".format(outfolder,modification,species_map[species])
            with open(outfpath,"w") as outfile:
                for positive in positives:
                    outfile.write("{0},positive\n".format(positive))
                for negative in negatives:
                    outfile.write("{0},negative\n".format(negative))

            # FASTA Output
            outfpath = "{0}/{1}_{2}.fasta".format(outfolder,modification,species_map[species])
            count = 1
            with open(outfpath,"w") as outfile:
                for positive in positives:
                    outfile.write(">seq{0}\n{1}\n".format(count,positive))
                    count += 1
                for negative in negatives:
                    outfile.write(">seq{0}\n{1}\n".format(count,negative))
                    count += 1
    
 
modifications = [
"Ubiquitination",
"Acetylation",
"Crotonylation",
"Glutarylation",
"Glycation",
"Succinylation",
]

species_list = [
    "Homo sapiens",
    "Mus musculus",
    "Saccharomyces cerevisiae (strain ATCC 204508 / S288c)"
]

species_map = {
    "Homo sapiens": "human",
    "Mus musculus": "mouse",
    "Saccharomyces cerevisiae (strain ATCC 204508 / S288c)": "yeast"
}

reverse_species_map = {value: key for key, value in species_map.items()}
    
left = 13
outfolder = "processed_data/{0}".format(2*left+1)
if not os.path.exists(outfolder):
    os.makedirs(outfolder)

#preparefasta(modifications,species_list, outfolder, left)
    
if left == 10:
    fname2info = {
        "1674144331.fas.1": ("human", "Ubiquitination", "30"),
        "1674150599.fas.1": ("human", "Crotonylation", "30"),
        "1674151139.fas.1": ("mouse", "Ubiquitination", "30"),
        "1674151906.fas.1": ("yeast", "Ubiquitination", "30"),
        "1674152320.fas.1": ("human", "Acetylation", "30"),
        "1674153844.fas.1": ("human", "Succinylation", "30"),

        "1674206343.fas.1": ("human", "Ubiquitination", "40"),
        "1674210803.fas.1": ("human", "Crotonylation", "40"),
        "1674212109.fas.1": ("mouse", "Ubiquitination", "40"),
        "1674220662.fas.1": ("yeast", "Ubiquitination", "40"),
        "1674375517.fas.1": ("human", "Succinylation", "40"),
        "1674375863.fas.1": ("human", "Glycation", "40"),
        "1674394117.fas.1": ("mouse", "Succinylation", "40"),
        "1674394291.fas.1": ("yeast", "Succinylation", "40"),

        "Crotonylation_human100.fasta": ("human", "Crotonylation", "100"),
        "Crotonylation_mouse100.fasta": ("mouse", "Crotonylation", "100"),
        "Crotonylation_yeast100.fasta": ("yeast", "Crotonylation", "100"),
        "Glycation_human100.fasta": ("human", "Glycation", "100"),
        "Glycation_mouse100.fasta": ("mouse", "Glycation", "100"),
        "Glycation_yeast100.fasta": ("yeast", "Glycation", "100"),
        "Succinylation_human100.fasta": ("human", "Succinylation", "100"),
        "Succinylation_mouse100.fasta": ("mouse", "Succinylation", "100"),
        "Succinylation_yeast100.fasta": ("yeast", "Succinylation", "100"),
        "Ubiquitination_mouse100.fasta": ("mouse", "Ubiquitination", "100"),
        "Ubiquitination_human100.fasta": ("human", "Ubiquitination", "100"),
        "Ubiquitination_yeast100.fasta": ("yeast", "Ubiquitination", "100"),
    }
elif left == 16:
    fname2info = {
        "succhuman04.fas.1": ("human", "Succinylation", "40"),
    }
elif left == 7:
    fname2info = {
        "1700463771.fas.1": ("human", "Succinylation", "40"),
    }
elif left == 13:
    fname2info = {
        "1700464031.fas.1": ("human", "Succinylation", "40"),
    }
        
        
for fname, tags in fname2info.items():
    print(fname, tags)
    species, modification, threshold = tags

    localcount = 0
    positives = set()
    negatives = set()
    seq2label = {}
    mainfname = "{0}.txt".format(modification)
    print(mainfname)
    with open(mainfname,"r") as infile:
        for line in infile:
            splitted = line.split("\t")
            assert(splitted[3] == modification)
            if splitted[5] != reverse_species_map[species]:
                continue
                
            location = int(splitted[2])
            sequence = splitted[6]
            assert(sequence[location-1] == "K")
                
            positions = []
            for ind,char in enumerate(sequence):
                if char == "K":
                    positions.append(ind+1)

            sequence = "O" * left + sequence + "O" * left
            localcount = 0
            for position in positions:
                part = sequence[position-1:position+2*left]
                if "O" in part:
                    continue

                assert(sequence[left+position-1] == "K")
                if position == location:
                    localcount += 1
                    positives.add(part)
                else:
                    negatives.add(part)
                assert(localcount <= 1)

    print("Positives:", len(positives))
    negatives = negatives - positives # Clean overlaps
    print("Negatives:", len(negatives))

    for seq in positives:
        seq2label[seq] = True
    for seq in negatives:
        seq2label[seq] = False
    
    outfpath = "{0}/cdhit/{1}".format(outfolder,fname)
    positives, negatives = set(), set()
    with open(outfpath,"r") as infile:
        for line in infile:
            seq = line.rstrip()
            if ">seq" in seq:
                continue

            if seq2label[seq]:
                positives.add(seq)
            else:
                negatives.add(seq)    

    print("After thresholding: ")
    print("Positives: ", len(positives))
    print("Negatives: ", len(negatives))

    outfpath2 = "{0}/output/{1}_{2}_cdhit{3}.csv".format(outfolder,modification,species,threshold)
    if os.path.exists(outfpath2):
        continue
    
    with open(outfpath2,"w") as outfile:
        for positive in positives:
            outfile.write("{0},positive\n".format(positive))
        for negative in negatives:
            outfile.write("{0},negative\n".format(negative))

    positives = list(positives)
    negatives = list(negatives)
    if len(positives) < len(negatives):
        negatives = random.sample(negatives, len(positives))
    elif len(positives) > len(negatives):
        positives = random.sample(positives, len(negatives))

    print("After sampling: ")
    print("Positives: ", len(positives))
    print("Negatives: ", len(negatives))
        
    equalfpath = "{0}/output/{1}_{2}_cdhit{3}_equal.csv".format(outfolder,modification,species,threshold)
    with open(equalfpath,"w") as outfile:
        for positive in positives:
            outfile.write("{0},positive\n".format(positive))
        for negative in negatives:
            outfile.write("{0},negative\n".format(negative))

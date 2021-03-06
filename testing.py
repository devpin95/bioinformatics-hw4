import pickle
import json
from model import Model
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("seq_file", help="a file containing a DNA sequence")
parser.add_argument("test_seqs", help="a file containing a list of test sequences")
args = parser.parse_args()

full_seq_file_name = args.seq_file
test_filename = args.test_seqs

# full_seq_file_name = "sequenceOfEcoliStrainM54.txt"
full_seq_file = open(full_seq_file_name)

full_seq = full_seq_file.read().split()  # read the file, and remove whitespace and newline chars
full_seq = ''.join(full_seq)  # rejoin the split with no space between elements

# test_filename = "test_seqs.txt"
test_seqs = []

testfile = open(test_filename)
for line in testfile:
    metadata = {}

    seq = line.split()
    metadata['query-start'] = int(seq[0]) - 1
    metadata['query-end'] = int(seq[1]) - 1
    metadata['targets'] = []

    for i in range(2, len(seq) - 1, 2):
        rangestr = seq[i] + seq[i + 1]
        rangestr = rangestr[1:len(rangestr)-1].split(',')
        metadata['targets'].append({'start': int(rangestr[0]) - 1, 'end': int(rangestr[1]) - 1})

    test_seqs.append(metadata)

model = Model(full_seq, {})

bmodelfile = 'model_bin'
model.load_model(bmodelfile)

# model.test(test_seqs[0]['query-start'], test_seqs[0]['query-end'])

for seq in test_seqs:
    model.test(seq['query-start'], seq['query-end'])

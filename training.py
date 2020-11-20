import re
from model import Model

full_seq_file_name = "sequenceOfEcoliStrainM54.txt"
training_seqs_file_name = "train_seqs.txt"

full_seq_file = open(full_seq_file_name)

full_seq = full_seq_file.read().split()  # read the file, and remove whitespace and newline chars
full_seq = ''.join(full_seq)  # rejoin the split with no space between elements

training_seqs_file = open(training_seqs_file_name)

sequence_metadata = []

for line in training_seqs_file:
    metadata = {}

    seq = line.split()
    seq = ' '.join(seq)
    seq = seq.split(' ')

    metadata['query_start'] = int(seq[0]) - 1  # subtract 1 because sequence starts at 1
    metadata['query_end'] = int(seq[1]) - 1  # subtract 1 because sequence starts at 1
    metadata['gene_coding_runs'] = []

    coding_runs = ' '.join(seq[2:])
    coding_runs = re.split('] |,', coding_runs)

    for i in range(0, len(coding_runs), 2):
        coding_runs[i] = coding_runs[i].replace('[', '')  # remove [ from the string
        coding_runs[i+1] = coding_runs[i+1].replace(']', '')  # remove ] from the string

        # add the start and end index of the coding region, subtracting 1 because the sequence count starts from 1
        metadata['gene_coding_runs'].append({'start': int(coding_runs[i]) - 1, 'end': int(coding_runs[i+1]) - 1})

    # append this gene metadata
    sequence_metadata.append(metadata)


model = Model(full_seq, sequence_metadata)
model.train()
model.save_model_json('model.json')
model.save_model_bin('model_bin')

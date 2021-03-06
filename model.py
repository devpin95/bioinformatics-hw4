import pickle
import json
import math
from pandas import *

NONCODING = 'noncoding'
START_CODON = 'start-codon'
STOP_CODONS = 'stop-codons'
INTERNAL_CODONS = 'internal-codons'
SUB_MODEL = 'submodel'
OCC = 'occurrences'
START = 'start'

state_transition_filename = "model_transitions.json"


class Model:
    model = {
        'start': {'a': 1, 't': 1, 'g': 1, 'c': 1},
        'noncoding': {
            'a': {'a': 0, 't': 0, 'g': 0, 'c': 0, 'transitions_to_end': 1, 'transitions_to_start_codon': {'a': 1, 't': 1, 'g': 1, 'c': 1}, 'occurrences': 0},
            't': {'a': 0, 't': 0, 'g': 0, 'c': 0, 'transitions_to_end': 1, 'transitions_to_start_codon': {'a': 1, 't': 1, 'g': 1, 'c': 1}, 'occurrences': 0},
            'g': {'a': 0, 't': 0, 'g': 0, 'c': 0, 'transitions_to_end': 1, 'transitions_to_start_codon': {'a': 1, 't': 1, 'g': 1, 'c': 1}, 'occurrences': 0},
            'c': {'a': 0, 't': 0, 'g': 0, 'c': 0, 'transitions_to_end': 1, 'transitions_to_start_codon': {'a': 1, 't': 1, 'g': 1, 'c': 1}, 'occurrences': 0}
        },
        'start-codon': {
            'first-nt': {'transitions_to_second_nt': 1.0},
            'second-nt': {'transitions_to_third_nt': 1.0},
            'third-nt': {'transition_count': 0, 'transition_codons': {}}
        },
        'internal-codons': {
            'stop-codons': {
                'taa': {'a': 1, 't': 1, 'g': 1, 'c': 1},
                'tag': {'a': 1, 't': 1, 'g': 1, 'c': 1},
                'tga': {'a': 1, 't': 1, 'g': 1, 'c': 1}
            }
        }
    }

    nt_list = ['a', 'c', 't', 'g']
    num_nt = 4

    def __init__(self, seq, metadata):
        self.seq = seq
        self.metadata = metadata
        self.query_count = len(metadata)
        self.state_transitions = None
        self.T1 = []
        self.T2 = []
        # self.state_order = ['noncoding', 'start-codon-first', 'start-codon-t', 'start-codon-g', 'internal-codons', 'stop-codons']

    def load_model(self, bfilename):
        with open(bfilename, 'rb') as model_file:
            self.model = pickle.load(model_file)

    def train(self):
        print('Training model...')
        self.train_start_transition()
        self.train_noncoding_region()
        self.train_start_codon_region()
        self.train_internal_codons()
        self.train_stop_codon_transitions()

        self.save_model_json('transition_counts.json')

        self.convert_model_to_prob()

    def convert_model_to_prob(self):
        self.start_to_prob()
        self.noncoding_to_prob()
        self.start_codon_to_prob()
        self.internal_codons_to_prob()
        self.end_codon_to_prob()

    def train_start_transition(self):
        for i in range(0, self.query_count):
            query = self.metadata[i]
            query_start_index = query['query_start']
            start_char = self.seq[query_start_index]
            self.model[START][start_char] += 1

    def train_noncoding_region(self):
        # train_noncoding_region
        # trains the transitions between single characters in the noncoding region of the training sequences
        # Keeps track of which characters ended the query range so that we can calculate the probability a character
        # moves to the end state of the model

        for i in range(0, self.query_count):
            mid_runs, end_runs = self.build_noncoding_seq(i)

            # iterate through the intergenic regions that do no end at the end of the query range
            for seq in mid_runs:
                # iterate through the sequence and count the number of transitions from one letter to the next
                for j in range(0, len(seq) - 1):
                    current_nt = seq[j]  # the current base
                    next_nt = seq[j+1]  # the next base in the sequence
                    self.model[NONCODING][current_nt][next_nt] += 1  # increment the number of transitions to the next base
                    self.model[NONCODING][current_nt][OCC] += 1  # increment the number of times we've seem this base

            # iterate through the intergenic regions that run up to the end of the query range
            for seq in end_runs:
                # iterate through the sequence and count the number of transitions from one letter to the next
                for j in range(0, len(seq) - 1):
                    current_nt = seq[j]  # the current base
                    next_nt = seq[j+1]  # the next base in the sequence
                    self.model[NONCODING][current_nt][next_nt] += 1  # increment the number of transitions to the next base
                    self.model[NONCODING][current_nt][OCC] += 1  # increment the number of times we've seem this base

                    # we are looking at the last pair of bases in the sequence
                    # count the number of times that the last character ended the sequence
                    if j == len(seq) - 2:
                        self.model[NONCODING][next_nt]['transitions_to_end'] += 1

        self.train_noncoding_transition_to_start_codon()

    def train_noncoding_transition_to_start_codon(self):
        # train_noncoding_transition_to_start_codon
        # trains the transitions from the last base in the noncoding region to the first base in the state codon

        # loop through the query reanges
        for i in range(0, self.query_count):

            # get the ranges of the coding regions
            gene_ranges = self.metadata[i]['gene_coding_runs']

            # go through each of the coding regions for this query
            for bound in gene_ranges:
                start = bound['start']  # the index of the first base in the coding region
                coding_start = self.seq[start]  # the first character in the coding region (first char of start codon)
                noncoding_end = self.seq[start-1]  # the last character in the noncoding region

                # increment the transition from the noncoding base to the first coding base
                self.model[NONCODING][noncoding_end]['transitions_to_start_codon'][coding_start] += 1

    def train_start_codon_region(self):
        # train_start_codon_region
        # trains the transition in the start codon
        # we already know that the first character always transitions to 't' and that the second character always
        # transitions to 'g'. We need to train the transitions from the last character to the first codon after
        # the start codon. Codons will be considered as groups of 3

        # generate possible list of codons the last nucleotide in the start codon can go to
        for i in range(0, self.num_nt):
            for j in range(0, self.num_nt):
                for k in range(0, self.num_nt):
                    codon = self.nt_list[i] + self.nt_list[j] + self.nt_list[k]
                    self.model[START_CODON]['third-nt']['transition_codons'][codon] = 1

        # iterate through the queries
        for i in range(0, self.query_count):

            # get the ranges of the coding regions in each query
            ranges = self.metadata[i]['gene_coding_runs']

            # count the number of times the last nucleotide in the start codon transitions to each possible codon
            for coding_range in ranges:
                start = coding_range['start'] + 3  # the index of the first character after the start codon
                codon = self.seq[start: start + 3]  # concatenate the next 3 characters

                # increment the transition from 'g' in the start codon to the codon that appeared after it in this seq
                # also increment the number of transitions we've seen so that we can later find the probability easier
                self.model[START_CODON]['third-nt']['transition_codons'][codon] += 1
                self.model[START_CODON]['third-nt']['transition_count'] += 1

    def train_internal_codons(self):
        # train_internal_codons
        # Train the transitions between codons in the coding region, including the stop codons

        # first, build the tree, connect every possible codon to all other codons
        self.build_internal_codon_tree()

        # iterate through the queries
        for i in range(0, self.query_count):
            # get the list of codons in each coding region
            codon_lists = self.build_internal_codons(self.metadata[i]['gene_coding_runs'])

            # loop through each list of codons
            for clist in codon_lists:
                # loop through each codon and see what comes after it, do not look at the last codon
                for j in range(0, len(clist) - 1):
                    current_codon = clist[j]  # current codon
                    next_codon = clist[j+1]  # codon to transition to

                    # increment the codon the current codon transitioned to
                    self.model[INTERNAL_CODONS][current_codon][next_codon] += 1

    def train_stop_codon_transitions(self):
        # train_stop_codon_transitions
        # train what letter to return to in the noncoding region after a stop codon is encountered

        # iterate through the queries
        for i in range(0, self.query_count):
            # get the coding ranges in each query
            codon_lists = self.metadata[i]['gene_coding_runs']

            # loop through the coding ranges for this query
            for clist in codon_lists:
                end = clist['end']  # index of the last character in the coding region

                # concatenate the last 3 chars of the seq (this is the stop codon)
                # and find what it transitions to back in the noncoding region
                end_codon = self.seq[end-2] + self.seq[end-1] + self.seq[end]
                next_base = self.seq[end+1]

                # increment the base the stop codon transitioned to
                self.model['internal-codons']['stop-codons'][end_codon][next_base] += 1

    def build_noncoding_seq(self, query_index):
        # build_noncoding_seq: query to examine
        # for a given query, concatenate all continuous runs in the noncoding region
        # for example, 1 2 3 4 5 a b c 6 7 8 e f g 9 h i is a sequence
        # number represent noncoding regions and letters coding regions
        # we need to return a list of sequences in the form:
        # [ '1 2 3 4 5', '6 7 8', '9' ]

        mid_seqs = []  # sequences that are not at the end of the sequence after the last coding region
        end_seqs = []  # sequences that are after the last coding region
        gene = self.metadata[query_index]  # collect the metadata for this query

        # concatenate the sequence from the start of the gene up to the first coding region
        lower_bound = gene['query_start']
        upper_bound = gene['gene_coding_runs'][0]['start']
        mid_seqs.append(self.seq[lower_bound:upper_bound])

        # get the number of coding regions that we need to split up
        num_runs = len(gene['gene_coding_runs'])

        # loop through the coding regions so that we can extract the noncoding regions between them
        for i in range(0, num_runs - 1):
            lower_bound = gene['gene_coding_runs'][i]['end'] + 1
            upper_bound = gene['gene_coding_runs'][i+1]['start']
            mid_seqs.append(self.seq[lower_bound:upper_bound])

        # sequence from the last coding region to the end of the gene
        lower_bound = gene['gene_coding_runs'][-1]['end'] + 1
        upper_bound = gene['query_end']

        # only append the seq if it is nonempty
        if lower_bound < upper_bound:
            end_seqs.append(self.seq[lower_bound:upper_bound])

        return mid_seqs, end_seqs

    def build_internal_codon_tree(self):
        # build_internal_codon_tree
        # forms all combinations of the bases 'a', 'c', 't', 'g' in groups of 3

        for i in range(0, self.num_nt):
            for j in range(0, self.num_nt):
                for k in range(0, self.num_nt):
                    codon = self.nt_list[i] + self.nt_list[j] + self.nt_list[k]
                    self.model[INTERNAL_CODONS][codon] = {}

                    for u in range(0, self.num_nt):
                        for m in range(0, self.num_nt):
                            for n in range(0, self.num_nt):
                                subcodon = self.nt_list[u] + self.nt_list[m] + self.nt_list[n]
                                self.model[INTERNAL_CODONS][codon][subcodon] = 1

    def build_internal_codons(self, ranges):
        # build_internal_codons: ranges
        # takes a list of ranges with start and end indexes
        # returns a list of all the codons in the range
        # codons are formed in groups of 3 beginning at the start index

        codon_seqs = []  # list to hold extracted codons from each range

        # loop through the ranges
        for coding_range in ranges:
            start = coding_range['start'] + 3  # index of the first character to consider, skipping the start codon
            end = coding_range['end']  # end of the coding range

            # concatenate the characters that we need to extract codons from
            # length of string will be a multiple of 3
            nt_seq = self.seq[start:end+1]

            codon_seq = []  # list to hold the codons for this coding range

            # loop through the concatenated string 3 at a time
            for i in range(0, len(nt_seq), 3):
                # concatenate groups of three to form the codon and add it to the list of codons extracted
                codon_seq.append(nt_seq[i] + nt_seq[i+1] + nt_seq[i+2])

            # add the codons for this list to the list of codons extracted from each range
            codon_seqs.append(codon_seq)

        return codon_seqs

    def start_to_prob(self):
        bsum = sum(self.model[START].values())

        for base in self.model[START]:
            self.model[START][base] /= bsum
            self.model[START][base] = math.log(self.model[START][base], 2)

    def internal_codons_to_prob(self):
        # internal_codons_to_prob
        # takes the raw codon counts and converts them to transition probabilities
        # probabilities are the number of times a codon is seen divided by the total number of codons

        # loop through the codon states
        for codon in self.model[INTERNAL_CODONS]:
            # do no look at the stop-codons def in the dict
            if codon != 'stop-codons':
                # sum the number of codons seen in the sequence
                codon_count = sum(self.model[INTERNAL_CODONS][codon].values())

                # for each of the condons we can transition to, divide by the total number of codons seen
                for trans_codon in self.model[INTERNAL_CODONS][codon]:
                    count = self.model[INTERNAL_CODONS][codon][trans_codon]
                    prob = count / codon_count
                    self.model[INTERNAL_CODONS][codon][trans_codon] = math.log(prob, 2)

    def noncoding_to_prob(self):
        # noncoding_to_prob
        # converts raw counts of noncoding characters to transition probabilites
        # for chars A, B, C, a transition probability from char A to char B is the number of times A transitioned to
        # B divided by the number of all transitions from A to A or B or C

        # iterate through each possible nucleotide
        for nucleotide in self.model[NONCODING]:
            # count the number of times the nucleotide has been seen
            count = self.model[NONCODING][nucleotide][OCC]

            # loop through each nucleotide the current nucleotide could transition to
            for transition_nucleotide in self.model[NONCODING][nucleotide]:
                # skip the other model metadata except the nucleotide metadata
                if transition_nucleotide in self.nt_list:
                    # divide each the total number of transitions by the count, include the transitions to the first
                    # char in the start codon
                    prob = self.model[NONCODING][nucleotide][transition_nucleotide] / count
                    self.model[NONCODING][nucleotide][transition_nucleotide] = math.log(prob, 2)
                    # self.model[NONCODING][nucleotide][transition_nucleotide] /= count
                    prob = self.model[NONCODING][nucleotide]['transitions_to_start_codon'][transition_nucleotide] / count
                    # self.model[NONCODING][nucleotide]['transitions_to_start_codon'][transition_nucleotide] /= count
                    self.model[NONCODING][nucleotide]['transitions_to_start_codon'][transition_nucleotide] = math.log(prob,2)

            # we also need to know how ofter this nucleotide ended the query
            self.model[NONCODING][nucleotide]['transitions_to_end'] /= count

    def start_codon_to_prob(self):
        # start_codon_to_prob
        # finds the transition probabilities from the last char of the start codon to the second codon

        # find the total number of transitions we've seen
        count = sum(self.model[START_CODON]['third-nt']['transition_codons'].values())

        # loop through each of the possible codons to transition to and divide by the count
        for codon in self.model[START_CODON]['third-nt']['transition_codons']:
            prob = self.model[START_CODON]['third-nt']['transition_codons'][codon] / count
            # self.model[START_CODON]['third-nt']['transition_codons'][codon] /= count
            self.model[START_CODON]['third-nt']['transition_codons'][codon] = math.log(prob, 2)

    def end_codon_to_prob(self):
        # end_codon_to_prob
        # converts raw counts of transitions from the stop codon to the first char in the noncoding region

        # loop through the possible stop codons
        for stop_codon in self.model['internal-codons']['stop-codons']:
            # for this stop codon, sum the number of transitions we've counted for each nucleotide a, c, t, and g
            count = sum(self.model['internal-codons']['stop-codons'][stop_codon].values())

            # for each of the nucleotides, divide the counted transitions by the total number of transitions
            for base_trans in self.model['internal-codons']['stop-codons'][stop_codon]:
                prob = self.model['internal-codons']['stop-codons'][stop_codon][base_trans] / count
                self.model['internal-codons']['stop-codons'][stop_codon][base_trans] = math.log(prob, 2)

    def save_model_bin(self, filename):
        # Save the model to binary file for use in prediction
        file = open(filename, 'wb')
        pickle.dump(self.model, file)

    def save_model_json(self, filename):
        # Save the model to json for debugging
        j = json.dumps(self.model, indent=4, sort_keys=True)
        f = open(filename, 'w')
        f.write(j)
        f.close()

    def test(self, start, end):
        if self.state_transitions is None:
            self.load_state_transitions()
        test_seq = self.seq[start:end]

        print(str(start) + "-" + str(end) + ": ", end='')

        # test_seq = 'gcatgacctaatagaat'

        self.build_tables(test_seq)

        self.viterbi(test_seq)

        # print(self.T1[7][len(test_seq)-1])
        # print(self.T2[7][len(test_seq) - 1])

        # print(DataFrame(self.T1))

        # columns = [char for char in test_seq]
        # print(DataFrame(self.T2, columns=columns))

        path = []
        distance = 0
        move_to = self.state_transitions['metadata']['state-space'].index('noncoding')
        for i in range(len(test_seq) - 1, 0, -1):
            path.insert(0, self.T2[move_to][i - distance])
            move_to = self.T2[move_to][i - distance]
            if move_to == 0:
                break
            if move_to == 6 or move_to == 5:
                distance += 2

        self.extract_ranges(start, end, path)

        self.T1 = []
        self.T2 = []

    def load_state_transitions(self):
        with open(state_transition_filename) as st_json:
            self.state_transitions = json.load(st_json)

    def build_tables(self, seq):
        for i in range(0, len(self.state_transitions['metadata']['state-space'])):
            self.T1.append([])
            self.T2.append([])
            for j in range(0, len(seq)):
                self.T1[i].append(0)
                self.T2[i].append(-1)

    # def set_first_state(self, seq):
    #     state = self.state_transitions['metadata']['state-order'].index("noncoding")
    #     self.probability_table[state][0] = 0
    #     self.state_changes_table[state][0]['from'] = 'start'

    def viterbi(self, seq):
        # state space: ["start", "noncoding", "start-codon-first", "start-codon-t", "start-codon-g", "internal-codons", "stop-codons", "end"]
        s = self.state_transitions['metadata']['state-space']
        pi = [-math.inf, self.model['start'][seq[0]], -math.inf, -math.inf, -math.inf, -math.inf, -math.inf, -math.inf]
        pi_a = [-1, 0, -1, -1, -1, -1, -1, -1]

        for i in range(0, len(s)):
            self.T1[i][0] = pi[i]
            self.T2[i][0] = pi_a[i]

        for base in range(1, len(seq)):
            for state in range(0, len(s)):
                # print(base, s[state], end='')
                self.T1[state][base], self.T2[state][base] = self.max_and_arg(state, s, seq, base)

    def max_and_arg(self, state, s, seq, base):
        max_prob = -math.inf
        arg_max = -1
        current_base = seq[base]
        prev_base = seq[base - 1]
        current_codon = ''
        prev_codon = ''

        codon_range = False
        if base >= 3:
            codon_range = True
            prev_codon = seq[base - 3] + seq[base - 2] + seq[base - 1]

        if base <= len(seq) - 4:
            current_codon = seq[base] + seq[base + 1] + seq[base + 2]

        prob_list = [-math.inf, -math.inf, -math.inf, -math.inf, -math.inf, -math.inf, -math.inf, -math.inf, -math.inf]

        for prev_state in range(0, len(s)):
            current_state_name = s[state]
            prev_state_name = s[prev_state]

            prob = -math.inf

            if current_state_name == 'start':
                # prob = self.T1[prev_state][base - 1] + self.model[current_state_name][prev_base]
                break

            elif current_state_name == 'noncoding':
                int_prob, sc = self.v_noncoding(current_base, prev_base, current_state_name, prev_state_name, prev_codon)
                if sc:
                    prob = self.T1[prev_state][base - 3] + int_prob
                else:
                    prob = self.T1[prev_state][base - 1] + int_prob

            elif current_state_name == 'start-codon-first':
                if base + 2 < len(seq):
                    next_2_bases = seq[base + 1] + seq[base + 2]
                    if next_2_bases != 'tg':
                        break

                prob = self.T1[prev_state][base - 1] + self.v_start_codon_first(current_base, prev_base, current_state_name, prev_state_name)

            elif current_state_name == 'start-codon-t':
                if base + 1 < len(s):
                    next_base = seq[base + 1]
                    if next_base != 'g':
                        break

                prob = self.T1[prev_state][base - 1] + self.v_start_codon_t(current_base, prev_base, current_state_name, prev_state_name)

            elif current_state_name == 'start-codon-g':
                if prev_base != 't':
                    break
                prob = self.T1[prev_state][base - 1] + self.v_start_codon_g(current_base, prev_base, current_state_name, prev_state_name)

            elif current_state_name == 'internal-codons':
                prev_2_bases = ''
                if codon_range:
                    if base - 2 > 0:
                        prev_2_bases = seq[base - 2] + seq[base - 1]

                    if prev_state_name == 'internal-codons':
                        prob = self.T1[prev_state][base - 3] + self.v_internal_codons(current_codon, prev_2_bases,
                                                                                      prev_codon, current_state_name,
                                                                                      prev_state_name)
                    else:
                        prob = self.T1[prev_state][base - 1] + self.v_internal_codons(current_codon, prev_2_bases,
                                                                                      prev_codon, current_state_name,
                                                                                      prev_state_name)

            elif current_state_name == 'stop-codons':
                if codon_range:
                    prob = self.T1[prev_state][base - 3] + self.v_stop_codons(current_codon, prev_codon, current_state_name, prev_state_name)

            elif current_state_name == 'end':
                break

            prob_list[prev_state] = prob

            if prob > max_prob:
                max_prob = prob
                arg_max = prev_state

        # print(prob_list)
        return max_prob, arg_max

    def v_noncoding(self, current_base, prev_base, current_state_name, prev_state_name, prev_codon):
        # Calculate transition probabilities into the noncoding submodel in state current_base
        # State space: start, noncoding, stop codons

        prob = -math.inf  # store the probability
        from_stop_codon = False

        stop_codons = self.state_transitions['metadata']['stop-codons']

        # get the states we could transition from into the noncoding submodel
        transition_space = set([pstate['submodel'] for pstate in self.state_transitions[current_state_name][current_base]['from']])

        if prev_state_name in transition_space:

            if prev_state_name == 'start':
                prob = -math.inf  # store the probability

            elif prev_state_name == 'stop-codons':
                if prev_codon in stop_codons:
                    prob = self.model[INTERNAL_CODONS][STOP_CODONS][prev_codon][current_base]
                    # prob = 0.0
                    from_stop_codon = True

            elif prev_state_name == 'noncoding':
                prob = self.model[NONCODING][prev_base][current_base]

        return prob, from_stop_codon

    def v_start_codon_first(self, current_base, prev_base, current_state_name, prev_state_name):

        prob = -math.inf

        transition_space = set([d['submodel'] for d in self.state_transitions['start-codon-first'][current_base]['from']])

        if prev_state_name in transition_space:
            if prev_state_name == 'noncoding':
                prob = self.model[NONCODING][prev_base]['transitions_to_start_codon'][current_base]

        return prob

    def v_start_codon_t(self, current_base, prev_base, current_state_name, prev_state_name):
        # Calculate transition probabilities into the start-codon-t submodel in state current_base
        # State space: start-codon-first

        prob = -math.inf  # start with a log-probability of -infinity for a probability of 0.0

        if current_base == 't':
            if prev_state_name == 'start-codon-first':
                prob = math.log(1.0, 2)

        return prob

    def v_start_codon_g(self, current_base, prev_base, current_state_name, prev_state_name):
        # Calculate transition probabilities into the start-codon-g submodel in state current_base
        # State space: start-codon-t

        prob = -math.inf  # start with a log-probability of -infinity for a probability of 0.0

        if current_base == 'g' and prev_base == 't':
            if prev_state_name == 'start-codon-t':
                prob = math.log(1.0, 2)

        return prob

    def v_internal_codons(self, current_codon, prev_base, prev_codon, current_state_name, prev_state_name):
        # Calculate transition probabilities into the internal-codon submodel in state current_base
        # State space: start, internal-codons, start-codon-g

        prob = -math.inf  # start with a log-probability of -infinity for a probability of 0.0

        stop_codons = self.state_transitions['metadata']['stop-codons']

        # get the states we could transition from into the noncoding submodel
        transition_space = set([pstate['submodel'] for pstate in
                            self.state_transitions[current_state_name]['from']])

        if prev_state_name in transition_space:
            if prev_state_name == 'internal-codons' and prev_codon not in stop_codons:
                if current_codon != '':
                    prob = self.model[INTERNAL_CODONS][prev_codon][current_codon]

            elif prev_state_name == 'start-codon-g':
                if current_codon != '' and prev_base == 'tg':
                    prob = self.model['start-codon']['third-nt']['transition_codons'][current_codon]

        return prob

    def v_stop_codons(self, current_codon, prev_codon, current_state_name, prev_state_name):
        # Calculate transition probabilities into the internal-codon submodel in state current_base
        # State space: start, internal-codons

        prob = -math.inf  # start with a log-probability of -infinity for a probability of 0.0

        stop_codons = self.state_transitions['metadata']['stop-codons']

        if current_codon in stop_codons:
            submodel = 'internal-codons'

            if prev_state_name == 'internal-codons':
                prob = self.model[submodel][prev_codon][current_codon]

        return prob

    def v_end(self, current_base, current_state_name, prev_state_name):
        prob = -math.inf

        if prev_state_name == 'noncoding':
            prob = self.model[NONCODING][current_base]['transitions_to_end']

        return prob

    def extract_ranges(self, start, end, path):
        normalized_path = []

        for i in range(0, len(path)):
            if path[i] == 5 or path[i] == 6:
                normalized_path.append(path[i])
                normalized_path.append(path[i])
                normalized_path.append(path[i])
            else:
                normalized_path.append(path[i])

        starts = []
        ends = []
        for i in range(0, len(normalized_path)):
            if normalized_path[i] == 2:
                starts.append(i)

        for starti in starts:
            for i in range(starti, len(normalized_path)):
                if normalized_path[i] == 6:
                    ends.append(i)
                    break

        if len(starts) == 0:
            print("no gene", end='')

        for i in range(0, len(starts)):
            start_index = start + starts[i]
            end_index = start + ends[i]

            print("[" + str(start_index) + ", " + str(end_index) + "]", end=' ')

        print()

        # print("STARTS: ", starts, "ENDS:", ends)

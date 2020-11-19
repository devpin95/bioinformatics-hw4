# {'a': {'a': 27675, 't': 26974, 'g': 22421, 'c': 22232, 'occurrences': 99302},
# 't': {'a': 18248, 't': 30094, 'g': 29028, 'c': 24297, 'occurrences': 101667},
# 'g': {'a': 23661, 't': 23942, 'g': 25855, 'c': 36272, 'occurrences': 109730},
# 'c': {'a': 29753, 't': 20608, 'g': 32441, 'c': 25049, 'occurrences': 107851}}

NONCODING = 'noncoding'
START_CODON = 'start-codon'
INTERNAL_CODONS = 'internal-codons'
SUB_MODEL = 'submodel'
TRANS_TO = 'to'
OCC = 'occurrences'


class Model:
    model = {
        'noncoding': {
            'a': {'a': 0, 't': 0, 'g': 0, 'c': 0, 'occurrences': 0, 'transitions_to_start_codon': 1, 'transitions_to_end': 1},
            't': {'a': 0, 't': 0, 'g': 0, 'c': 0, 'occurrences': 0, 'transitions_to_start_codon': 1, 'transitions_to_end': 1},
            'g': {'a': 0, 't': 0, 'g': 0, 'c': 0, 'occurrences': 0, 'transitions_to_start_codon': 1, 'transitions_to_end': 1},
            'c': {'a': 0, 't': 0, 'g': 0, 'c': 0, 'occurrences': 0, 'transitions_to_start_codon': 1, 'transitions_to_end': 1}
        },
        'start-codon': {
            'first-nt': {'a': 0, 't': 0, 'g': 0, 'c': 0, 'transitions_to_internal_codons': 0.0},
            'second-nt': {'a': 0.0, 't': 1.0, 'g': 0.0, 'c': 0.0, 'transitions_to_internal_codons': 0.0},
            'third-nt': {'a': 0.0, 't': 0.0, 'g': 1.0, 'c': 0.0, 'transitions_to_internal_codons': 1.0, 'transition_count': 0}
        },
        'internal-codons': {
            'stop-codons': ['taa', 'tag', 'tga'],
        },
        'submodel': {
            'start': {
                'to': {
                    'noncoding': 1.0
                }
            },
            'noncoding': {
                'to': {
                    'start-codon': 0.5,
                    'end': 0.5
                }
            },
            'start-codon': {
                'to': {
                    'internal-codons': 1.0
                }
            },
            'internal-codons': {
                'to': {
                    'noncoding': 1.0
                }
            }
        }
    }

    def __init__(self, seq, metadata):
        self.seq = seq
        self.metadata = metadata
        self.query_count = len(metadata)

    def train(self):
        self.train_noncoding_region()
        print(self.model[NONCODING])

        self.train_start_codon_region()
        print(self.model[START_CODON])

        self.train_internal_codons()
        print(self.model[INTERNAL_CODONS])

    def train_noncoding_region(self):
        current_query = 0

        for i in range(0, self.query_count):
            mid_runs, end_runs = self.build_noncoding_seq(i)

            for seq in mid_runs:
                for j in range(0, len(seq) - 1):
                    current_nt = seq[j]
                    next_nt = seq[j+1]
                    self.model[NONCODING][current_nt][next_nt] += 1
                    self.model[NONCODING][current_nt][OCC] += 1

                    if j == len(seq) - 2:
                        self.model[NONCODING][next_nt]['transitions_to_start_codon'] += 1

            for seq in end_runs:
                for j in range(0, len(seq) - 1):
                    current_nt = seq[j]
                    next_nt = seq[j+1]
                    self.model[NONCODING][current_nt][next_nt] += 1
                    self.model[NONCODING][current_nt][OCC] += 1

                    if j == len(seq) - 2:
                        self.model[NONCODING][next_nt]['transitions_to_end'] += 1

        print('')

    def train_start_codon_region(self):
        for i in range(0, self.query_count):
            ranges = self.metadata[i]['gene_coding_runs']

            for coding_range in ranges:
                start = coding_range['start']
                self.model[START_CODON]['first-nt'][self.seq[start]] += 1

        # check if any nt is 0, if so add 1 to smooth the results
        for nt in self.model[START_CODON]['first-nt']:
            if self.model[START_CODON]['first-nt'][nt] == 0:
                self.model[START_CODON]['first-nt'][nt] = 1

        # generate possible list of codons the last nucleotide in the start codon can go to
        nt_list = ['a', 'c', 't', 'g']
        for i in range(0, len(nt_list)):
            for j in range(0, len(nt_list)):
                for k in range(0, len(nt_list)):
                    codon = nt_list[i] + nt_list[j] + nt_list[k]
                    self.model[START_CODON]['third-nt'][codon] = 1

        for i in range(0, self.query_count):
            ranges = self.metadata[i]['gene_coding_runs']

            for coding_range in ranges:
                start = coding_range['start']+3
                codon = self.seq[start: start + 3]
                self.model[START_CODON]['third-nt'][codon] += 1
                self.model[START_CODON]['third-nt']['transition_count'] += 1

    def train_internal_codons(self):
        self.build_internal_codon_tree()



    def build_noncoding_seq(self, query_index):
        print(query_index)
        mid_seqs = []
        end_seqs = []
        gene = self.metadata[query_index]

        # sequence from the start of the gene up to the start codon
        lower_bound = gene['query_start']
        upper_bound = gene['gene_coding_runs'][0]['start']
        mid_seqs.append(self.seq[lower_bound:upper_bound])

        print(lower_bound, upper_bound)

        num_runs = len(gene['gene_coding_runs'])

        for i in range(0, num_runs - 1):
            lower_bound = gene['gene_coding_runs'][i]['end'] + 1
            upper_bound = gene['gene_coding_runs'][i+1]['start']
            mid_seqs.append(self.seq[lower_bound:upper_bound])
            print(lower_bound, upper_bound)

        # sequence from the last coding nt to the end of the gene
        lower_bound = gene['gene_coding_runs'][-1]['end'] + 1
        upper_bound = gene['query_end']

        if lower_bound < upper_bound:
            end_seqs.append(self.seq[lower_bound:upper_bound])
            print('end', lower_bound, upper_bound, self.seq[upper_bound])

        # print(lower_bound, upper_bound)
        print('\n')

        return mid_seqs, end_seqs

    def build_internal_codon_tree(self):
        nt_list = ['a', 'c', 't', 'g']
        for i in range(0, len(nt_list)):
            for j in range(0, len(nt_list)):
                for k in range(0, len(nt_list)):
                    codon = nt_list[i] + nt_list[j] + nt_list[k]
                    self.model[INTERNAL_CODONS][codon] = {}

                    for u in range(0, len(nt_list)):
                        for m in range(0, len(nt_list)):
                            for n in range(0, len(nt_list)):
                                subcodon = nt_list[u] + nt_list[m] + nt_list[n]
                                self.model[INTERNAL_CODONS][codon][subcodon] = 1
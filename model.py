# {'a': {'a': 27675, 't': 26974, 'g': 22421, 'c': 22232, 'occurrences': 99302},
# 't': {'a': 18248, 't': 30094, 'g': 29028, 'c': 24297, 'occurrences': 101667},
# 'g': {'a': 23661, 't': 23942, 'g': 25855, 'c': 36272, 'occurrences': 109730},
# 'c': {'a': 29753, 't': 20608, 'g': 32441, 'c': 25049, 'occurrences': 107851}}

NONCODING = 'noncoding'
START_CODON = 'start-codon'
INTERNAL_CODONS = 'internal-codons'
OCC = 'occurrences'


class Model:
    model = {
        'noncoding': {
            'a': {'a': 0, 't': 0, 'g': 0, 'c': 0, 'occurrences': 0},
            't': {'a': 0, 't': 0, 'g': 0, 'c': 0, 'occurrences': 0},
            'g': {'a': 0, 't': 0, 'g': 0, 'c': 0, 'occurrences': 0},
            'c': {'a': 0, 't': 0, 'g': 0, 'c': 0, 'occurrences': 0}
        },
        'start-codon': {},
        'internal-codons': {
            'stop-codons': ['taa', 'tag', 'tga'],
            'a': {'t': 1.0},
            't': {'t': 0.5, 'g': 0.5},
            'g': {'t': 0.5, 'internal-codons': 0.5},
            'c': {'t': 1.0}
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
        self.gene_count = len(metadata)

    def train(self):
        self.train_noncoding_region()
        print(self.model[NONCODING])

    def train_noncoding_region(self):
        current_gene = 0
        current_gene_start = self.metadata[current_gene]['gene_start']

        for i in range(0, self.gene_count):
            seqs = self.build_noncoding_seq(i)
            for seq in seqs:
                for j in range(0, len(seq) - 1):
                    current_nt = seq[j]
                    next_nt = seq[j+1]
                    self.model[NONCODING][current_nt][next_nt] = self.model[NONCODING][current_nt][next_nt] + 1
                    self.model[NONCODING][current_nt][OCC] = self.model[NONCODING][current_nt][OCC] + 1

    def build_noncoding_seq(self, gene_index):
        print(gene_index)
        seqs = []
        gene = self.metadata[gene_index]

        # sequence from the start of the gene up to the start codon
        lower_bound = gene['gene_start']
        upper_bound = gene['gene_coding_runs'][0]['start']
        seqs.append(self.seq[lower_bound:upper_bound])

        print(lower_bound, upper_bound)

        num_runs = len(gene['gene_coding_runs'])

        for i in range(0, num_runs - 1):
            lower_bound = gene['gene_coding_runs'][i]['end'] + 1
            upper_bound = gene['gene_coding_runs'][i+1]['start']
            seqs.append(self.seq[lower_bound:upper_bound])
            print(lower_bound, upper_bound)

        # sequence from the last coding nt to the end of the gene
        lower_bound = gene['gene_coding_runs'][-1]['end'] + 1
        upper_bound = gene['gene_end']
        seqs.append(self.seq[lower_bound:upper_bound])

        print(lower_bound, upper_bound)
        print('\n')

        return seqs

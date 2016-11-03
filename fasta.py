class Reader():
    """Simple reader of dna sequences in FASTA format.

    """

    def __init__(self, filename, verbose=False):
        """

        :param filename: location of file with FASTA-formated sequences
        """

        self.ordered_ids = []
        self.sequences = dict()
        if verbose:
            print('Reading file: ', end='')

        with open(filename, 'r') as fin:
            current_seq = None
            current_seq_id = None
            for line in fin:
                if line[0] == ">":
                    # line with FASTA label
                    # save previous sequence
                    if current_seq_id != None:
                        self._add_seq(current_seq_id, current_seq)
                    current_seq = ""
                    current_seq_id = line[1:-1]
                elif len(line) > 1:
                    # dna sequence
                    current_seq += line[:-1]
            # save last sequence
            self._add_seq(current_seq_id, current_seq)
            if verbose:
                print(len(self.sequences), "sequences loaded.")

    def _add_seq(self, seq_id, seq):
        if seq_id in self.sequences:
            print("ERROR: Duplicate sequence ids.")
            return

        self.sequences[seq_id] = seq
        self.ordered_ids.append(seq_id)

    def __len__(self):
        return len(self.sequences)

    def __getitem__(self, key):
        return self.sequences[key]

    def __contains__(self, item):
        return item in self.sequences

    def __iter__(self):
        return self.sequences.__iter__()


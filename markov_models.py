

class HMM():
    def __init__(self, transition_matrix, emission_matrix):
        self.transition_matrix = transition_matrix
        self.emission_matrix = emission_matrix

    def path_prob(self, path):
        prob = 1 / len(self.transition_matrix)
        prev = path[0]
        for curr in path[1:]:
            prob *= self.transition_matrix[prev][curr]
            prev = curr

        return prob

    def output_prob_given_path(path, output):
        prob = 1
        for o, s in zip(output, path):
            prob *= self.emission_matrix[s][o]

        return prob

    def most_probably_path_given_output(output):
        pass


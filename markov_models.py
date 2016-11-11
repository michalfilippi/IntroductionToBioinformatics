

class HMM:
    """Class for representation of hidden markov models.
    """

    def __init__(self, transition_matrix, emission_matrix):
        """
        :param transition_matrix: dictionary of dictionaries representing
        probability matrix for transitions
        :param emission_matrix: dictionary of dictionaries representing
        probability matrix for emissions
        """

        self.transition_matrix = transition_matrix
        self.emission_matrix = emission_matrix

    def path_prob(self, path):
        """Calculates probability of given path using probability matrix for
        transitions.

        :param path: path over states
        :return: probability of given path
        """

        prob = 1 / len(self.transition_matrix)
        prev = path[0]
        for curr in path[1:]:
            prob *= self.transition_matrix[prev][curr]
            prev = curr

        return prob

    def output_prob_given_path(self, path, output):
        """Calculates probability of given output assuming path over states as
        given in 'path' using probability matrix for emissions.

        :param path: path over states
        :param output: output of HMM
        :return: probability of given output
        """

        prob = 1
        for o, s in zip(output, path):
            prob *= self.emission_matrix[s][o]

        return prob

    def most_probable_path_given_output(self, output):
        """Calculates the most probable path over states assuming seen output as
        given in 'output'.
        Algorithm is called Viterbi Algorithm.

        :param output: output of HMM
        :return: tuple (path, prob), where path is the most probable path over
        states and prob is it's probability
        """
        probs = dict()
        prev = dict()
        states = list(self.transition_matrix.keys())

        for state in states:
            emiss_prob = self.emission_matrix[state][output[0]]
            probs[(state, 0)] = 1 / len(self.transition_matrix) * emiss_prob

        for i, o in enumerate(output):
            if i == 0:
                continue

            for state in states:
                emiss_prob = self.emission_matrix[state][o]
                predecessors = []
                for prev_state in states:
                    trans_prob = self.transition_matrix[prev_state][state]
                    prob = probs[(prev_state, i - 1)] * trans_prob * emiss_prob
                    predecessors.append((prob, (prev_state, i - 1)))
                best_predecessor = max(predecessors)

                probs[(state, i)], prev[(state, i)] = best_predecessor

        # traceback most probable path
        reversed_path = []
        end_nodes = [(probs[(state, len(output) - 1)], (state, len(output) - 1))
                     for state in states]

        prob, current = max(end_nodes)
        reversed_path.append(current[0])

        while current[1] != 0 :
            current = prev[current]
            reversed_path.append(current[0])

        return ''.join(reversed(reversed_path)), prob
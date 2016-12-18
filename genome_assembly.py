

class DeBrujinGraph:
    def __init__(self, k_mers):
        """Creates De Brujin graph from a list of strings.
        https://en.wikipedia.org/wiki/De_Bruijn_graph

        :param k_mers: list of strings
        """

        self.vertices = set()
        self.edges = dict()

        for k_mer in k_mers:
            s_1 = k_mer[:-1]
            s_2 = k_mer[1:]

            self.vertices.add(s_1)
            self.vertices.add(s_2)

            if s_1 not in self.edges:
                self.edges[s_1] = []

            if s_2 not in self.edges[s_1]:
                self.edges[s_1].append(s_2)

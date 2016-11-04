def n_dimensional_iterator(*lengths):
    """Generator of all positions in len(lengths)-dimensional cube,
    where n-dimension has length lengths[n].

    :param lengths: lengths of dimensions
    :return: generator
    """
    if len(lengths) == 0:
        yield []
        return
    for i in range(lengths[0]):
        for pos in n_dimensional_iterator(*lengths[1:]):
            yield [i] + pos



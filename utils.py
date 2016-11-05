def n_dim_cube_iterator(*lengths):
    """Generator of all positions in len('lengths')-dimensional cube,
    where n-dimension has length 'lengths'[n].

    :param lengths: lengths of dimensions
    :return: generator
    """

    if len(lengths) == 0:
        yield tuple()
        return

    for i in range(lengths[0]):
        for pos in n_dim_cube_iterator(*lengths[1:]):
            yield (i,) + pos


def n_dim_cube_predecessors(position):
    """Generates all preceding positions to 'position'. Preceding position is
    position with all coordinates equal or lesser by one then original position.
    First returned position is identical position.

    :param position: position in n-dim cube
    :return: generator of all preceding positions
    """

    if len(position) == 0:
        yield tuple()
        return

    for pos in n_dim_cube_predecessors(position[1:]):
        yield (position[0],) + pos
        if position[0] - 1 >= 0:
            yield (position[0] - 1,) + pos


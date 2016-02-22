from . import sg_align


def align_pair(reference, reference_rc, pair):
    """
    Align a pair of markers to the forward reference sequence. The reverse
    complement is used to align the second element of the pair (which is also
    reverse complemented).

    :arg reference: Reference sequence to align to.
    :type reference: str
    :arg reference_rc: Reverse complement of the reference sequence.
    :type reference_rc: str
    :arg pair: A pair (forward, reverse) of markers to align.
    :type pair: list[str, str]

    :returns: A tuple (score, position) of the best alignment.
    :rtype: tuple(int, int)
    """
    right_alignment = sg_align.align(reference_rc, pair[1])
    left_alignment = sg_align.align(reference, pair[0])

    return ((left_alignment.distance, left_alignment.position),
        (right_alignment.distance, len(reference) - right_alignment.position))

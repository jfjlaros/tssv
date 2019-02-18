from . import sg_align_i


def align_pair(reference, reference_rc, pair, indel_score=1):
    """Align a pair of markers to the forward reference sequence. The reverse
    complement is used to align the second element of the pair (which is also
    reverse complemented).

    :arg str reference: Reference sequence to align to.
    :arg str reference_rc: Reverse complement of the reference sequence.
    :arg list pair: A pair (forward, reverse) of markers to align.
    :arg int indel_score: Penalty score for insertions and deletions per
        nucleotide.

    :returns tuple: A tuple (score, position) of the best alignment.
    """
    right_alignment = sg_align_i.align_i(reference_rc, pair[1], indel_score)
    left_alignment = sg_align_i.align_i(reference, pair[0], indel_score)

    #raise ValueError(str(left_alignment))
    return (
        (left_alignment['distance'], left_alignment['position']),
        (right_alignment['distance'], len(reference) - right_alignment['position']))

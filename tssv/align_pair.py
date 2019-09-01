from .sg_align import align, align_sse


def align_pair(reference, reference_rc, pair, indel_score=1, method_sse=False):
    """Align a pair of markers to the forward reference sequence. The reverse
    complement is used to align the second element of the pair (which is also
    reverse complemented).

    :arg str reference: Reference sequence to align to.
    :arg str reference_rc: Reverse complement of the reference sequence.
    :arg list pair: A pair (forward, reverse) of markers to align.
    :arg int indel_score: Penalty score for insertions and deletions per
        nucleotide.
    :arg bool method_sse: Use SSE2 alignment implementation.

    :returns tuple: A tuple (score, position) of the best alignment.
    """
    if method_sse:
        left = align_sse(reference, pair[0], indel_score)
        right = align_sse(reference_rc, pair[1], indel_score)
    else:
        left = align(reference, pair[0], indel_score)
        right = align(reference_rc, pair[1], indel_score)

    return (
        (left['distance'], left['position']),
        (right['distance'], len(reference) - right['position']))

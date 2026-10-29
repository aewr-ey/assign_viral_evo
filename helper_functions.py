import math

def jukes_cantor(reference_sequence: str, distant_sequence: str) -> float:
    """The Jukes-Cantor correction for estimating genetic distances
    calculated with Hamming distance.
    Should return genetic distance with the same unit as if not corrected.

    Parameters
    ----------
    reference_sequence: str
        A string of nucleotides in a sequence used as a reference
        in an alignment with other (e.g. AGGT-GA)
    distant_sequence: str
        A string of nucleotides in a sequence after the alignment
        with a reference (e.g. AGC-AGA)

    Returns
    -------
    float
        The Jukes-Cantor corrected genetic distance using Hamming distance.
        For example 1.163.

    """
    mismatch = 0 
    total = 0
    hamming_proportion = 0

    ## hamming proportion calc - p 
    for i in range(len(reference_sequence)):
        ref = reference_sequence[i]
        dist = distant_sequence[i]

        if ref == '-' or dist == '-':
            continue

        total = total + 1 

        if ref != dist:
            mismatch = mismatch + 1

    hamming_proportion = mismatch/total 

    # jukes-cantor correction - d
    distance = -3/4 * math.log(1- (4/3*(hamming_proportion)))

    # apply distance for all the nucleotides 
    return distance * total





def kimura_two_parameter(reference_sequence: str, distant_sequence: str) -> float:
    """The Kimura Two Parameter correction for estimating genetic distances
    calculated with Hamming distance.
    Should return genetic distance with the same unit as if not corrected.

    Parameters
    ----------
    reference_sequence: str
        A string of nucleotides in a sequence used as a reference
        in an alignment with other (e.g. AGGT-GA)
    distant_sequence: str
        A string of nucleotides in a sequence after the alignment
        with a reference (e.g. AGC-AGA)

    Returns
    -------
    float
        The Kimura corrected genetic distance using Hamming distance.
        For example 1.196.

    """
    raise NotImplementedError()

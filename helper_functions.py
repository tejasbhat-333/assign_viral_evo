import math

PURINES = {"A", "G"}
PYRIMIDINES = {"C", "T"}
NUCLEOTIDES = PURINES | PYRIMIDINES


def _compared_pairs(reference_sequence: str, distant_sequence: str) -> list:
    """Pairs of nucleotides at positions without indels or ambiguous characters."""
    pairs = []
    for ref_nt, dist_nt in zip(reference_sequence.upper(), distant_sequence.upper()):
        if ref_nt in NUCLEOTIDES and dist_nt in NUCLEOTIDES:
            pairs.append((ref_nt, dist_nt))
    return pairs


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
    pairs = _compared_pairs(reference_sequence, distant_sequence)
    length = len(pairs)
    if length == 0:
        return 0.0

    mismatches = sum(1 for ref_nt, dist_nt in pairs if ref_nt != dist_nt)
    p = mismatches / length                      # proportion of differing sites
    if p >= 0.75:
        return float("inf")                      

    distance_per_site = -0.75 * math.log(1 - 4 * p / 3)
    return abs(distance_per_site * length)       # back to the unit of the Hamming distance 


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
    pairs = _compared_pairs(reference_sequence, distant_sequence)
    length = len(pairs)
    if length == 0:
        return 0.0

    transitions = sum(1 for ref_nt, dist_nt in pairs
                      if ref_nt != dist_nt
                      and ((ref_nt in PURINES) == (dist_nt in PURINES)))
    transversions = sum(1 for ref_nt, dist_nt in pairs
                        if (ref_nt in PURINES) != (dist_nt in PURINES))
    p = transitions / length                     # proportion of transitions
    q = transversions / length                   # proportion of transversions

    if 1 - 2 * p - q <= 0 or 1 - 2 * q <= 0:
        return float("inf")                      

    distance_per_site = -0.5 * math.log(1 - 2 * p - q) - 0.25 * math.log(1 - 2 * q)
    return abs(distance_per_site * length)
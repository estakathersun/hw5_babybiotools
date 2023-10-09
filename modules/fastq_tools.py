def gc_filter(seq: str, gc_bounds: tuple | int):
    """
Filters a sequence based on GC content bounds.
Args:
    - seq (str): The nucleotide sequence.
    - gc_bounds (tuple or int): A tuple containing the minimum and maximum
    acceptable percent of GC nucleotides in the sequence.
    A tuple of two integers representing the lower and upper bounds for GC
    content.
    If an integer is provided, it is treated as the maximum bound and the
    minimum bound is set to 0.
Returns:
    - bool: True if the GC content of the sequence is within the acceptable
    bounds, False otherwise.
"""
    gc_count = seq.upper().count('G') + seq.upper().count('C')
    gc_percent = gc_count / len(seq) * 100
    if isinstance(gc_bounds, int):
        gc_bounds = (0, gc_bounds)
    return gc_bounds[0] <= gc_percent <= gc_bounds[1]


def length_filter(seq: str, length_bounds: tuple | int):
    """
Filters a sequence based on length bounds.
Args:
    - seq (str): The nucleotide sequence.
    - length_bounds (tuple or int): A tuple containing the minimum and maximum
    acceptable length of the sequence.
    If an integer is provided, it is treated as the maximum length and the
    minimum length is set to 0.
Returns:
    - bool: True if the length of the sequence is within the acceptable bounds,
      False otherwise.
"""
    if isinstance(length_bounds, int):
        length_bounds = (0, length_bounds)
    return length_bounds[0] <= len(seq) <= length_bounds[1]


def quality_filter(seq: str, qual_seq: str, quality_threshold: int):
    """
Filters a FASTQ sequence based on quality threshold.
Args:
    - seq (str): The nucleotide sequence.
    - qual_seq (str): The string of ASCII symbols where each symbol represents
    quality score for each base in the sequence.
    - quality_threshold (int): The minimum acceptable quality score for the
    sequence.
Returns:
    - bool: True if the average quality score of the sequence is greater than
    or equal to the quality threshold, False otherwise.
"""
    q_total = 0
    for symbol in qual_seq:
        q_total += ord(symbol) - 33
    q_mean = q_total / len(seq)
    return quality_threshold <= q_mean


def run_fastq_filter(seqs: dict,
                     gc_bounds: tuple | int = (0, 100),
                     length_bounds: tuple | int = (0, 2**32),
                     quality_threshold: int = 0):
    """ Filters a dictionary of FASTQ sequences based on GC content, sequence
    length, and quality threshold.
Args:
    - seqs (dict): A dictionary of FASTQ sequences where the keys are the sequence
    IDs and the values are tuples of the form (sequence, quality scores).
    - gc_bounds (tuple | int): A tuple of two integers representing the lower
    and upper bounds for GC content.
    Alternatively, a single integer can be provided as a shorthand for setting
    the upper bound to that value and the lower bound to 0.
    - length_bounds (tuple | int): A tuple of two integers representing
    the lower and upper bounds for sequence length.
    Alternatively, a single integer can be provided as a shorthand for setting
    the upper bound to that value and the lower bound to 0.
    - quality_threshold (int): The minimum acceptable quality score for each
    base in the sequence.
Returns:
    - dict: A filtered dictionary of FASTQ sequences where the keys are the
    sequence IDs and the values are tuples of the form (sequence, quality scores).

Example:
    seqs = {
        "@seq1": ("ATCGATCG", "!##$%&'"),
        "@seq2": ("AGCTAGCTAGCT", "!!!!@@@@####")
    }
    filtered_seqs = run_fastq_filter(seqs, gc_bounds = (30, 70), length_bounds = 10, quality_threshold = 20)
"""
    filtered_seqs = {}
    for idx, value in seqs.items():
        seq = value[0]
        qual_seq = value[1]
        if (
            gc_filter(seq, gc_bounds) is True &
            length_filter(seq, length_bounds) is True &
            quality_filter(seq, qual_seq, quality_threshold) is True
        ):
            filtered_seqs[idx] = value
    return filtered_seqs

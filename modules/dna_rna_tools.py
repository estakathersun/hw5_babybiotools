COMPLEMENTS_DNA = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
                   'a': 't', 't': 'a', 'g': 'c', 'c': 'g'}
COMPLEMENTS_RNA = {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G',
                   'a': 'u', 'u': 'a', 'g': 'c', 'c': 'g'}


def is_nucleic_acid(seq: str) -> bool:
    """
    Checks if the provided nucleotide sequences contain only acceptable nucleic acid types.
Args:
    - seq (str): string containing nucleotide sequence to be checked.
Returns:
    - bool: True if a sequence is acceptable nucleic acid types, False otherwise.
Raises:
    - ValueError: If a sequence contains both 'T' and 'U' nucleic acid types.
    """
    seq_upper = seq.upper()
    t_in_seq, u_in_seq = 'T' in seq_upper, 'U' in seq_upper
    if t_in_seq and u_in_seq:
        raise ValueError('Unacceptable nucleic acid type')
    return True


def transcribe(seq: str) -> str:
    """
Transcribes DNA sequence into RNA sequence.
Args:
    - seq (str): string containing the DNA sequence to be transcribed into RNA.
Returns:
    - str: returns a string representing the transcribed sequence.
    """
    return seq.replace('T', 'U').replace('t', 'u')


def reverse(seq: str) -> str:
    """
Reverses DNA or RNA sequence.
Args:
    - seq (str): string containing the DNA or RNA sequence to be reversed.
Returns:
   - str: returns a string representing the reversed sequence.
   """
    return seq[::-1]


def complement(seq: str) -> str:
    """
Complements DNA or RNA sequence.
Args:
       - seq (str): string containing the DNA or RNA sequence to be complemented.
Returns:
   - str: returns a string representing the complemented sequence.
    """
    complemented_seq = ''
    if seq.upper().count('U') > 0:
        for s in seq:
            complemented_seq += COMPLEMENTS_RNA.get(s)
    else:
        for s in seq:
            complemented_seq += COMPLEMENTS_DNA.get(s)
    return complemented_seq


def reverse_complement(seq: str) -> str:
    """
Complements and reverses DNA or RNA sequence.
Args:
   - seq (str): string containing the DNA or RNA sequence to be complemented and reversed.
Returns:
   - str: returns a string representing the complemented and reversed sequence.
    """
    return reverse(complement(seq))


DNA_RNA_OPERATIONS = {'reverse': reverse, 'transcribe': transcribe,
                      'complement': complement, 'reverse_complement': reverse_complement}

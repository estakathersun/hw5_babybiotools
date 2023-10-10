import modules.aminoacid_tools as aa
import modules.fastq_tools as fq
import modules.dna_rna_tools as na
from typing import Union


def run_dna_rna_tools(*args: str) -> Union[list, str]:
    """
Runs a specified DNA or RNA operation on one or more DNA or RNA sequences.
Args:
    - *args: Variable-length argument list containing the DNA or RNA sequences
    to be operated on, followed by the operation to be performed.
        - The arguments can be either strings representing DNA sequences or RNA
        sequences.
        - The final argument must be a string representing the operation to be
        performed, and must be one of the keys in the DNA_RNA_OPERATIONS dictionary.
Returns:
    - str or list: If only one sequence is provided, returns a single string
    representing the output sequence after the specified operation is performed.
     If multiple sequences are provided, returns a list of strings representing
     the output sequences for each input sequence after the specified operation
     is performed.
Example:
    seq1 = "ATCG"
    seq2 = "GCUA"
    run_dna_rna_tools(seq1, "complement") -> "TAGC"
    run_dna_rna_tools(seq1, seq2, "reverse_complement") -> ["CGAT", "UAGC"]
"""
    seqs, operation = args[0:-1], args[-1]
    output_seqs = []
    for seq in seqs:
        if na.is_nucleic_acid(seq) is True:
            output_seqs.append(na.DNA_RNA_OPERATIONS[operation](*seqs))
    if len(seqs) == 1:
        return output_seqs[0]
    else:
        return output_seqs


def run_aminoacid_tools(*seqs: str, operation: str) -> str:
    """Runs a specified operation on one or more amino acid sequences.
 Args:
    -  *seqs: Variable-length argument list containing the amino acid sequences to be operated on.
         - The arguments must be strings representing amino acid sequences.
     - operation: A string representing the operation to be performed.
         - The operation must be one of the keys in the OPERATION_DICT dictionary.
 Returns:
     str: returns a string representing the input sequences and result of specified operation performance.
 Raises:
    - ValueError: If the operation value is not specified or is not one of the supported operations.
    - TypeError: If any of the input sequences are not valid amino acid sequences.

 Example:
     seq1 = "MSRSLLLRFLLFLLLLPPLP"
     seq2 = "MSRSLLLRFLLFLLLLPPLP"
     run_aminoacid_tools(seq1, seq2, operation="calculate_molecular_weight") ->
     """
    if operation == '':
        raise ValueError('Operation value is not specified')
    if operation not in aa.OPERATION_DICT:
        raise ValueError(f'Incorrect operation value\nSupported operations: {list(aa.OPERATION_DICT.keys())}')
    output = ''
    for seq in seqs:
        aa.is_peptide(seq)
        output += aa.OPERATION_DICT[operation](seq)
        output += '\n'
    return output


def run_fastq_filter(seqs: dict,
                     gc_bounds: Union[tuple, int] = (0, 100),
                     length_bounds: Union[tuple, int] = (0, 2 ** 32),
                     quality_threshold: int = 0) -> dict:
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
    filtered_seqs = run_fastq_filter(seqs, gc_bounds = (30, 70),
    length_bounds = 10, quality_threshold = 20)
"""
    filtered_seqs = {}
    for seq_id, value in seqs.items():
        seq = value[0]
        qual_seq = value[1]
        if (
                fq.gc_filter(seq, gc_bounds) is True &
                fq.length_filter(seq, length_bounds) is True &
                fq.quality_filter(seq, qual_seq, quality_threshold) is True
        ):
            filtered_seqs[seq_id] = value
    return filtered_seqs

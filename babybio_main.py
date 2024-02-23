import os
from typing import Union, Dict


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
        if is_nucleic_acid(seq) is True:
            output_seqs.append(DNA_RNA_OPERATIONS[operation](*seqs))
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
    if operation not in OPERATION_DICT:
        raise ValueError(f'Incorrect operation value\nSupported operations: {list(OPERATION_DICT.keys())}')
    output = ''
    for seq in seqs:
        is_peptide(seq)
        output += OPERATION_DICT[operation](seq)
        output += '\n'
    return output


def run_fastq_filter(input_path: str,
                     output_filename: str = None,
                     gc_bounds: Union[tuple, int] = (0, 100),
                     length_bounds: Union[tuple, int] = (0, 2 ** 32),
                     quality_threshold: int = 0):
    """
    Reads a FASTQ file, filters sequences based on GC content, sequence
    length and quality threshold, and writes the filtered sequences to
    a new file.
    Args:
        - input_path (str): A string representing the path to the input FASTQ
        file.
        - output_filename (str): A string representing the name of the output file.
        If None, the original filename will be used as output filename.
    as the input filename with the extension removed
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
    """
    seqs = read_fastq_file(input_path)
    filtered_seqs = {}
    for seq_id, value in seqs.items():
        seq = value[0]
        qual_seq = value[1]
        if (
                gc_filter(seq, gc_bounds) is True
                and length_filter(seq, length_bounds) is True
                and quality_filter(seq, qual_seq, quality_threshold) is True
        ):
            filtered_seqs[seq_id] = value
    if output_filename is None:
        # remove extension
        output_filename = os.path.basename(input_path)[:-6]
    write_output_file(filtered_seqs, output_filename)


#### fastq_tools ####

def gc_filter(seq: str, gc_bounds: tuple | int) -> bool:
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


def length_filter(seq: str, length_bounds: tuple | int) -> bool:
    """
    Filters a sequence based on length bounds.
    Args:
        - seq (str): The nucleotide sequence.
        - length_bounds (tuple or int): A tuple containing the minimum and
        maximum acceptable length of the sequence.
        If an integer is provided, it is treated as the maximum length and the
        minimum length is set to 0.
    Returns:
        - bool: True if the length of the sequence is within the acceptable
        bounds, False otherwise.
    """
    if isinstance(length_bounds, int):
        length_bounds = (0, length_bounds)
    return length_bounds[0] <= len(seq) <= length_bounds[1]


def quality_filter(seq: str, qual_seq: str, quality_threshold: int) -> bool:
    """
    Filters a FASTQ sequence based on quality threshold.
    Args:
        - seq (str): The nucleotide sequence.
        - qual_seq (str): The string of ASCII symbols where each symbol
        represents quality score for each base in the sequence.
        - quality_threshold (int): The minimum acceptable quality score for the
        sequence.
    Returns:
        - bool: True if the average quality score of the sequence is greater
        than or equal to the quality threshold, False otherwise.
    """
    q_total = 0
    for symbol in qual_seq:
        q_total += ord(symbol) - 33
    q_mean = q_total / len(seq)
    return quality_threshold <= q_mean


def read_fastq_file(input_path: str) -> dict:
    """read FASTQ file to dictionary {'name': ('seq', 'quality')}"""
    if not os.path.isfile(input_path):
        raise ValueError(f'{input_path} is not a file!')
    with open(input_path) as input_file:
        line = input_file.readline()
        if line == '':
            raise ValueError(f'{input_path} is empty!')
        # на случай наличия посторонних строк до смысловых записей:
        while not line.startswith('@'):
            line = input_file.readline()
        reads_list = []
        fastq_dict = {}
        while line != '':
            reads_list.append(line.strip())
            if len(reads_list) == 4:
                seqname = reads_list[0].split()[0]
                seq = reads_list[1]
                quality = reads_list[3]
                fastq_dict[seqname] = seq, quality
                reads_list = []
            line = input_file.readline()
    return fastq_dict


def write_output_file(filtered_dict: dict, output_filename: str):
    """write FASTQ dictionary {'name': ('seq', 'quality')} to FASTQ file"""
    output_dir = 'fastq_filtrator_results'
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    elif not os.path.isdir(output_dir):
        raise ValueError(f'{output_dir} is not a directory!')
    ext = '.fastq'
    # если output_filename будет передан с расширением:
    if output_filename[-6:] == ext:
        ext = ''
    output_path = os.path.join(output_dir, f'{output_filename}{ext}')
    if os.path.isfile(output_path):
        raise FileExistsError(f'{output_path} already exists')
    with open(output_path, mode='w') as out_file:
        div = ''
        for key, value in filtered_dict.items():
            out_file.write(div + key + '\n' + value[0] + '\n+\n' + value[1])
            div = '\n'


#### dna_rna_tools ####

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

#### aminoacid_tools ####

ALL_AMINOACIDS = {
    'A', 'R', 'N', 'D', 'C', 'H', 'G', 'Q', 'E', 'I',
    'L', 'K', 'M', 'P', 'S', 'Y', 'T', 'W', 'F', 'V'
}

AA_WEIGHT = {
    'G': 57.051, 'A': 71.078, 'S': 87.077, 'P': 97.115, 'V': 99.131,
    'T': 101.104, 'C': 103.143, 'I': 113.158, 'L': 113.158, 'N': 114.103,
    'D': 115.087, 'Q': 128.129, 'K': 128.172, 'E': 129.114, 'M': 131.196,
    'H': 137.139, 'F': 147.174, 'R': 156.186, 'Y': 163.173, 'W': 186.210
}


def is_peptide(seq: str) -> bool:
    """Check whether the incoming sequence is an aminoacid"""
    return set(seq).issubset(ALL_AMINOACIDS)  # if set(seq) <= all_aminoacids


def find_cleavage_sites(seq: str, motif: list) -> list:
    """Find cleavage sites for motif-specific proteases.
            Arguments:
            - seq - string sequence to be analyzed
            - motif - subsequence to be found in a sequence. Subsequence is specified as
            list of lists. Each nested list means more than one possible aminoacid at a
            single position (checked by OR condition).
            Returns:
            - list of cleavage sites coordinates (C-end aminoacid of *potentially*
            cleaved sequence)
            """
    cleavage_sites = []
    seq_idx = 0
    while seq_idx < len(seq):
        motif_idx = 0
        chars_at_motif_idx = motif[motif_idx]
        seq_char = seq[seq_idx]
        if seq_char in chars_at_motif_idx:
            motif_idx += 1
            while motif_idx < len(motif):
                chars_at_motif_idx = motif[motif_idx]
                seq_char = seq[seq_idx + motif_idx]
                if seq_char in chars_at_motif_idx:
                    motif_idx += 1
                else:
                    break
            if motif_idx == len(motif):
                cleavage_sites.append(seq_idx + motif_idx)
        seq_idx += 1
    return cleavage_sites


def get_cleavage_sites(seq: str) -> str:
    """
        Return amount and coordinates of cleavage sites for proteases,
        specified in motif_dict
        Args:
            - seq (str): aminoacid sequence to be analyzed
        Returns:
            - output (str): string containing info about amount and position of cleavage sites for each protease
            """
    # this dict can be expanded with other motifs according to the item format
    motif_dict = {
        'Caspase 3': [['D'], ['M'], ['Q'], ['D']],
        'Caspase 6': [['V'], ['E'], ['H', 'I'], ['D']],
        'Caspase 7': [['D'], ['E'], ['V'], ['D']],
        'Enterokinase': [['D', 'E'], ['D', 'E'], ['D', 'E'], ['K']]
    }
    output = f'{seq}\n'
    for motif_name, motif_value in motif_dict.items():
        sites = find_cleavage_sites(seq, motif_value)
        output += f'{len(sites)} protease cleavage site(s) for {motif_name}: {sites}\n'
    return output


def calculate_aa_percentage(seq: str) -> str:
    """
    Calculates the percentage of amino acids in the entered amino acid sequence
    Arguments:
        - seq (str): amino acid sequences to be analyzed
    Return:
        - str: a string with the percentage of each amino acid
    """
    amino_acid_counts: Dict[str, int] = {}  # dict to store count of each amino acid
    for amino_acid in seq:
        if amino_acid in amino_acid_counts:
            amino_acid_counts[amino_acid] += 1
        else:
            amino_acid_counts[amino_acid] = 1
    total_amino_acids = len(seq)
    amino_acid_percentages = {}  # dict to store each amino acid and its %
    for amino_acid, count in amino_acid_counts.items():
        percentage = round(((count / total_amino_acids) * 100), 2)
        amino_acid_percentages[amino_acid] = percentage
    return f'Amino acids percentage of the sequence {seq}: {amino_acid_percentages}'


def calculate_molecular_weight(seq: str) -> str:
    """
    Calculates the molecular weight of entered amino acid sequence
    Arguments:
        - seq (str): amino acid sequences to be analyzed
    Return:
        - str: a string with the molecular weight value for amino acid sequence
    """
    weight = 18.02  # for the H and OH at the termini
    for amino_acid in seq:
        weight += AA_WEIGHT[amino_acid]
    return f'Molecular weight of the sequence {seq}: {round(weight, 2)} Da'


OPERATION_DICT = {
    'get_cleavage_sites': get_cleavage_sites,
    'calculate_molecular_weight': calculate_molecular_weight,
    'calculate_aa_percentage': calculate_aa_percentage,
}

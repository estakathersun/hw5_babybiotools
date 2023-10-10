from typing import Dict

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

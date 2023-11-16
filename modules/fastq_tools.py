import os


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

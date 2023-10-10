# BabyBioTools

## Overview

A library containing modules for performing simple operations with aminoacid, nucleotide and fastq sequences.
Contains the next modules:

* [DNA and RNA processing](#dnarnatools)
* [FASTQ sequences filtering](#fastqtools)
* [aminoacid sequences analysis](#aminoacidtools)

## System requirements

- [Python](https://www.python.org/downloads/) version 3.9 or higher

## Modules

To use BabyBioTools modules import one or more modules described below into your `*.py` script, for example:

```python
from babybio_main import run_fastq_filter
```

> Modules scripts can be found in `modules/` dir.

### dna_rna_tools

```python
from babybio_main import run_dna_rna_tools
```

Runs a specified DNA or RNA operation on one or more DNA or RNA sequences.

Supported operations:

* `transcribe` - converts incoming DNA sequences to RNA
* `reverse` - returns reversed sequences
* `complement` - returns DNA or RNA sequences complementary to the input ones
* `reverse_complement` - returns reversed and complemented sequences

Arguments:

* `*args` - variable-length argument list containing the DNA or RNA sequences to be operated on, followed by the
  operation to be performed.
  The arguments can be either strings representing DNA sequences or RNA sequences.
* *The final argument* must be a string representing the operation to be performed, and must be one of the supported
  operations.

If only one sequence is provided, returns a single string representing the output sequence after the specified operation
is performed. If multiple sequences are provided, returns a list of strings representing the output sequences for each
input sequence after the specified operation is performed.

Example:

```python
seq1 = "ATCG"
seq2 = "GCUA"
run_dna_rna_tools(seq1, "complement") -> "TAGC"
run_dna_rna_tools(seq1, seq2, "reverse_complement") -> ["CGAT", "UAGC"]
```

### fastq_tools

```python
from babybio_main import run_fastq_filter
```

Filters a dictionary of FASTQ sequences based on GC content, sequence length, and quality threshold.
Takes as input a dictionary of FASTQ sequences and options for filtration.

Arguments:

* `seq` - a dictionary of FASTQ sequences where the keys are the sequence IDs and the values are tuples of the
  form (`sequence`, `quality scores`).
* `gc_bounds` - a tuple of two integers representing the lower and upper bounds for GC content. Alternatively, a single
  integer can be provided as a shorthand for setting the upper bound to that value and the lower bound to 0. `(0, 100)`
  as default.
* `length_bounds` - a tuple of two integers representing the lower and upper bounds for sequence length. Alternatively,
  a single integer can be provided as a shorthand for setting the upper bound to that value and the lower bound to
  0. `(0, 2**32)` as default.
* `quality_threshold` - the minimum acceptable quality score for each base in the sequence. `0` as default.

Returns an equivalent dictionary of filtered FASTQ sequences where the keys are the sequence IDs and the values are
tuples of the form (`sequence`, `quality scores`).

Example:

```python
# input 
seqs = {"seq1": ("ATCGATCG", "!##$%&'"), 
        "seq2": ("AGCTAGCTAGCT", "!!!!@@@@####")}
filtered_seqs = run_fastq_filter(seqs, gc_bounds = (10, 70), length_bounds = 10, quality_threshold = 1)
print(filtered_seqs)

# output
{'seq1': ('ATCGATCG', "!##$%&'")}
```

### aminoacid_tools

```python
from babybio_main import run_aminoacid_tools
```

Runs a specified operation on one or more amino acid sequences.

Arguments:

* `*seqs` - variable-length argument list containing the aminoacid sequences to be operated on. The arguments must be
  strings representing amino acid sequences.
* `operation` - a string representing the operation to be performed.
  - The operation must be one of the supported operations.

Returns a string representing the input sequences and result of specified operation performance.

Raises:

* ValueError: If the operation value is not specified or is not one of the supported operations.
* TypeError: If any of the input sequences are not valid amino acid sequences.

**Supported operations:**

* `calculate_aa_percentage` - calculate the percentage of amino acids in a sequence

```python
run_aminoacid_tools('ARG', operation='calculate_percentage')  # input

Amino acids percentage of the sequence ARG: {A: 33.33, R: 33.33, G: 33.33}  # output
```

* `calculate_molecular_weight` - calculate the molecular weight of the input amino acid sequences based on the mass of
  each amino acid residue.

Reference values for the masses of amino acid residues are taken
from [The University of Washington's Proteomics Resource (UWPR)](https://proteomicsresource.washington.edu/protocols06/masses.php)
and rounded to three decimal places. The calculations took into account the presence of *H* and *OH* groups at the
termini of the sequences.
The input is a string with amino acid sequence. The output is a string with the molecular weight of sequence in
Daltons (the result is rounded to two decimal places):

```python
run_aminoacid_tools('ARG', operation='calculate_molecular_weight')  # input
Molecular weight of the sequence ARG: 302.33 Da  # output
```

* `get_cleavage_sites` - return amount and coordinates of cleavage sites for motif-specific proteases (casp3, casp6,
  casp7, enterokinase).

The function finds traces of motifs in amino acid sequence that can be recognized by next site-specific proteases:
caspases 3, 6, 7 and enterokinase, then calculates amount of each protease's site and its coordinates. The coordinate is
a position of amino acid in C-end of potentially cleaved peptide. Motifs that can be recognized by these proteases were
taken from [PeptideCutter](https://web.expasy.org/peptide_cutter/peptidecutter_enzymes.html) documentation.
> Some of the proteases permit more than one possible amino acids at a single position. For example, Caspase 6 protease
> motif is defined as `V, E, H or I, D`. It means there are two possible sequences `...VEHD...` and `...VEID...`
> containing the protease motif. Please note the function implementation supports 'OR' condition only.
> Optionally you can add motifs for site-specific proteases of your interest by editing `aminoacid_tools.py` script (by
> adding them in `motif_dict` dictionary)

Example:

```python
# input
run_aminoacid_tools('ESDMQDMDSGISLDNDDEEKMQ', operation='get_cleavage_sites')

# output
ESDMQDMDSGISLDNDDEEKMQ
1 protease cleavage site(s) for Caspase 3: [6]
0 protease cleavage site(s) for Caspase 6: []
0 protease cleavage site(s) for Caspase 7: []
1 protease cleavage site(s) for Enterokinase: [20]
```

## Authors

This library was written by Kseniia Matveeva (@estakathersun) as a completing of 5th homework task of the *Python BI
2023* course.
> *P.S.: with tears and love* ️🥲❤️







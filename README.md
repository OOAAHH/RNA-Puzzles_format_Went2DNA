# RNA-Puzzles Format

Utilities for generating RNA/DNA reference PDB templates and checking submitted
PDB files against those templates.

## Files

- `rna_puzzles_format.py`
  Generate a standard-format reference PDB from a FASTA file.
- `format_check.py`
  Compare a submitted PDB file against a reference PDB template. If more than
  one format error is found, the script writes `<submission>.format_check.txt`.

## Requirements

- Python 3
- `Bio.PDB` from Biopython for `format_check.py`

## FASTA header format

`rna_puzzles_format.py` expects headers in this form:

```text
>name chain [DNA|RNA] [other optional fields]
```

Examples:

```text
>RNA1 A RNA
UGCGAUGAGAAGAAGAGUAUUAAGGAUUUACUAUGAUUAGCGACUCUAGGAUAGUGAAAG
CUAGAGGAUAGUAACCUUAAGAAGGCACUUCGAGCA

>DNA1 B DNA
ATCG
```

Notes:

- DNA and RNA are both supported.
- DNA output uses canonical residue names `DA/DC/DG/DT`.
- If the sequence contains `T` only, it is treated as DNA.
- If the sequence contains `U` only, it is treated as RNA.
- If the sequence is ambiguous (for example only `A/C/G`), declare `DNA` or
  `RNA` explicitly in the header.
- A single chain must not mix `T` and `U`.

## Usage

Generate a reference template:

```bash
python3 rna_puzzles_format.py example/2gdi.fa > out.pdb
python3 rna_puzzles_format.py example/2gdi.fa 1 > out.pdb
python3 rna_puzzles_format.py example/2gdi.fa 5 > out.pdb
python3 rna_puzzles_format.py example/dna.fa 1 > dna.out.pdb
```

The checked-in files under `example/` use `1` model to keep the repository
compact. If you omit the second argument, the script still writes `5` models.

Check a submitted structure:

```bash
python3 format_check.py example/2gdi.pdb example/2gdi.out.pdb
```

## Example files

See the `example/` directory for RNA and DNA examples:

- `example/2gdi.fa`
- `example/2gdi.pdb`
- `example/2gdi.out.pdb`
- `example/2gdi.pdb.format_check.txt`
- `example/dna.fa`
- `example/dna.out.pdb`

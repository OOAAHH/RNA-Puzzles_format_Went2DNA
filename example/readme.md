Example files for the current Python 3 RNA/DNA scripts.

The checked-in `.out.pdb` files use `1` model to keep the examples compact.
`rna_puzzles_format.py` still defaults to `5` models if you omit the second
argument.

- `2gdi.fa`
  RNA FASTA example with an explicit `RNA` declaration in the header.
- `2gdi.out.pdb`
  Reference template generated from `2gdi.fa` with one model:
  `python3 ../rna_puzzles_format.py 2gdi.fa 1 > 2gdi.out.pdb`
- `2gdi.pdb`
  Experimental PDB structure used as an example submission input.
- `2gdi.pdb.format_check.txt`
  Sample failure report produced by `format_check.py` for `2gdi.pdb` versus
  `2gdi.out.pdb`. Regenerating this file requires Biopython.
- `dna.fa`
  Minimal DNA FASTA example with an explicit `DNA` declaration in the header.
- `dna.out.pdb`
  Reference template generated from `dna.fa` with one model:
  `python3 ../rna_puzzles_format.py dna.fa 1 > dna.out.pdb`

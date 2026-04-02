#!/usr/bin/python

#===========================================================
#Copyright(c)2013, IBMC, CNRS
#All rights reserved.
#NAME:          rnatemplate.py
#ABSTRACT:      input a RNA/DNA sequence (fasta format), output the standard PDB format
#DATE:          Tue Sep 17 15:08:40 2013
#Usage:
#VERSION:       0.02
#AUTHOR:        Miao Zhichao
#CONTACT:       chichaumiau AT gmail DOT com
#NOTICE: This is free software and the source code is freely
#available. You are free to redistribute or modify under the
#conditions that (1) this notice is not removed or modified
#in any way and (2) any modified versions of the program are
#also available for free.
#              ** Absolutely no Warranty **
#===========================================================

import sys

Usage = """rnatemplate.py usage:

input a RNA/DNA sequence (fasta format), output the standard PDB format

./rnatemplate.py fasta.file number_of_model(optional) >output.pdb

Header format:
>name chain [DNA|RNA] [other optional fields]

Examples:
>RNA1 A RNA
UGCGAUGAGAAGAAGAGUAUUAAGGAUUUACUAUGAUUAGCGACUCUAGGAUAGUGAAAG
     CUAGAGGAUAGUAACCUUAAGAAGGCACUUCGAGCA
>DNA1 B DNA
ATCG
"""

RNA_RESIDUE_ATOMS = {
    'A': [ 'P', 'OP1', 'OP2', "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'", 'N9', 'C8', 'N7', 'C5', 'C6', 'N6', 'N1', 'C2', 'N3', 'C4' ],
    'G': [ 'P', 'OP1', 'OP2', "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'", 'N9', 'C8', 'N7', 'C5', 'C6', 'O6', 'N1', 'C2', 'N2', 'N3', 'C4' ],
    'U': [ 'P', 'OP1', 'OP2', "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'", 'N1', 'C2', 'O2', 'N3', 'C4', 'O4', 'C5', 'C6' ],
    'C': [ 'P', 'OP1', 'OP2', "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'", 'N1', 'C2', 'O2', 'N3', 'C4', 'N4', 'C5', 'C6' ],
}

DNA_RESIDUE_ATOMS = {
    'A': [ 'P', 'OP1', 'OP2', "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "C1'", 'N9', 'C8', 'N7', 'C5', 'C6', 'N6', 'N1', 'C2', 'N3', 'C4' ],
    'G': [ 'P', 'OP1', 'OP2', "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "C1'", 'N9', 'C8', 'N7', 'C5', 'C6', 'O6', 'N1', 'C2', 'N2', 'N3', 'C4' ],
    'T': [ 'P', 'OP1', 'OP2', "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "C1'", 'N1', 'C2', 'O2', 'N3', 'C4', 'O4', 'C5', 'C7', 'C6' ],
    'C': [ 'P', 'OP1', 'OP2', "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "C1'", 'N1', 'C2', 'O2', 'N3', 'C4', 'N4', 'C5', 'C6' ],
}


class TemplateInputError(Exception):
    pass


SUPPORTED_BASES = set('ACGTU')


def parse_declared_polymer_type(header_fields):
    for field in header_fields[2:]:
        normalized = field.strip().upper().rstrip(':;,')
        if normalized in ('DNA', 'RNA'):
            return normalized
        if normalized.startswith('TYPE=') or normalized.startswith('POLYMER='):
            _, value = normalized.split('=', 1)
            if value in ('DNA', 'RNA'):
                return value
    return None


def readfasta(fp):
    chains = []
    seqs = []
    polymer_types = []
    seq = ''
    with open(fp) as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if len(line) < 1:
                continue
            if line[0] == '#':
                continue
            if line[0] == '>':
                fields = line.split()
                if len(fields) < 2:
                    raise TemplateInputError('Invalid FASTA header: %s' % line)
                if len(chains) > len(seqs):
                    seqs.append(seq)
                    seq = ''
                chains.append(fields[1][0])
                polymer_types.append(parse_declared_polymer_type(fields))
                continue
            if not chains:
                raise TemplateInputError('Sequence data found before first FASTA header')
            seq += line.upper()
    if not chains:
        raise TemplateInputError('No FASTA records found')
    if len(chains) > len(seqs):
        seqs.append(seq)
    return chains, seqs, polymer_types


def infer_polymer_type(seq, chain, declared_type=None):
    invalid_bases = sorted(set(base for base in seq if base not in SUPPORTED_BASES))
    if invalid_bases:
        label = 'base' if len(invalid_bases) == 1 else 'bases'
        raise TemplateInputError(
            'Unsupported %s %s in chain %s'
            % (label, ','.join(invalid_bases), chain)
        )
    has_t = ('T' in seq)
    has_u = ('U' in seq)
    if has_t and has_u:
        raise TemplateInputError('Chain %s contains both T and U; cannot mix DNA and RNA bases in one chain' % chain)
    if declared_type == 'DNA':
        if has_u:
            raise TemplateInputError('Chain %s is declared as DNA but contains U' % chain)
        return 'DNA'
    if declared_type == 'RNA':
        if has_t:
            raise TemplateInputError('Chain %s is declared as RNA but contains T' % chain)
        return 'RNA'
    if has_t:
        return 'DNA'
    if has_u:
        return 'RNA'
    raise TemplateInputError(
        'Chain %s contains only A/C/G and is ambiguous; please declare DNA or RNA in the FASTA header, e.g. >name %s DNA'
        % (chain, chain)
    )


def get_residue_atoms(base, polymer_type):
    if polymer_type == 'DNA':
        atoms = DNA_RESIDUE_ATOMS.get(base)
    else:
        atoms = RNA_RESIDUE_ATOMS.get(base)
    if atoms is None:
        raise TemplateInputError('Unsupported base %s for %s chain' % (base, polymer_type))
    return atoms


def get_residue_name(base, polymer_type):
    if polymer_type == 'DNA':
        return 'D%s' % base
    return base


def format_atom_field(atom_name):
    if len(atom_name) >= 4:
        return atom_name[:4]
    return ' %-3s' % atom_name


def format_atom_line(serial, atom_name, residue_name, chain, residue_seq):
    element = atom_name[0]
    return (
        'ATOM  %5d %4s %3s %1s%4d       0.000   0.000   0.000  1.00  0.00           %1s\n'
        % (serial, format_atom_field(atom_name), residue_name, chain, residue_seq, element)
    )


def format_ter_line(serial, residue_name, chain, residue_seq):
    return 'TER   %5d      %3s %1s%4d                      \n' % (serial, residue_name, chain, residue_seq)


def prepare_model(chains, seqs, declared_polymer_types):
    serial = 1
    output = ''
    for chain, seq, declared_type in zip(chains, seqs, declared_polymer_types):
        if not seq:
            raise TemplateInputError('Chain %s has no sequence' % chain)
        polymer_type = infer_polymer_type(seq, chain, declared_type)
        residue_seq = 0
        for base in seq:
            residue_seq += 1
            residue_name = get_residue_name(base, polymer_type)
            for atom_name in get_residue_atoms(base, polymer_type):
                output += format_atom_line(serial, atom_name, residue_name, chain, residue_seq)
                serial += 1
        output += format_ter_line(serial, residue_name, chain, residue_seq)
        serial += 1
    return output


def format_pdb(fp, num=5):
    chains, seqs, declared_polymer_types = readfasta(fp)
    output = ''
    for model_index in range(num):
        output += 'MODEL       %2d                                              \n' % (model_index + 1)
        output += prepare_model(chains, seqs, declared_polymer_types)
        output += 'ENDMDL                                                      \n'
    output += 'END                                                        \n'
    sys.stdout.write(output + '\n')


if __name__ == '__main__':
    try:
        if len(sys.argv) < 2:
            sys.stdout.write(Usage)
            sys.exit(0)
        if len(sys.argv) > 2:
            format_pdb(sys.argv[1], int(sys.argv[2]))
        else:
            format_pdb(sys.argv[1])
    except TemplateInputError as error:
        sys.stderr.write('ERROR: %s\n' % error)
        sys.exit(1)

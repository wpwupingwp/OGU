#!/usr/bin/python3

import re
from functools import lru_cache

from Bio.Seq import Seq


@lru_cache(maxsize=None)
def gene_rename(old_name: str, genbank_format=False) -> (str, str):
    """
    Old doc:
        Different name of same gene will cause data to be splited to numerous
        files instead of one and some data may be dropped.

        For chloroplast genes, the author summarized various kinds of
        annotation error of gene name or synonyms and try to use regular
        expression to fix it.

        Ideally, use BLAST to re-annotate sequence is the best(and slow) way to
        find the correct name. This function only offers a "hotfix" which is
        enough.
    Rename plastid genes.
    May be dangerous.
    Will cache results.
    Args:
        old_name: old gene name
        genbank_format: use style like "trnH-GUG" or "trnHgug"
    Returns:
        new_name(str): new name, if fail, return old name
        gene_type(str): gene types, guessed from name
    """
    lower = old_name.lower()
    # (trna|trn(?=[b-z]))
    s = re.compile(r'(\d+\.?\d?)(s|rrn|rdna)')
    if lower.startswith('trn'):
        pattern = re.compile(r'([atcgu]{3})')
        prefix = 'trn'
        aa_letter = 'X'
        try:
            anticodon = Seq(re.search(pattern, lower[3:]).group(1))
        except Exception:
            return old_name, 'bad_name'
        # rna editing? trnI-CAU
        if anticodon == 'cau' and lower.startswith('trni'):
            aa_letter = 'I'
        # for trnfM-CAU
        elif lower.startswith('trnfm'):
            prefix = 'trnf'
            aa_letter = 'M'
        else:
            aa_letter = anticodon.reverse_complement().translate().upper()
            #anticodon = anticodon.transcribe()
        if genbank_format:
            new_name = f'{prefix}{aa_letter}-{anticodon.upper()}'
        else:
            new_name = f'{prefix}{aa_letter}{anticodon.lower()}'
        gene_type = 'tRNA'
    elif lower.startswith('rrn'):
        pattern = re.compile(r'(\d+\.?\d?)')
        try:
            number = re.search(pattern, lower).group(1)
        except Exception:
            return old_name, 'bad_name'
        new_name = 'rrn{}'.format(number)
        gene_type = 'rRNA'
    elif re.search(s, lower) is not None:
        new_name = 'rrn{}'.format(re.search(s, lower).group(1))
        gene_type = 'rRNA'
    else:
        pattern = re.compile(r'[^a-z]*'
                             '(?P<gene>[a-z]+)'
                             '[^a-z0-9]*'
                             '(?P<suffix>[a-z]|[0-9]+)')
        match = re.search(pattern, lower)
        try:
            gene = match.group('gene')
            suffix = match.group('suffix')
        except Exception:
            return old_name, 'bad_name'
        new_name = '{}{}'.format(gene, suffix.upper())
        # capitalize last letter
        if len(new_name) > 3:
            s = list(new_name)
            if s[-1].isalpha():
                new_name = '{}{}'.format(
                    ''.join(s[:-1]), ''.join(s[-1]).upper())
        gene_type = 'normal'
    if len(lower) >= 15:
        gene_type = 'suspicious_name'
    return new_name, gene_type


def plastid_rename():
    """
    Use name database.
    """
    pass


def safe_average(x):
    """
    Safe average.
    """
    if len(x) == 0:
        return 0
    else:
        return sum(x) / len(x)


def safe_path(old):
    """
    Remove illegal character in file path or name.
    """
    return re.sub(r'\W', '_', old)

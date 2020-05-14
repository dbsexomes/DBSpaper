"""
Normalization of variants.
"""


def normalize(pos, ref, alt):
    """Normalize a variant so that it is parsimonious and left-aligned

    >>> normalize(4, 'GCAT', 'GTGC')
    (5, 'CAT', 'TGC')
    >>> normalize(5, 'CATG', 'TGCG')
    (5, 'CAT', 'TGC')
    >>> normalize(4, 'GCATG', 'GTGCG')
    (5, 'CAT', 'TGC')
    >>> normalize(5, 'CAT', 'TGC')
    (5, 'CAT', 'TGC')
    >>> normalize(8, 'CA', '.')
    (8, 'CA', '.')
    >>> normalize(3, 'GCACA', 'GCA')
    (3, 'GCA', 'G')
    >>> normalize(2, 'GGCA', 'GG')
    (3, 'GCA', 'G')
    >>> normalize(4, 'TGAC', 'TTAC')
    (5, 'G', 'T')
    """

    if 1 == len(ref) == len(alt):
        return pos, ref, alt
    else:
        ref, alt = right_trim(ref, alt)
        return left_trim(pos, ref, alt)


def right_trim(r, a):
    """Given two strings `r` and `a` remove any common suffix
    between the characters

    >>> right_trim('AGCT', 'ACAT')
    ('AGC', 'ACA')
    >>> right_trim('AGAT', 'ACAT')
    ('AG', 'AC')
    >>> right_trim('TGAT', 'AGAT')
    ('T', 'A')
    >>> right_trim('AGGA', 'AA')
    ('AGG', 'A')
    >>> right_trim('AGGA', 'A')
    ('AGGA', 'A')
    >>> right_trim('AA', 'AGGA')
    ('A', 'AGG')
    >>> right_trim('A', 'T')
    ('A', 'T')
    """
    while r[-1] == a[-1]:
        if len(r) > 1 and len(a) > 1:
            r = r[:-1]
            a = a[:-1]
        else:
            break
    return r, a


def left_trim(pos, r, a):
    """Given two strings `r` and `a` remove any common prefix
    between the characters. The two strings passed to this
    should have *previously been trimmed* using right_trim

    >>> left_trim(1, 'A', 'G')
    (1, 'A', 'G')
    >>> left_trim(1, 'ATGC', 'AT')
    (2, 'TGC', 'T')
    >>> left_trim(1, 'AAT', 'AT')
    (1, 'AAT', 'AT')

    """
    if len(r) == len(a):
        clen = 1
    else:
        clen = 2
    while len(r) > 1 and len(a) > 1:
        if r[:clen] == a[:clen]:
            r = r[1:]
            a = a[1:]
            pos += 1
        else:
            break
    return pos, r, a


def primitives(pos, ref, alt):
    """Decompose a variant into it's primitives

    >>> primitives(25, 'ATGC', 'AACC')
    [(26, 'T', 'A'), (27, 'G', 'C')]
    """
    if len(ref) != len(alt):
        return [(pos, ref, alt)]
    p = []
    for idx, (r1, a1) in enumerate(zip(ref, alt)):
        if r1 != a1:
            p.append((pos + idx, r1, a1))
    return p


def variant_type(ref, alt):
    """Classify a variant as SNP, MNP, COMPLEX, INS or DEL

    >>> try:
    ...     variant_type('A', 'A')
    ... except ValueError:
    ...     pass
    >>> variant_type('A', 'C')
    'SNP'
    >>> variant_type('ATAC', 'ATTC')
    'SNP'
    >>> variant_type('ATA', 'GTG')
    'COMPLEX'
    >>> variant_type('ATG', 'GCT')
    'MNP'
    >>> variant_type('C', 'CA')
    'INS'
    >>> variant_type('CAGAG', 'CAG')
    'DEL'
    >>> variant_type('TGAC', 'GA')
    'COMPLEX'
    """
    if ref == alt:
        raise ValueError('Not a variant')

    _, ref, alt = normalize(0, ref, alt)
    nr = len(ref)
    na = len(alt)
    if nr == na:
        if nr == 1:
            vtype = 'SNP'
        else:
            distinct = sum(a != b for (a, b) in zip(ref, alt))
            if distinct == 1:  # cannot be reached
                vtype = 'SNP'
            elif distinct == nr:
                vtype = 'MNP'
            else:
                vtype = 'COMPLEX'
    else:
        if 'INS' in alt:
            vtype = 'INS'
        elif 'DEL' in alt:
            vtype = 'DEL'
        elif nr == 1:
            if na > 1 and (ref[0] == alt[0]):
                vtype = 'INS'
            else:
                vtype = 'COMPLEX'
        else:
            if na == 1 and (ref[0] == alt[0]):
                vtype = 'DEL'
            else:
                vtype = 'COMPLEX'
    return vtype


if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose=1)

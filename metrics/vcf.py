"""
.. module:: vcf
    :platform: Unix, Windows, MacOSX
    :synopsis: Parsing VCF 4.1 format files


Parser for VCF 4.1 format files. The format itself is documented_ at
the 1000genomes project site

.. _documented: http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41
"""
import sys
import os
from itertools import izip_longest as izip
import anyopen


# dictionary VCF data type to python functions
def toint(v):
    """Check and convert a value to an integer"""
    return int(v) if not v in '.AG' else v


def tofloat(v):
    """Check and convert a value to a floating point number"""
    return float(v) if v != '.' else None


def tostr(v):
    """Check and convert a value to a string"""
    return v if isinstance(v, str) else str(v)


def tochar(v):
    """Check and convert value to a character"""
    if isinstance(v, str) and len(v) == 1:
        return v
    else:
        raise ValueError('Expected a character, but found %s' % v)


TYPEMAP = {'float': tofloat,
           'integer': toint,
           'flag': bool,
           'string': str,  # short cut since we only have strings
           'character': tochar,
           None: None
          }

RTYPEMAP = {}
for key, value in TYPEMAP.items():
    if key is not None:
        RTYPEMAP[value] = key.capitalize()


class DotDict(dict):

    def __getattr__(self, v):
        return self.get(v, False)

    def __setattr__(self, v, item):
        self[v] = item



class VCFParser:
    """Parser for VCF 4.0 files"""

    def __init__(self, filename, sampleids=None, strict=False):
        self.filename = filename
        self.strict = strict
        self.stream = anyopen.openfile(filename, 'rt')
        self.meta = {'INFO': {},
                     'FORMAT': {},
                     'FILTER': {},
                     'ALT': {},
                     'SAMPLE': {},
                     'contig': []}
        self.samples = []
        self._headerlines = []
        self.parsemeta()
        if sampleids:
            self._sampleids = [self.samples.index(s) for s in sampleids]
        else:
            self._sampleids = list(range(len(self.samples)))

    def parsemeta(self):
        while 1:
            line = self.stream.readline()
            if line[:2] == '##':
                self._headerlines.append(line)
                self._parsemetarecord(line)
            elif line[:6].upper() == '#CHROM':
                self._parseheader(line)
                break

    def _parsemetarecord(self, line):
        """Parse meta information from the VCF headers.
        """
        key, value = line[2:].split('=', 1)
        if key not in self.meta:
            self.meta[key] = value.strip()
        else:
            d = self.meta[key]
            if isinstance(d, list):
                d.append(value)
            elif isinstance(d, dict):
                id, desc = self._parsemvalue(value.strip())
                d[id] = desc
            else:
                self.meta[key] = [d]
                self.meta[key].append(value)

    def _parsemvalue(self, v):
        if v[:1] != '<' or v[-1:] != '>':
            raise ValueError('Not a valid meta value:', v)

        v = v[1:-1]
        d = {}
        for item in self._splitw(v):
            if not item: continue
            key, value = item.split('=', 1)
            if key == 'Type':
                value = TYPEMAP[value.lower()]
            d[key] = value

        return d['ID'], (d.get('Number', None), d.get('Type', None),
                         d.get('Description', 'Unknown'))

    def _splitw(self, v, sep=','):
        words = []
        in_quote = False
        prevc = ''
        w = ''
        for c in v:
            if c == sep:
                if in_quote:
                    w += c
                else:
                    words.append(w)
                    w = ''
            elif c == '"' and prevc != '\\':
                if in_quote:
                    words.append(w)
                    w = ''
                    in_quote = False
                else:
                    in_quote = True
            else:
                w += c
            prevc = c
        if w.strip():
            words.append(w)
        return words

    def _parseheader(self, line):
        args = line[1:].split('\t')
        if len(args) < 8:
            raise ValueError('Header does not have all standard fields')
        args[-1] = args[-1].strip()
        self.samples = args[9:]

    def __next__(self):
        line = self.stream.readline()
        if not line:
            raise StopIteration
        return self.parseline(line)

    next = __next__  # for python 2.x

    def __iter__(self):
        return self

    def add_meta_info(self, id, number, type, description):
        self.meta['INFO'][id] = (str(number), TYPEMAP[type.lower()], \
                                 description)

    def add_meta(self, section, id, **kw):
        if not kw:
            self.meta[section] = id
        else:
            v = kw.get('Number', None)
            w = kw.get('Type', None)
            if w is not None:
                w = w.lower()
            x = kw.get('Description', 'Unknown Field')
            if section not in self.meta:
                self.meta[section] = {}
            self.meta[section][id] = (v, TYPEMAP[w], x)

    def parseline(self, line):
        rec = line.split('\t')
        rec[-1] = rec[-1].rstrip()
        vcfrec = DotDict()
        vcfrec.chrom = rec[0]
        vcfrec.pos = int(rec[1])
        vcfrec.id = [el.strip() for el in rec[2].split(';')]
        vcfrec.ref = rec[3].strip().upper()
        vcfrec.alt = rec[4].split(',')
        vcfrec.qual = tofloat(rec[5])
        if rec[6] != '.':
            vcfrec.filter = rec[6].split(';')
        else:
            vcfrec.filter = []
        vcfrec.info = rec[7]
        vcfrec.genotypes = rec[8:]
        vcfrec._numalt = rec[4].count(',') + 1
        vcfrec._line = line
        vcfrec._gtparsed = False
        vcfrec._infoparsed = False
        return vcfrec

    def parseinfo(self, vcfrec):
        if vcfrec._infoparsed:
            return
        dinfo = DotDict()
        rec = vcfrec.info
        if not isinstance(rec, str):
            return
        section = self.meta['INFO']
        for feat in rec.split(';'):
            # short cut for booleans. For non-boolean cases at this
            # point feat is of the form key=value, while for booleans
            # feat is just the key name. So, we can check if the feat
            # is present in the meta dictionary and if so, set the
            # value to true and move on
            if feat in section:
                dinfo[feat] = True
                continue
            try:
                key, value = feat.split('=', 1)
            except ValueError:
                dinfo[feat] = True
                continue
            try:
                n, Type, _ = section[key]
            except KeyError:
                if not self.strict:
                    n = '1'
                    Type = tostr
                    section[key] = ('1', tostr, 'Unknown field')
                else:
                    raise KeyError("Unknown INFO field %s" % key)
            if n == '1':
                value = Type(value)
            else:
                value = [Type(f) for f in value.split(',')]
                if n != '.' and self.strict:
                    if n == 'A':
                        n = vcfrec['_numalt']
                    elif n == 'G':
                        n = (vcfrec['_numalt'] * (vcfrec['_numalt'] + 1)) / 2
                    else:
                        n = int(n)
                    if len(value) != n:
                        raise ValueError("Need %d values for %s. Got %d" %
                                         ((n, key, len(value))))
            dinfo[key] = value

        vcfrec.info = dinfo
        vcfrec._infoparsed = True

    def parsegenotypes(self, vcfrec, _izip=izip):
        if vcfrec._gtparsed:
            return
        rec = vcfrec.genotypes
        fields = rec.pop(0).split(':')
        if len(rec) != len(self.samples):
            raise ValueError('%s:\nGenotype data not available for all ' + \
                             'samples' % self.rec)
        section = self.meta['FORMAT']
        for idx in self._sampleids:
            vcfrec[self.samples[idx]] = gt = DotDict()
            for fname, val in _izip(fields, rec[idx].split(':')):
                if val is None:
                    gt[fname] = val
                    continue
                n, Type, _ = section[fname]
                if n == '1':
                    gt[fname] = Type(val)
                    continue
                val = [Type(f) for f in val.split(',')]
                if n != '.' and self.strict:
                    if n == 'G':
                        ntot = vcfrec['_numalt']
                        n = (ntot * (ntot + 1)) // 2
                    elif n == 'A':
                        n = vcfrec['_numalt']
                    else:
                        n = int(n)
                    if n != len(val):
                        raise ValueError(
                                'Expected %s values for %s. Got %d' %
                                (n, fname, len(val)))
                gt[fname] = val
        vcfrec._gtparsed = True

    def writeheader(self, stream, original=False):
        metad = {'INFO': ('ID', 'Number', 'Type', 'Description'),
                 'FORMAT': ('ID', 'Number', 'Type', 'Description'),
                 'FILTER': ('ID', 'Name'),
                }

        if original:
            stream.writelines(self._headerlines)
        else:
            stream.write('##fileformat=%s\n' % self.meta.get('fileformat',
                                                             'VCFv4.1'))
            for key in self.meta.keys():
                if key == 'fileformat': continue
                v = self.meta[key]
                if isinstance(v, str):
                    stream.write('##%s=%s%s' % (key, v, os.linesep))
                elif isinstance(v, list):
                    for item in v:
                        stream.write('##%s=%s%s' % (key, item.rstrip(),
                                                    os.linesep))
                elif isinstance(v, dict):
                    for key2 in v:
                        vstr = ['##%s=<ID=%s' % (key, key2)]
                        Number, Type, Descr = v[key2]
                        if Number is not None:
                            vstr.append('Number=%s' % Number)
                        if Type is not None:
                            if not isinstance(Type, str):
                                Type = RTYPEMAP[Type]
                            vstr.append('Type=%s' % Type)
                        if Descr is not None:
                            vstr.append('Description="%s">%s' % (Descr,
                                                                 os.linesep))
                        stream.write(','.join(vstr))
                else:
                    raise TypeError('Unknown meta data type %s, %s' % (key, v))

        h = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO',
             'FORMAT']
        for idx in self._sampleids:
            h.append(self.samples[idx])
        stream.write('\t'.join(h))
        stream.write(os.linesep)

    def _writeinfoheader(self, stream):
        INFO = self.meta['INFO']
        for key in INFO:
            number, typ, descr = INFO[key]
            typ = RTYPEMAP[typ]
            stream.write('##INFO=<ID=%s,Number=%s,Type=%s,Description="%s">\n'\
                          % (key, number, typ, descr))

    def write(self, stream, rec, original=False):
        if original:
            stream.write(rec['_line'])
            return

        if not rec.filter:
            rec.filter = ['.']

        args = [rec.chrom, str(rec.pos), ';'.join(rec.id),
                rec.ref, ','.join(rec.alt), str(rec.qual or '.'),
                ';'.join(rec.filter)]
        args.append(self._getinfostring(rec))
        args.extend(self._getgtstring(rec))
        stream.write('\t'.join(args))
        stream.write(os.linesep)

    def _getinfostring(self, rec):
        if not rec._infoparsed:
            return rec.info
        else:
            section = self.meta['INFO']
            infargs = list(rec.info.keys())
            infargs.sort()
            atmp = []
            for key in infargs:
                n, _, _ = section[key]
                if n == '0':
                    atmp.append(key)
                else:
                    if n == '1':
                        v = rec.info[key]
                        if v is None:
                            v = '.'
                        atmp.append('%s=%s' % (key, v))
                    else:
                        output = []
                        for el in rec.info[key]:
                            if el is None:
                                output.append('.')
                            else:
                                output.append(str(el))
                        atmp.append('%s=%s' % (key, ','.join(output)))

            return ';'.join(atmp)

    def _getinfolastpos(self):
        'Return the index of the last INFO field in the header list'
        infolist = []
        for indx, item in enumerate(self._headerlines):
            if item.startswith('##INFO'):
                infolist.append(item)
        return self._headerlines.index(infolist[-1])

    def insert_infoheader(self, el):
        'Insert an element in the header list of type info'
        lastpos = self._getinfolastpos()
        self._headerlines.insert(lastpos + 1, el)

    def _getgtstring(self, rec):
        if not rec._gtparsed:
            return rec.genotypes
        first = 1
        args = []
        for idx in self._sampleids:
            v = rec.get(self.samples[idx])
            if first:
                formats = list(v.keys())
                formats.sort()
                if 'GT' in formats:
                    formats.remove('GT')
                    formats.insert(0, 'GT')
                args.append(':'.join(formats))
                first = 0
            ftmp = []
            section = self.meta['FORMAT']
            for key in formats:
                n, _, _ = section[key]
                item = v[key]
                if item is None:
                    if key == 'GT':
                        item = './.'
                    else:
                        item = '.'
                if n == '1':
                    ftmp.append(str(item))
                else:
                    ftmp.append(','.join(str(el) for el in item))
            if ftmp[0] in ('.', './.'):
                args.append(ftmp[0])
            else:
                args.append(':'.join(ftmp))

        return args

    def add_info(self, vcfrec, el, val):
        vcfrec['info'][el] = val


def run():
    import sys
    from time import time
    n = 0
    t0 = time()
    parser = VCFParser(sys.argv[1])
    for rec in parser.iterate():
        parser.parseinfo()
        n += 1
        if not n % 10000:
            print('Read %d records in %g seconds' % (n, time() - t0))
    print('Total of %d records parsed' % n)
    print(rec)


if __name__ == "__main__":
    run()


"""
.. module:: anyopen
    :platform: Unix, Windows, MacOSX
    :synopsis: Transparent opening of compressed and uncompressed files


Helper to open files of various types.  This module defines a function
*openfile* which is passed the name of the file and returns a file object
opened for reading.  The function can open files compressed using gzip,
bzip2, unix compress in addition to plain text uncompressed files.
"""
"""
.. module:: anyopen
    :platform: Unix, Windows, MacOSX
    :synopsis: Transparent opening of compressed and uncompressed files


Helper to open files of various types.  This module defines a function
*openfile* which is passed the name of the file and returns a file object
opened for reading.  The function can open files compressed using gzip,
bzip2, unix compress in addition to plain text uncompressed files.
"""
import gzip
import bz2
try:
    import lzma
except ImportError:
    lzma = None
import os
import tempfile
import subprocess
import sys
from urllib2 import urlopen


def openfile(filename, mode='rt'):
    """Open a file for reading or writing

    Many of the files that are used in bioinformatics are compressed due
    of their huge size. Some are compressed using gzip and others using
    unix compress. It becomes tedious for programs that process the data in
    these files to check for the file format and then use the appropriate
    method to open the file.  Using this function eliminates the need for
    such programs to check for the file format. The function handles files
    compressed using gzip, bzip2 and unix compress in addition to uncompressed
    files. As a bonus files can also be opened for reading using http or
    ftp protocols. Finally specifying the filename as '-' returns stdin or
    stdout depending on whether the mode is 'r' or 'w' respectively

    Args:
        filename (str): Name or URL of the file to be opened
        mode (str): The file open mode 'rb', 'wb', 'rt', etc.

    Returns:
        A file object opened for reading or writing

    Examples

    >>> infile = openfile('/var/data/uniprot/uniprot_sprot.dat.gz')
    >>> for line in infile:
    ...     pass

    >>> infile = openfile('/var/data/enzyme/enzyme.dat.Z')
    >>> for line in infile:
    ...     pass

    """
    func = {'.gz': gzip.open,
            '.bz2': bz2.BZ2File,
            '.Z': _open_zcat}
    if lzma:
        func['.xz'] = lzma.LZMAFile

    if not isinstance(filename, str):
        if 'r' in mode and hasattr(filename, 'readline'):
            return filename
        if 'w' in mode and hasattr(filename, 'write'):
            return filename
        raise TypeError('Expected a filename or an already opened stream')

    if filename == '-':
        if 'r' in mode:
            return sys.stdin
        else:
            return sys.stdout

    if 'w' in mode or 'a' in mode or 'x' in mode:
        opener = func.get(os.path.splitext(filename)[1], open)
        tmpname = filename
    else:
        if filename.startswith('http://') or filename.startswith('ftp://') or \
           filename.startswith('file://'):
            tmpname = __geturl(filename)
        else:
            tmpname = filename
        fobj = open(tmpname, 'rb')
        data = fobj.read(6)
        fobj.close()

        if data[:2] == b'\037\213':
            opener = func['.gz']
        elif data[:3] == b'BZh':
            opener = func['.bz2']
        elif data[:2] == b'\037\235':
            opener = func['.Z']
        elif data == b'\xfd7zXZ\x00':
            try:
                opener = func['.xz']
            except KeyError:
                raise NotImplementedError('Module lzma not available')
        else:
            opener = open

    stream = opener(tmpname, mode)

    return stream


def _open_cat(progname, filename):
    return subprocess.Popen('%s %s' % (progname, filename), shell=True,
                            stdout=subprocess.PIPE, close_fds=True,
                            universal_newlines=True).stdout


def _open_zcat(filename):
    return _open_cat('zcat', filename)


def readfile(filename, splitlines=True):
    """Read data from a file.

    This function reads all data from

        -- file compressed using gzip, bzip2 or unix compress,
        -- uncompressed file

    Args:
        filename (str): Name of file or URL
        splitlines (bool): If specified and True then the data is split into
           lines

    Returns:
        A string containing all the data if *splitlines* is False (the default)
        or a list of lines if *splitlines* is True.

    """

    tmpname = filename
    if isinstance(filename, str):
        if filename.startswith('http://') or filename.startswith('ftp://') or \
           filename.startswith('file://'):
            tmpname = __geturl(filename)

    fobj = openfile(tmpname, 'rt')
    data = fobj.stream.read()
    fobj.stream.close()

    if tmpname != filename:
        os.unlink(tmpname)

    if splitlines:
        return data.splitlines(True)
    else:
        return data


def __geturl(url):
    # private method to fetch a file given an URL

    fd, tmpname = tempfile.mkstemp()
    fobj = urlopen(url)
    while 1:
        data = fobj.read(1024000)
        if not data:
            break
        os.write(fd, data)
    os.close(fd)
    fobj.close()
    return tmpname

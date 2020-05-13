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
from urllib.request import urlopen


def openfile(filename, mode='rt', encoding=None):
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
        encoding (str): Specify encoding for file (optional)

    Returns:
        A file object opened for reading or writing

    Examples

    >>> infile = openfile('/var/data/uniprot/uniprot_sprot.dat.gz')
    >>> for line in infile:
    ...     pass

    >>> infile = openfile('/var/data/enzyme/enzyme.dat.Z')
    >>> for line in infile:
    ...     pass

    Notes:

    In python 3.x opening a file in binary mode will return a byte stream.
    To convert this to text, call decode on the returned byte string.
    Alternatively passing the mode string at 'rt', the default, to the openfile
    function will produce a text stream as output.
    """
    func = {'.gz': gzip.open,
            '.bz2': bz2.BZ2File,
            '.Z': _open_zcat}
    if lzma:
        func['.xz'] = lzma.open

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

    if 'b' in mode:
        stream = opener(tmpname, mode)
    else:
        stream = opener(tmpname, mode, encoding=encoding)

    return stream


def _open_cat(progname, filename, binary):
    un = not binary
    return subprocess.Popen('%s %s' % (progname, filename), shell=True,
                            stdout=subprocess.PIPE, close_fds=True,
                            universal_newlines=un).stdout


def _open_zcat(filename, binary):
    return _open_cat('zcat', filename, binary)


def readfile(filename, splitlines=True, encoding=None):
    """Read data from a file.

    This function reads all data from

        -- file compressed using gzip, bzip2 or unix compress,
        -- uncompressed file

    Args:
        filename (str): Name of file or URL
        splitlines (bool): If specified and True then the data is split into
           lines
        encoding (str): Optional argument specifying the encoding of the file

    Returns:
        A string containing all the data if *splitlines* is False (the default)
        or a list of lines if *splitlines* is True.

    """

    tmpname = filename
    if isinstance(filename, str):
        if filename.startswith('http://') or filename.startswith('ftp://') or \
           filename.startswith('file://'):
            tmpname = __geturl(filename)

    with openfile(tmpname, 'rt', encoding=encoding) as stream:
        data = stream.read()

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

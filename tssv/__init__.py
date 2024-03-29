from argparse import FileType
from os.path import exists
from pkg_resources import get_distribution

from .align_pair import align_pair
from .sg_align import align


class ProtectedFileType(FileType):
    def __call__(self, string):
        if 'w' in self._mode and exists(string):
            raise IOError('failed to create "{}": file exists.'.format(string))
        return super(ProtectedFileType, self).__call__(string)


def _get_metadata(name):
    pkg = get_distribution('tssv')

    for line in pkg.get_metadata_lines(pkg.PKG_INFO):
        if line.startswith('{}: '.format(name)):
            return line.split(': ')[1]

    return ''


_copyright_notice = 'Copyright (c) {} <{}>'.format(
    _get_metadata('Author'), _get_metadata('Author-email'))

usage = [_get_metadata('Summary'), _copyright_notice]


def doc_split(func):
    return func.__doc__.split('\n\n')[0]


def version(name):
    return '{} version {}\n\n{}\nHomepage: {}'.format(
        _get_metadata('Name'), _get_metadata('Version'), _copyright_notice,
        _get_metadata('Home-page'))

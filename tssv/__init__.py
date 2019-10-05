from argparse import FileType
from os.path import abspath, dirname, exists

from configparser import ConfigParser

from .align_pair import align_pair
from .sg_align import align


config = ConfigParser()
with open('{}/setup.cfg'.format(dirname(abspath(__file__)))) as handle:
    config.read_file(handle)

_copyright_notice = 'Copyright (c) {} {} <{}>'.format(
    config.get('metadata', 'copyright'),
    config.get('metadata', 'author'),
    config.get('metadata', 'author_email'))

usage = [config.get('metadata', 'description'), _copyright_notice]


class ProtectedFileType(FileType):
    def __call__(self, string):
        if 'w' in self._mode and exists(string):
            raise IOError('failed to create "{}": file exists.'.format(string))
        return super(ProtectedFileType, self).__call__(string)


def doc_split(func):
    return func.__doc__.split('\n\n')[0]


def version(name):
    return '{} version {}\n\n{}\nHomepage: {}'.format(
        config.get('metadata', 'name'),
        config.get('metadata', 'version'),
        _copyright_notice,
        config.get('metadata', 'url'))

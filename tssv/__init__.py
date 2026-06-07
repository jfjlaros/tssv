from argparse import FileType
from importlib.metadata import PackageNotFoundError, metadata
from os.path import exists
from re import split
from typing import Callable

from .align_pair import align_pair
from .sg_align import align


class ProtectedFileType(FileType):
    def __call__(self, string):
        if 'w' in self._mode and exists(string):
            raise IOError('failed to create "{}": file exists.'.format(string))
        return super(ProtectedFileType, self).__call__(string)


def _extract(key: str, delim: str = r'[^\s\S]', index: int = 0) -> str:
    try:
        value = metadata(__package__).get(key, '')
    except PackageNotFoundError:
        return '<NO DATA>'
    return split(delim, value)[index]


def doc_split(func: Callable) -> str:
    return func.__doc__.split('\n\n')[0]


_project = _extract('Name')
_version = _extract('Version')
_year = '2010-2026'
_author = _extract('Author-email', r'"', 1)
_email = _extract('Author-email', r'<|>', 1)
_description = _extract('Summary')
_copyright = f'Copyright (c) {_year} by {_author} <{_email}>'
_url = _extract('Project-URL', r', ', 1)
_info = f'{_project} version {_version}\n\n{_copyright}\nHomepage: {_url}'

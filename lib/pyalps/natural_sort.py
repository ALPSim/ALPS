#  Copyright Haruhiko Matsuo 2011.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

import doctest
import re

def natural_sort(str_list, lower_case=False):
    """
    Sort a given string list in natural order,

    >>> natural_sort(['file2b', 'file11A', 'file2B', 'file1a'])
    ['file1a', 'file2B', 'file2b', 'file11A']

    Sort lower-case to upper-case,

    >>> natural_sort(['file2b', 'file11A', 'file2B', 'file1a'], True)
    ['file1a', 'file2b', 'file2B', 'file11A']
    """
    if lower_case:
        atoi = lambda text: int(text) if text.isdigit() else text.lower()
    else:
        atoi = lambda text: int(text) if text.isdigit() else text

    natural_key = lambda key: [atoi(s) for s in re.split('([0-9]+)',key)]
    return sorted(str_list, key=natural_key)

if __name__ == '__main__':
    doctest.testmod()

from __future__ import absolute_import
# ****************************************************************************
#
# ALPS Project: Algorithms and Libraries for Physics Simulations
#
# ALPS Libraries
#
# Copyright (C) 2016 by Michele Dolfi <dolfim@phys.ethz.ch>
#
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the “Software”),
# to deal in the Software without restriction, including without limitation
# the rights to use, copy, modify, merge, publish, distribute, sublicense,
# and/or sell copies of the Software, and to permit persons to whom the
# Software is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included
# in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS
# OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.
#
# ****************************************************************************


## The purpose of this script is to import the compiled modules both in the
## installed directory (relative to the script) and in the build dicrectory
## while testing (absolute modules, available via PYTHONPATH)

try:
    from . import pyalea_c
    from . import pymcdata_c
    from . import pyngsapi_c
    from . import pyngsbase_c
    from . import pyngshdf5_c
    from . import pyngsobservable_c
    from . import pyngsobservables_c
    from . import pyngsparams_c
    from . import pyngsrandom01_c
    from . import pyngsresult_c
    from . import pyngsresults_c
    from . import pytools_c
except ImportError:
    import pyalea_c
    import pymcdata_c
    import pyngsapi_c
    import pyngsbase_c
    import pyngshdf5_c
    import pyngsobservable_c
    import pyngsobservables_c
    import pyngsparams_c
    import pyngsrandom01_c
    import pyngsresult_c
    import pyngsresults_c
    import pytools_c

import sys
for k in list(locals().keys()):
    if k.endswith('_c'):
        sys.modules["{}.cxx.{}".format(__package__, k)] = locals()[k]

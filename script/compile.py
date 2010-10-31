#! python

import compileall
compileall.compile_dir("lib", force=1)
compileall.compile_dir("vistrails", force=1)

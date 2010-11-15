import sys, os, glob, subprocess

# The python interpreter running this script
pycmd = sys.executable

# Store current working directory
cwd = os.getcwd()

# Find .py files in subdirectories
pyfiles = glob.glob('[d-z]*/*.py')

# Test all Python tutorials
logfile = open('pytest.log', 'w')
for f in pyfiles:
    if f.find('build_lattice.py') >= 0: # This is no tutorial and needs cmd line arguments
        continue
    print f + ':',
    dir = os.path.dirname(f)
    fn = os.path.basename(f)
    os.chdir(os.path.join(cwd,dir))
    process = subprocess.Popen([pycmd,fn], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output = process.communicate()
    logfile.write(f+':'+'\n')
    errors = []
    for k in output:
        logfile.write(k)
    logfile.write('===============================================================\n')
    print 'returned',process.returncode

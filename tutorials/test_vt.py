import sys, os, subprocess, glob

vtapp = '/Applications/VisTrails/Vistrails.app/Contents/MacOS/vistrails'

# Find .vt files
vtfiles = glob.glob('*/*.vt')

# Extract workflow tags from vt files
workflows = []
for vt in vtfiles:
    cmd = ['unzip', '-c', vt, 'vistrail']
    xmltrail = subprocess.Popen(cmd, stdout=subprocess.PIPE).communicate()[0]
    for line in xmltrail.splitlines():
        if line.find('key="__tag__"') == -1:
            continue
        tagstart = line.find('value="')+len('value="')
        tagend = line.find('"', tagstart+1)
        tag = line[tagstart:tagend]
        if tag != 'cannot prune':   # this seems to be some auto-generated tag
            workflows.append( (vt,tag) )
            #print os.path.basename(vt) + ':"' + tag + '"'
    
# Test all tagged workflows
logfile = open('vttest.log', 'w')
for workflow in workflows:
    (fn,tag) = workflow
    fn = os.path.join(os.getcwd(),fn)
    if not os.path.exists(fn):
        print fn,'does not exist!'
    print fn + ':"' + tag + '" ',
    cmd = [vtapp, '-b', fn+':'+tag]
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output = process.communicate()
    logfile.write(fn+':'+tag+'\n')
    for k in output:
        logfile.write(k)
    logfile.write('===============================================================\n')
    print 'returned',process.returncode


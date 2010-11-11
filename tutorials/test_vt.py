import sys, os, subprocess

workflows = [('mc-02-susceptibilities/mc-02-susceptibilities.vt','quantum_ladder'),('ed-02-gaps/ed-02-gaps.vt','Flexible_system_sizes')]
vt = '/Applications/VisTrails/Vistrails.app/Contents/MacOS/vistrails'

logfile = open('log.log', 'w')
for workflow in workflows:
    fn = os.path.join(os.getcwd(), workflow[0])
    if not os.path.exists(fn):
        print fn,'does not exist!'
    cmd = [vt, '-b', fn+':'+workflow[1]]
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output = process.communicate()
    logfile.write(workflow[0]+':'+workflow[1]+'\n')
    for k in output:
        logfile.write(k)
    print 'Return code of',workflow,process.returncode


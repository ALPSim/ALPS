import sys, shutil, tempfile
import os.path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pyalps
import pyalps.plot
from subprocess import check_call


basename = 'sim_a'
data = pyalps.loadMPSIterations(pyalps.getResultFiles(prefix=basename), what=['Local Magnetization'])

dirname = tempfile.mkdtemp()
print '..temp dir:', dirname

data = pyalps.flatten(data)
ymax = 0.
for d in data:
    d.y = d.y[0]
    ymax = max(ymax, max(d.y))
    # formatting
    d.props['line']  = '-x'
    d.props['color'] = 'b'

counter = 0
for d in sorted(data, key = lambda d: d.props['iter_t']):
    sw = d.props['iter_t']
    t = d.props['dt'] * (sw+1)
    print '..plot iteration', int(sw)
    # print group
    plt.figure()
    plt.title('$t = %s$' % t)
    pyalps.plot.plot([d])
    plt.xlabel('x')
    plt.ylabel('Local magnetization')

    plt.ylim(-1, +1)
    
    plt.savefig(dirname+'/_anim.%08d.png' % counter)
    plt.close()
    counter += 1

oname = basename +'.density.mp4'
print '..generating', oname
check_call(['ffmpeg',
            '-i', dirname+'/_anim.%08d.png',
            '-c:v','libx264','-profile:v','high','-crf','23','-pix_fmt','yuv420p','-y','-r','30',
            oname])

shutil.rmtree(dirname)
    
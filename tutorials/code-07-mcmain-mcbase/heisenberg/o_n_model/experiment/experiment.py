import subprocess32 as sp
import sys, os, argparse, pyalps, re
from matplotlib import pyplot as plt
import numpy as np

class commandline_interface(object):
    def __init__(self):
        self.parser = argparse.ArgumentParser(description='Run and analyze heisenberg simulation experiment')
        self.parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.0')
        self.subparsers = self.parser.add_subparsers(title='subcommands')

        self.sp_run = self.subparsers.add_parser('run', help='run the complete experiment')
        self.sp_analyze = self.subparsers.add_parser('analyze', help='''run just the analysis step on already 
                                                     created simulation results''')
        self.sp_clean = self.subparsers.add_parser('clean', help='''remove all intermediary and result files for the given
                                                   input file''')


        self.sp_run.add_argument('-j', metavar='N', help='run N simulations in parallel (requires whe GNU parallel program.)')
        self.sp_run.add_argument('--no-analysis', action='store_true', help='only run the simulations, do not analyze results')
        self.sp_run.add_argument('program', help='the simulation executable')
        self.sp_run.add_argument('infile', default='param.txt')

        self.sp_analyze.add_argument('infile', default='param.txt')
        
        self.sp_clean.add_argument('infile', default='param.txt', metavar='infile')

    def add_run_func(self, run_func):
        self.sp_run.set_defaults(func=run_func)

    def add_analyze_func(self, an_func):
        self.sp_analyze.set_defaults(func=an_func)

    def add_clean_func(self, cl_func):
        self.sp_clean.set_defaults(func=cl_func)

def run_exp(args):
    '''create individual input files and run our simulation on them.
    '''
    sp.check_output(['parameter2xml', args.infile])         # create an xml input file for each run
    runs = [f for f in os.listdir(os.getcwd()) if re.match('{}.task\d*.in.xml'.format(args.infile), f)]
    if args.j:                                      # in case the -j option is specified we want to use parallel
        prog_call = ['parallel', '-j', args.j, args.program, ':::']
        prog_call.extend(runs)
        run_out = sp.check_output(prog_call)
    else:                                           # good old sequential way of doing things
        for run_infile in runs:
            run_out = sp.check_output([args.program, run_infile])
            print '--- run {} ---'.format(run_infile)
            print run_out
    if not args.no_analysis:                        # run analysis if not otherwise specified
        analyze(args)

def get_vector_data(result_files, what):
    '''read data from a VectorObservable into a multidimensional numpy array,
    discarding all info except "mean" and "error".
    '''
    from pyalps.alea import MCScalarData as msd
    data_set = pyalps.loadMeasurements(result_files, what=what)
    return np.array([[msd(i.mean, i.error) for i in j[0].y] for j in data_set])

def get_scalar_data(result_files, what):
    '''read data from a ScalarObservable into a numpy array,
    discarding all info except "mean" and "error"
    '''
    from pyalps.alea import MCScalarData as msd
    data_set = pyalps.loadMeasurements(result_files, what=what)
    return np.array([msd(i[0].y[0].mean, i[0].y[0].error) for i in data_set])

def get_vector_mean(result_files, what):
    '''read only means from a VectorObservable into a multidimensional numpy array.
    '''
    data_set = pyalps.loadMeasurements(result_files, what=what)
    return np.array([i[0].y.mean for i in data_set])

def get_chi(result_files):
    '''collect information and perform calculations for the magnetic susceptibility plot.
    '''
    T   = np.array([p['T'] for p in pyalps.loadProperties(result_files)])
    V   = np.array([p['L']**3 for p in pyalps.loadProperties(result_files)])
    mag  = get_vector_data(result_files, what='Magnetization')
    mm   = np.sum(np.vectorize(pyalps.math.sq)(mag), 1)
    m2   = get_scalar_data(result_files, what='Magnetization^2')

    beta = 1. / T
    chi = beta * V * (m2 - mm)

    chi_dt = [('T', 'f8'), ('chi', 'object')]
    result = np.array(zip(T, chi), dtype = chi_dt)
    result.sort(order='T')
    return result

def get_corr(result_files):
    '''collect information and perform calculations for the correlation function plot.
    '''
    T    = np.array([p['T'] for p in pyalps.loadProperties(result_files)])
    corr = get_vector_data(result_files, what = 'Correlations')
    dist = get_vector_mean(result_files, what = 'Distances')
    dist_corr_dt = [('dist', 'f8'), ('corr', 'O')]

    cbins = [{i:[] for i in np.unique(j)} for j in dist]
    result = [i for i in range(len(result_files))]
    for i in range(len(result_files)):                                      # for each result file
        corr_data = np.array(zip(dist[i], corr[i]), dtype = dist_corr_dt)   # group distance and correlation together
        for cj in corr_data:
            cbins[i][cj['dist']].append(cj['corr'])                         # group all correlations belongning to the same distance
        for j in cbins[i]:
            cbins[i][j] = np.abs( np.sum(cbins[i][j]) / float(len(cbins[i][j])) ) # take the mean for every distance

        # group distance, mean correlation into a numpy array 
        result[i] = np.array([(d, c) for (d, c) in cbins[i].iteritems()], dtype = dist_corr_dt)
        result[i].sort(order='dist')                                        # sort in order of ascending distance
    return result, T

def analyze(args):
    print 'running analysis'
    runs = pyalps.getResultFiles(prefix = args.infile)

    chi_data = get_chi(runs)
    T = chi_data['T']
    chi = [c.mean for c in chi_data['chi']]
    chi_err = [c.error for c in chi_data['chi']]

    plt.figure()
    plt.errorbar(T, chi, yerr = chi_err)
    yl = plt.ylim()
    plt.vlines(1. / 0.693035, yl[0], yl[1])
    plt.ylim(yl)
    xl =(min(T), max(T))
    plt.xlim(xl)
    plt.xlabel(r'Temperature $T$')
    plt.ylabel(r'Magnetic Susceptibility $\chi$')
    plt.savefig('{}.plot_chi.pdf'.format(args.infile))

    idx = [1, len(runs)/3, -2]
    T = np.array([p['T'] for p in pyalps.loadProperties(runs)])
    runsorted = np.array(zip(T, runs), dtype=[('T', 'f8'), ('run', object)])
    runsorted.sort(order = 'T')
    corr_data, T = get_corr([runsorted[i]['run'] for i in idx])

    plt.figure()
    for i in range(len(corr_data)):
        temp = T[i]
        corr = corr_data[i]
        x = corr['dist']
        y = [i.mean for i in corr['corr']]
        ye = [i.error for i in corr['corr']]
        plt.errorbar(x, y, yerr = ye, label = 'T = {}'.format(temp))
    plt.ylim((-0.01, 1.01))
    plt.xlabel(r'Distance $r$ in Lattice Units')
    plt.ylabel(r'Pair Correlation Function $g(r)$')
    plt.legend()
    plt.savefig('{}.plot_corr.pdf'.format(args.infile))

def clean(args):
    resfiles = pyalps.getResultFiles(prefix = args.infile)
    for resf in resfiles:
        os.remove(resf)
    infiles = [f for f in os.listdir(os.getcwd()) if re.match('{}.*in.*'.format(args.infile), f)]
    for inf in infiles:
        os.remove(inf)

cli = commandline_interface()
cli.add_run_func(run_exp)
cli.add_analyze_func(analyze)
cli.add_clean_func(clean)

if __name__ == '__main__':
    args = cli.parser.parse_args()
    args.func(args)

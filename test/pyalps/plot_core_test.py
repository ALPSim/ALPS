# ALPS Project: https://alps.comp-phys.org/
# SPDX-License-Identifier: MIT
from pyalps.dataset import DataSet
from pyalps.floatwitherror import FloatWithError
from pyalps.plot_core import convertToText, makeGracePlot, makeGnuplotPlot


def dataset(x, y, label=None):
    d = DataSet()
    d.x = x
    d.y = y
    if label is not None:
        d.props['label'] = label
    return d


def data_lines(output):
    return [l for l in output.splitlines() if l and l[0] not in '@#&e' and not l.startswith('plot')]


def test_multi_point_errors():
    # Error arrays with more than one element must not be compared with == None.
    d = dataset([1.0, 2.0], [FloatWithError(1.0, 0.1), FloatWithError(2.0, 0.2)])
    grace = makeGracePlot([d])
    assert '@type xydy' in grace
    assert data_lines(grace) == ['1.0\t1.0\t0.1', '2.0\t2.0\t0.2']
    gnuplot = makeGnuplotPlot([d])
    assert '# X Y DY' in gnuplot
    assert data_lines(gnuplot) == ['1.0\t1.0\t0.1', '2.0\t2.0\t0.2']


def test_grace_x_errors():
    x = [FloatWithError(1.0, 0.1), FloatWithError(2.0, 0.2)]
    grace = makeGracePlot([dataset(x, [10.0, 20.0])])
    assert '@type xydx' in grace
    assert data_lines(grace) == ['1.0\t10.0\t0.1', '2.0\t20.0\t0.2']
    grace = makeGracePlot([dataset(x, [FloatWithError(10.0, 1.0), FloatWithError(20.0, 2.0)])])
    assert '@type xydxdy' in grace
    assert data_lines(grace) == ['1.0\t10.0\t0.1\t1.0', '2.0\t20.0\t0.2\t2.0']


def test_gnuplot_errors_per_dataset():
    # Each dataset is written in its own format, not the last dataset's.
    with_errors = dataset([1.0], [FloatWithError(5.0, 0.5)], 'with errors')
    plain = dataset([2.0], [7.0], 'plain')
    gnuplot = makeGnuplotPlot([with_errors, plain])
    assert '# X Y DY \n1.0\t5.0\t0.5\n' in gnuplot
    assert '# X Y \n2.0\t7.0\n' in gnuplot


def test_gnuplot_plot_line():
    unlabeled = [dataset([1.0], [3.0]), dataset([2.0], [4.0])]
    plot = [l for l in makeGnuplotPlot(unlabeled).splitlines() if l.startswith('plot')][0]
    assert plot == 'plot  "-" using 1:2 notitle , "-" using 1:2 notitle '
    # An empty dataset has no error bars and must not break the next one.
    gnuplot = makeGnuplotPlot([dataset([], []), dataset([1.0], [3.0])])
    assert '1.0\t3.0' in gnuplot


def test_grace_empty_dataset():
    # An empty dataset must not make the next one reuse its set number.
    grace = makeGracePlot([dataset([], []), dataset([1.0], [3.0])])
    targets = [l for l in grace.splitlines() if l.startswith('@target')]
    assert targets == ['@target G0.S0', '@target G0.S1']
    assert data_lines(grace) == ['1.0\t3.0']


def test_mismatched_lengths():
    # x and y of different lengths must raise instead of failing midway or dropping data.
    for x, y in [([1.0, 2.0, 3.0], [1.0, 2.0]), ([1.0, 2.0], [1.0, 2.0, 3.0]), ([1.0], [])]:
        for make_plot in (convertToText, makeGracePlot, makeGnuplotPlot):
            try:
                make_plot([dataset(x, y)])
            except ValueError:
                pass
            else:
                raise AssertionError('%s accepted x and y of lengths %d and %d' % (make_plot.__name__, len(x), len(y)))


if __name__ == '__main__':
    test_multi_point_errors()
    test_grace_x_errors()
    test_gnuplot_errors_per_dataset()
    test_gnuplot_plot_line()
    test_grace_empty_dataset()
    test_mismatched_lengths()

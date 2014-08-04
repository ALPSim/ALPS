
import numpy as np
from pyngsaccumulator_c import count_accumulator, mean_accumulator, error_accumulator, binning_analysis_accumulator, max_num_binning_accumulator

def checkResult(r):
	print r
	r *= 2
	print r
	print -r
	print r + 2
	print 2 + r
	print r + r
	print r - r
	print r * r
	print r / r
	print np.sin(r)
	print r.sin()
	print r.cos()
	print r.tan()
	print r.sinh()
	print r.cosh()
	print r.tanh()
	# print r.asin()
	# print r.acos()
	# print r.atan()
	print r.abs()
	print r.sqrt()
	print r.log()
	# print r.sq()
	# print r.cb()
	# print r.cbrt()

a = count_accumulator()
print a
a(1)
print a
checkResult(a.result())

b = mean_accumulator()
b(10)
print b
checkResult(b.result())

c = error_accumulator()
c(8)
c(12)
print c
checkResult(c.result())

d = binning_analysis_accumulator()
for i in range(1000):
	d(float(i))
print d
checkResult(d.result())

e = max_num_binning_accumulator()
for i in range(1000):
	e(float(i))
print e
r = e.result()
print r
# r *= 2
print r
print -r
# print r + 2
# print 2 + r
# print r + r
# print r - r
# print r * r
# print r / r

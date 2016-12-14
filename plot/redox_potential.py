from numpy import *
from scipy.optimize import *
import matplotlib.pyplot as plt
import sys

'''This program is to stimulate the function of Activity-Potential, and
calculate the pivot potential.
Data format:
X   Y   stdev.p(y)'''

guess = (1, 1, -1)
t = arange(-.4, -0.2, 1e-3)
half = 0.5
# 50% activity
sample = 'Example'


def logistics(x, a, b, c):
    return a/(1+b*exp(-1*c*x))


def fit(x):
    return logistics(x, a, b, c) - half

if len(sys.argv) > 1:
    with open(sys.argv[1], 'r') as data:
        raw = data.read()
        raw_1 = raw.split(sep='\n')
        raw_1.pop(-1)
        raw_1.pop(0)
        raw_2 = [i.split(sep='\t') for i in raw_1]
        sample = raw_2[0][0]
        # first line, first row
        xdata = [float(i[1]) for i in raw_2]
        ydata = [float(i[2]) for i in raw_2]
        yerror = [float(i[3]) for i in raw_2]
else:
    raise ValueError()

popt, pcov = curve_fit(logistics, xdata, ydata, p0=guess, maxfev=500000)
a, b, c = popt
print('function: y={0}/(1+{1}*exp(-1*{2}*x))'.format(a, b, c))
print('Standard deviation error: A {0:.3f} B {1:g} C {2:.3f}'.format(
    *sqrt(diag(pcov))))

pivot = fsolve(fit, -0.3)
plt.ylabel('Activity (%)')
plt.xlabel('Potential/V')
plt.plot(t, logistics(t, a, b, c), color='k',
         label='pivot={0:.3f}mV'.format(pivot[0] * 1000))
plt.legend(loc='upper right')
plt.errorbar(xdata, ydata, yerr=yerror, fmt='o', color='k', ecolor='k')
plt.show()

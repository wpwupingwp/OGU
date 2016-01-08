from numpy import *
from scipy.optimize import *
import matplotlib.pyplot as plt
from matplotlib import rcParams
import sys

'''Data format:
Sample  x   y   stdev.p(y)'''

guess = (1, 1e10, -5)
t = arange(-7.5,-3,1e-3)
ic = 0.5
#result from SPSS 21, R square equals 0.998
a = 1.023
b = 7.172e10
c = -4.6
sample = 'Example'
def logistics(x, a, b, c):
    return a/(1+b*exp(-1*c*x))

def fit(x):
    return logistics(x, a, b, c) - ic

if len(sys.argv) >1:
    with open(sys.argv[1], 'r') as data:
        raw = data.read()
        raw_1 = raw.split(sep='\n')
        raw_1.pop(-1)
        raw_1.pop(0)
        raw_2 = [i.split(sep='\t') for i in raw_1]
        sample = raw_2[0][0]
        xdata = [float(i[1]) for i in raw_2]
        ydata = [float(i[2]) for i in raw_2]
        yerror = [float(i[3]) for i in raw_2]
else:
    xdata = [-7,-6.3,-6,-5.7,-5.3,-5.1,-5,-4.7,-4.4]
    ydata = [1.04,0.99,0.93,0.82,0.33,0.20,0.13,0.03,0.02]
    yerror = [0.05,0.09,0.07,0.09,0.02,0.02,0.03,0.01,0.02]

popt, pcov = curve_fit(logistics, xdata, ydata, p0=guess, maxfev=50000)
a, b, c = popt
print('function: y={0}/(1+{1}*exp(-1*{2}*x))'.format(a, b, c))
print('Standard deviation error: A {0:.3f} B {1:g} C {2:.3f}'.format(*sqrt(diag(pcov))))
ic50_log10 = fsolve(fit, -5)
ic50 = 10**(6+ic50_log10)
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Tahoma']
plt.ylabel('Activity (% of control)')
plt.xlabel('Inhibitor (Log μM)')
plt.plot(t, logistics(t, a, b, c), label='IC50={0:.3f} μM'.format(ic50[0]))
plt.title(sample)
plt.legend(loc='upper right')
plt.errorbar(xdata, ydata, yerr=yerror, fmt='o', color='c', ecolor='r')
plt.show()
print(ic50_log10)
print('IC50 of this component is {0:.3f} μM'.format(ic50[0]))

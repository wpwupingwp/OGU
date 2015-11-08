from numpy import *
from scipy.optimize import *
import matplotlib.pyplot as plt

xdata = [-7,-6.3,-6,-5.7,-5.3,-5.1,-5,-4.7,-4.4]
ydata = [1.04,0.99,0.93,0.82,0.33,0.20,0.13,0.03,0.02]
t = arange(-8,-3,1e-3)
ic = 0.5
#result from SPSS 21, R square equals 0.998
spss_a = 1.023
spss_b = 7.172e10
spss_c = -4.6

def logistics(x, a, b, c):
    return a/(1+b*exp(-c*x))

def fit(x):
    return logistics(x, spss_a, spss_b, spss_c) - ic

plt.plot(t, logistics(t, spss_a, spss_b, spss_c), xdata, ydata, 'r^')
plt.show()
ic50 = 10**(6+fsolve(fit, -5))
print('IC50 of this component is {0} μM'.format(ic50[0]))
a, b = curve_fit(logistics, xdata, ydata, p0=(1, 1e10, -5))
print(a, b)

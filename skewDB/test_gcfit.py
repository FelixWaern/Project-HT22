#from IPython.display import set_matplotlib_formats
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()

import matplotlib
import math
import numpy as np
#import seaborn as sb
import matplotlib.pyplot as plt
plt.rcParams['figure.figsize'] = [9.5, 7]
import pandas
from random import shuffle
#from scipy.stats import poisson
#import statsmodels.api as sm
#import statsmodels.formula.api as smf

# set this to where you placed https://berthub.eu/antonie/gcskewdb.csv.gz and unpacked https://berthub.eu/antonie/gcfits.tar.bz2
#prefix="C:\Ashwini\Applied bioinformatics\gcfits"
m=pandas.read_csv("https://skewdb.org/view/gcskewdb.csv")
#print(m)
chromoname="NC_001263.1" 
fitted = pandas.read_csv(r"C:\Ashwini\Applied bioinformatics\gcfits\NC_001263.1_fit.csv")

us=m[m.name==chromoname]
plt.figure()
plt.plot(fitted.pos, fitted.gc2skew, label="Cumulative GC skew")
plt.plot(fitted.pos, fitted.predgc2skew, label="Fitted GC skew")
plt.show()

"""leshift=us["shift"].item()
print(leshift)
if leshift > 0:
    plt.axvline(leshift, ls=':', color='black')
else:
    plt.axvline(fitted.pos.max() + leshift, ls=':', color='black')

plt.xlabel("Locus")
plt.ylabel("Skew")
plt.legend()
plt.title(chromoname + " alpha1 " + str(round(us["alpha1gc"].item(),3)) + " alpha2 " + str(round(us["alpha2gc"].item(),3)) + 
          " div " + str(round(us["div"].item(),3)))
plt.grid()
plt.show()"""

"""# how is the split between leading/lagging strand distributed?
plt.figure()
plt.hist(m["div"], bins=50, density=True)
plt.grid()
plt.title("Division between leading/lagging strand")
plt.show()"""
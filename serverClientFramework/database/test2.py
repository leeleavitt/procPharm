# pip install rpy2
# Path to R must be corrected.
import os
os.environ["R_HOME"] = r"C:\Program Files\R\R-4.0.3"
os.environ["PATH"]   = r"C:\Program Files\R\R-4.0.3\bin\x64" + ";" + os.environ["PATH"]
import matplotlib.pyplot as plt
import rpy2.robjects as robjects
robjects.r['load'](r"C:\Users\rishi\Documents\gitrepos\procPharm\extras\RD.200309.30.m.m3.p1.Rdata")
a = robjects.r['RD.200309.30.m.m3.p1']
plt.plot(a[0][0:1], a[0][0:5])
plt.show()
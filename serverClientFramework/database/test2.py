from imaging_results import Imaging

img = Imaging(r"C:\Users\Arikill\Documents\gitrepos\procPharm\extras\RD.200309.30.m.m3.p1.Rdata", 200309, 30, "m", "m3", "p1", "drg", "")
img.organizeCellData()
img.writeToDB()

# pip install rpy2
# Path to R must be corrected.
# import os
# os.environ["R_HOME"] = r"C:\Program Files\R\R-4.0.3"
# os.environ["PATH"]   = r"C:\Program Files\R\R-4.0.3\bin\x64" + ";" + os.environ["PATH"]
# import numpy as np
# import matplotlib.pyplot as plt
# import rpy2.robjects as robjects
# data_object = robjects.r['load'](r"C:\Users\rishi\Documents\gitrepos\procPharm\extras\RD.200309.30.m.m3.p1.Rdata")[0]
# data = robjects.r[data_object]
# print(data[2][12])
# a = robjects.r['RD.200309.30.m.m3.p1']
# names = a[2][0]
# areas = a[2][1]
# position_x = a[2][2]
# position_y = a[2][3]
# circularity = a[2][4]
# perimeter = a[2][5]
# mean_bf_start = a[2][6]
# mean_bf_end = a[2][7]
# mean_gfp_start = a[2][8]
# mean_gfp_end = a[2][9]
# mean_dapi = a[2][10]
# mean_cy5_end = a[2][11]

# time = np.reshape(np.asarray(a[0][0], dtype=np.float32), [1, len(a[0][0])])
# print(time)
# cells = [{}]
# values = np.zeros([len(a[0][1:]), len(a[0][0])], dtype=np.float32)
# for i in range(len(a[0][1:])):
#     values[i, :] = a[0][i]
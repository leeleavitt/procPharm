import os
os.environ["R_HOME"] = r"C:\Program Files\R\R-4.0.3"
os.environ["PATH"]   = r"C:\Program Files\R\R-4.0.3\bin\x64" + ";" + os.environ["PATH"]
import rpy2.robjects as robjects
import numpy as np
from pyArango.connection import *

class Imaging:
    def __init__(self, filename):
        data_object = robjects.r['load'](filename)[0]
        data = robjects.r[data_object]
        self.names = data[2][0]
        self.areas = np.float32(data[2][1])
        self.position_x = np.float32(data[2][2])
        self.position_y = np.float32(data[2][3])
        self.circularity = np.float32(data[2][4])
        self.perimeter = np.float32(data[2][5])
        self.bf_end = np.float32(data[2][6])
        self.bf_start = np.float32(data[2][7])
        self.gfp_end = np.float32(data[2][8])
        self.gfp_start = np.float32(data[2][9])
        self.dapi = np.float32(data[2][10])
        self.cy5_end = np.float32(data[2][11])
        self.time = np.reshape(np.asarray(data[0][0], dtype=np.float32), [1, len(data[0][0])])
        self.values = np.reshape(np.asarray(data[0][1:], dtype=np.float32), [len(data[0][1:]), len(data[0][0])])
        pass

    def organizeCellData(self):
        self.cells = [{}]*len(self.names)
        for index, name in enumerate(self.names):
            self.cells[index] = {
                "name": name, 
                "area": self.areas[index],
                "x": self.position_x[index],
                "y": self.position_y[index],
                "bright-field": {
                    "start": self.bf_start[index],
                    "end": self.bf_end[index]
                },
                "gfp": {
                    "start": self.gfp_start[index],
                    "end": self.gfp_end[index]
                },
                "dapi": self.dapi[index],
                "cy5": {
                    "end": self.cy5_end[index]
                },
                "times": self.time,
                "ratio-values": self.values[index, :]
            }
        pass

    def writeToDB(self):
        conn = Connection(username="root", password="roselab2")
        db = conn["test"]
        collection = db.createCollection(name="cells")
        for cell in self.cells:
            doc = collection.createDocument()
            for key in cell:
                doc[key] = cell[key]
            doc.save()
        pass
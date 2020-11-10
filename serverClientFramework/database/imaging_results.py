import os
os.environ["R_HOME"] = r"C:\Program Files\R\R-4.0.3"
os.environ["PATH"]   = r"C:\Program Files\R\R-4.0.3\bin\x64" + ";" + os.environ["PATH"]
import rpy2.robjects as robjects
import numpy as np
from pyArango.connection import *

class Imaging:
    def __init__(self, filename, date, age, gender, microscope, slide, tissue="drg", info=""):
        self.date = np.int(date)
        self.age = float(age)
        self.gender = str(gender)
        self.microscope = str(microscope)
        self.slide = str(slide)
        self.tissue = str(tissue)
        self.info = str(info)
        data_object = robjects.r['load'](filename)[0]
        data = robjects.r[data_object]
        self.names = data[2][0]
        self.areas = data[2][1]
        self.position_x = data[2][2]
        self.position_y = data[2][3]
        self.circularity = data[2][4]
        self.perimeter = data[2][5]
        self.bf_end = data[2][6]
        self.bf_start = data[2][7]
        self.gfp_end = data[2][8]
        self.gfp_start = data[2][9]
        self.dapi = data[2][10]
        self.cy5_end = data[2][11]
        self.time = np.asarray(data[0][0], dtype=np.float32)
        self.values = np.asarray(data[0][1:], dtype=np.float32)
        pass

    def organizeCellData(self):
        self.cells = [{}]*len(self.names)
        for index, name in enumerate(self.names):
            self.cells[index] = {
                "date": self.date,
                "age": self.age,
                "gender": self.gender,
                "tissue": self.tissue,
                "info": self.info,
                "name": name, 
                "area": float(self.areas[index]),
                "x": float(self.position_x[index]),
                "y": float(self.position_y[index]),
                "bright-field": {
                    "start": float(self.bf_start[index]),
                    "end": float(self.bf_end[index])
                },
                "gfp": {
                    "start": float(self.gfp_start[index]),
                    "end": float(self.gfp_end[index])
                },
                "dapi": float(self.dapi[index]),
                "cy5": {
                    "end": float(self.cy5_end[index])
                },
                "times": self.time.tolist(),
                "ratio-values": self.values[index, :].tolist()
            }
        pass

    def writeToDB(self):
        conn = Connection(username="root", password="roselab2")
        try:
            db = conn.createDatabase(name="test")
        except:
            print("DB already exists")
        finally:
            db = conn["test"]
        try:
            collection = db.createCollection(name="cells")
        except:
            print("collection already exists")
        finally:
            collection = db["cells"]
        for cell in self.cells:
            doc = collection.createDocument()
            for key in cell:
                doc[key] = cell[key]
            doc.save()
        pass
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 21 16:08:56 2018

@author: DarthMaul
"""

#This is how to import thetime info into a format that r can use

import os
from nd2reader import ND2Reader
import glob

main_dir=os.getcwd()
#main_dir="Z:/Matias/181219.22.m.m3.p1"
os.chdir(main_dir)

video_names=glob.glob('video*')

print(30 * '-')
print("Video File Selection")
print(30 * '-')

for i in range(len(video_names)):
    print(i, video_names[i])

choice=input("Enter your choice: ")
choice=int(choice)

import numpy

video_data=ND2Reader(video_names[choice])

arr_val2=numpy.arange(1,len(video_data.timesteps)+1)

foo=numpy.column_stack((arr_val2, video_data.timesteps))

arr_to_export=numpy.asarray([ [arr_val2], [video_data.timesteps] ])

numpy.savetxt("time.info.txt",foo,delimiter="\t",header='index\tTime[s]')


#w.pack()
#selection=w.mainloop()


#video_data=ND2reader(

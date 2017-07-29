"""
Author: Chukwubuikem Ume-Ugwa
Email: chubiyke@gmail.com

MIT License

Copyright (c) 2017 CleverChuk

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""


import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D    
from dataparser import *

from multiprocessing import Pool, Manager
import os
import time

manager = Manager()
heightD = manager.dict() # holds values for minimum height of each particle

TSIZE = 8 # data type size in bytes
N_OFFSETS = 44 # number of data
FCOLOR = genColor(N_OFFSETS, manager)

# Dimension of the simulation bed
xsize = 78
ysize = 112
zsize = 104

hOut = "HeightData"
def startChild(fname):
	# DISALLOWED IN PYTHON
	iam.fn = fname
	dictn = iam.manager.dict()
	mylist = iam.manager.list()

	pool = Pool()
	# passing offset multiplier to the producer task
	pool.map(iam.producer, [i for i in range(1 , iam.N_OFFSETS)], 1)

	# Feeds task from producers into the list
	for i, j in self.dictn.items():
		mylist.append(j[0])


	# single process to handle plotting
	proc = Process(target=iam.consumer, args=(mylist, ))
	proc.start()
	proc.join()



def multijob(fname):
	"""
	Handles reading and plotting of data in file with name fname
	"""

	print("Starting multijob from process: %d" % os.getpid())
	
	fig = plt.figure()
	axis = Axes3D(fig)
	heightL = manager.list()
	
	axis = Axes3D(fig)
	axis.set_xlim([0,ysize])
	axis.set_ylim([0,ysize])
	axis.set_zlim([0,ysize])
	axis.view_init(elev = 40, azim = 50)

	coords = manager.list()
	rho = readsingle(fname)

	for i in range(1, N_OFFSETS):
		
		eta_s = readsingle(fname, i * TSIZE)
#		eta_s = process(rho, filter_eta(eta_s))
		coords.append(getcoords(eta_s, xsize, ysize, zsize))
		heightL.append(max(coords[-1][-2]) - min(coords[-1][-2]))

	
	writtable(hOut,str(heightL).strip('[]'))
		
	plot(coords, fig, axis, count = "ALL", fcolor = FCOLOR, multijob = (True,fname))		
	print("Finished multijob from process: %d" % os.getpid())




if __name__ == "__main__":
	print("Starting mutiple jobs in a process task")
	import timeit, sys
	start_time = timeit.default_timer()
	if(os.path.exists(hOut)):
		os.remove(hOut)

	
	pool = Pool()	
	files = list()

	MAXCOUNT = 4
	STEP = 2
	START = 0
	FNAME = "fullT{0}.dat"
	

	## file with filesname to work on
	for i in range(START, MAXCOUNT, STEP):		
		files.append(FNAME.format(i))

	
	pool.map(multijob, files, 1)
	elapsed = timeit.default_timer() - start_time
	
	print("total time %d seconds" % elapsed)
	print("Finished multiple job in a process task")














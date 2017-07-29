import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib as mpl
import matplotlib.animation as animation
from pylab import *
import pickle as pk
from multiprocessing import Process as pc
        
from dataparser import *
from multiprocessing import Process, Queue, cpu_count, Pool, Manager
from threading import Thread

xsize = 78
ysize = 112
zsize = 104

fig = plt.figure()
ax = Axes3D(fig)
ax.set_xlim([0,ysize])
ax.set_ylim([0,ysize])
ax.set_zlim([0,ysize])
ax.view_init(elev=40,azim=50)

"""
        Parallel Implemettation
"""
# constants
fn = None 
TSIZE = 8 # data type size in bytes
N_OFFSETS = 44 # number of data


# Creates process safe objects
manager = Manager()
ncores = cpu_count()
mylist = manager.list()
dictn = manager.dict()

#MyManager.register("MyFigure",MyFigure)
#mManager = MyManager()
#mManager.start()
#figure = mManager.MyFigure(ysize)


rho = None
import os, random

def producer(start):
        """
	handles making coordinate for each particle
	"""
	print ("starting producer job in process: %d" % os.getpid())

	eta_s = readsingle(fn, start * TSIZE)
	eta_s = process(rho,filter_eta( eta_s))

	coords = getcoords(eta_s, xsize, ysize, zsize)
	key = os.getpid() + random.randint(1,1000)
	if key in dictn.keys():
		key += 1

	dictn[key] = coords
        print ("finishing producer job in process: %d" % os.getpid())

def consumer(q, proxy=None ):
	"""
	handles plotting the data
	"""
        print ("starting consumer job")
	if not proxy == None:
		plot(q, proxy[0], proxy[1], start = 0, count = "ALL")
	else:
		plot(q, fig, ax, start = 0, count="ALL")
#		figure.plot(q)
        print ("finishing consumer job")

	


if __name__ == "__main__":
	print("Starting Parallel Tasks")

	import timeit
	start_time = timeit.default_timer()

	for num in range(0,4,2):
		fn = "fullT%d.dat" % num
		print("Beginning job on file: %s" %fn)
		rho = readsingle(fn)
		# Gets producer task workers. Defualts to number of CPU cores
		pool = Pool()

		# passing offset multiplier to the producer task
		pool.map(producer, [i for i in range(1 , N_OFFSETS)], 1)

		# Feeds task from producers into the list
		for i, j in dictn.items():
			mylist.append(j)

		# Creates a consumer process to handle plotting
		# single process to handle plotting
		proc = Process(target=consumer, args=(mylist, ))
		proc.start()
		proc.join()

		print("Finished job on file: %s" %fn)



	elapsed = timeit.default_timer() - start_time
	
	print("total time %d seconds" % elapsed)
	print("Parallel Tasks Completed")
















"""
Authors:
	Chukwubuikem Ume-Ugwa
	Georgina Obehi
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
import numpy as np
from multiprocessing.managers import BaseManager
import random
from multiprocessing import Manager
 
def getnumzeros(vector):
	"""
	Counts the number of zeroes in vector
	"""
        count = 0
        for i in vector:
                if str(i) == "0.0" or str(i) == "0":
                        count += 1
        return count

def readRow(filename):
    """
    Return the first line in file
    """
    with open(filename,'r') as file:
        row  = file.readline()
        return row


def getProduct(row):
    """
    Returns an array of the products of first element and the rest
    of the elements in row
    """
    temp = row.split('\t')
    data = []

    for i in temp[1:]:
        data.append(float(temp[0])*float(i))

    return data



def parseData(filename):
    """
    Return and m by n matrix of data from a text file and number of lines read
    """
    parsed = []
    numLines = 0
    with open(filename,'r') as file:
        
        while True:        
            row  = file.readline()
            
            if not row:
                break;
            parsed.append(getProduct(row));
            numLines += 1

        file.close()
    return (parsed,numLines)



def unpack(matrix):
    """
    unpacks the matrix so that each column represent a single particle.
    returns the new matrix
    """
    numP = len(matrix[0])
    unpacked = list()

    for i in range(numP):
        unpacked.append(list())

    for j in matrix:
        for k in range(numP):
            unpacked[k].append(j[k])

    return unpacked
            
def genColor(n = 1, manager = None):
        """
        Generates n hexadecimal color codes.
        """
	temp = set()
	random.seed(random.randint(100000, 5000000000))
        r = lambda: random.randint(0,16777216)

        while True:
		temp.add(('#%06x' % r()))
		if len(temp) == n:
			break
        
	if manager == None:
		return list(temp)
	else:
	         return manager.list(temp)
    


"""
 Data class to represent C structures
"""

"""
struct { double rho; double eta;};
"""
from ctypes import *
import struct
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D 
class Particle(Structure):
        """
        Represents the C structure
        """
        _fields_ = [
                ('rho', c_double),
                ('eta', c_double),
                ]
        def __new__(self, cstruct = None):
                return self.from_buffer_copy(cstruct)                   

def decode(filename):
        """
        Reads a single C structure from file
        """
        data = Particle()          
        with open(filename, 'rb') as file:
                file.readinto(data)
        return data

def getfmt(n,dtype ="d"):
        """
	Generates the binary format n dtype
	defualts to double
        @param n: number of particles
        @param dtype: data type
        returns the format string used to decode the binary data
        """
        fmt = ""
        fmt += dtype
        for i in range(n):
                fmt += dtype
        return fmt

def readBinary(fname,fmt,bsize):
        """     
	Reads data from a data file with name fname using the binary
	format fmt. Reads exactly data of size bsize until EOF

        @param fname: filename
        @param fmt: format string represent how the data is in memory
        @param bsize: size of the data in bytes to read
        
        returns a matrix of data
        """
        def process(data):
                num = len(data)
                pdata = list()
                for i in range(1, num):
                        pdata.append(float(data[0]) * float(data[i]))
                return pdata

        import struct as st
        data = list()
        with open(fname,'rb') as file:
                while True:             
                        rbuffer = file.read(bsize)
                        if not rbuffer:
                                break
                        if len(rbuffer) < bsize:
                                continue
                        data.append(process(st.unpack(fmt,rbuffer)))
        return data
n  = 0
def getcoords(data, xsize, ysize, zsize):
        """
	Generates the x,y,z coordinate for data
        @param data: data to plot
        @param xsize:
        @param ysize:
        @param zsize:
        """
     	global n
	coords = []
        if type(data[0]) == list: #checks if data needs to be unpacked
				# raw data from readBinary has to be unpacked
                pdata = unpack(data)                                  

                for V in pdata:
                        n += 1
                        size = len(V)
                        countz = getnumzeros(V)
                        countnz = size - countz

                                
                        Xd = [0 for x in range(size)]
                        Yd = [0 for x in range(size)]
                        Zd = [0 for x in range(size)]
                        Vd = [0 for x in range(size)]
                        
                        print("Particle number %d" % n); Vd = V
                        
                        ad = np.zeros(size)
                        bd = np.zeros(size)
                        cd = np.zeros(size)
                        vd = np.zeros(size)

                        counti = 0
                        for z in range(zsize):
                                for y in range(ysize):
                                        for x in range(xsize):
                                                Zd[counti] = x
                                                Yd[counti] = z
                                                Xd[counti] = y
                                                counti = counti+1
        # Adjusted the limit since Roh * Eta is less than 0.1
                        lim = 0.1
                        countl = 0
                        for i in range(0,size):
                                if Vd[i] > lim:
                                        ad[i] = Xd[i]
                                        bd[i] = Yd[i]
                                        cd[i] = Zd[i]
                                        vd[i] = Vd[i]
                                        countl = countl+1
                                else:
                                        ad[i] = 1000.0
                                        bd[i] = 1000.0
                                        cd[i] = 1000.0
                                        vd[i] = 1000.0

                        countth = np.count_nonzero(vd==1000.0)

                        val = np.array([1000.0])

                        x = np.setdiff1d(ad,val,assume_unique = True)
                        y = np.setdiff1d(bd,val,assume_unique = True)
                        z = np.setdiff1d(cd,val,assume_unique = True)
                        v = np.setdiff1d(vd,val,assume_unique = True)
                        
                        coords.append([x,y,z,v])
        else:
                V = data                 
                n += 1
		size = len(V)
                countz = getnumzeros(V)
                countnz = size - countz

                        
                Xd = [0 for x in range(size)]
                Yd = [0 for x in range(size)]
                Zd = [0 for x in range(size)]
                Vd = [0 for x in range(size)]
                
		print("Particle number %d" % n); Vd = V
                
                ad = np.zeros(size)
                bd = np.zeros(size)
                cd = np.zeros(size)
                vd = np.zeros(size)

                counti = 0
                for z in range(zsize):
                        for y in range(ysize):
                                for x in range(xsize):
                                        Zd[counti] = x
                                        Yd[counti] = z
                                        Xd[counti] = y
                                        counti = counti+1
# Adjusted the limit since Roh * Eta is less than 0.1
                lim = 0.1
                countl = 0
                for i in range(0,size):
                        if Vd[i] > lim:
                                ad[i] = Xd[i]
                                bd[i] = Yd[i]
                                cd[i] = Zd[i]
                                vd[i] = Vd[i]
                                countl = countl+1
                        else:
                                ad[i] = 1000.0
                                bd[i] = 1000.0
                                cd[i] = 1000.0
                                vd[i] = 1000.0

                countth = np.count_nonzero(vd==1000.0)

                val = np.array([1000.0])

                x = np.setdiff1d(ad,val,assume_unique = True)
                y = np.setdiff1d(bd,val,assume_unique = True)
                z = np.setdiff1d(cd,val,assume_unique = True)
                v = np.setdiff1d(vd,val,assume_unique = True)
                
                coords = [x,y,z,v]
                
        return coords
    


def plot(coords, fig, ax, start = 0, count = 10,fcolor = None,  multijob = (False,"filename")):
	"""
	plots the data in coords on ax starting from start and ending at count
	"""
	pcount = len(coords) - 1
	
	if fcolor == None:	
		fcolor = genColor(len(coords))

	if len(coords) == 4:	
		    x,y,z,v = coords;# print("sizes: x =%d y =%d z = %d" % (len(x[0]),len(y[0]), len(z[0])))
		    ax.scatter(x, y, z, s = 100*v, marker="o", facecolor = fcolor[pcount], linewidths = 1.2, edgecolor = "black")

	elif type(count) is str and count.upper() == "ALL": 
			for point in coords:
				    x,y,z,v = point;# print("sizes: x =%d y =%d z = %d" % (len(x[0]),len(y[0]), len(z[0])))
				    ax.scatter(x, y, z, s = 100*v, marker=".", facecolor = fcolor[pcount], linewidths = 1.2, edgecolor = "black")
                   		    pcount -= 1

	else:
			for i in range(start,count):
				    i = 0;
				    x,y,z,v = coords[i]
				    ax.scatter(x, y, z, s = 100*v, marker=".", facecolor = fcolor[pcount], linewidths = 1.2, edgecolor = "black")
				    pcount -= 1
	if multijob[0]:
		fig.savefig(multijob[1]+".png")
	else:
		num = random.randint(0,10)
		fig.savefig("figure%d.png" % num)
	print("Done... from plot")
        

def readsingle(file, offset = 0, ndata = 44, dtype='d'):
	"""
	reads ndata  of type dtype from file using offset as the byte offset
	to control read location.
	"""

        dtypes = {'c':1, 'i':4, 'd':8}
        tsize = dtypes[dtype]        
        data = list()
        
        if not hasattr(file,'read'): #checks if file has been opened for reading
                with open(file,"rb") as file:
                        count = 0
                        while True:                        
                                if count == 0:
                                        file.seek(offset)
                                        buffer = file.read(tsize)                                
                                else:
                                        file.seek(count * ndata * tsize + offset)
                                        buffer = file.read(tsize)

                                count += 1                        
                                if not buffer:
                                        break
                                
                                data.append(struct.unpack(dtype, buffer)[0]) # first element of the tuple
        else:
                count = 0
                while True:                        
                        if count == 0:
                                file.seek(offset)
                                buffer = file.read(tsize)                                
                        else:
                                file.seek(count * ndata * tsize + offset)
                                buffer = file.read(tsize)

                        count += 1                        
                        if not buffer:
                                break
                        
                        data.append(struct.unpack(dtype, buffer)[0]) # first element of the tuple
                file.close()
        
                
        print("Done... from readsingle")
        return data
                        
def process(rho, eta):
	"""
	returns the element by element  product of two vectors 
	"""
        for i in xrange(len(eta)):
                eta[i] = eta[i] * rho[i]
        return eta



def filter_eta(eta, K = 0.0041):
	for i in range(len(eta)):
		if(eta[i] >= K):
			eta[i] = 1.0
		else:
			eta[i] = 0.0
	return eta
def writtable(filename, writtable):
	"""
	writes writtable data to file with name filename
	"""
	with open(filename,'a') as file:
		file.write(writtable+"\n")
	        
def readable(filename):
	"""
	reads data from file with name filename
	"""
	buffer = []
	with open(filename,"r") as file:
		while True:
			data = file.readline()
			
			if not data:
				break
			buffer.append(data.strip().split(','))
	return buffer

def plotLine(t, y, numParticles = 43):
	numPoints_per_Particle = len(t)

	for i in xrange(numParticles):	
		temp = []

		for j in xrange(numPoints_per_Particle):
			temp.append(float(y[j][i]))

		plt.plot(t,temp, linewidth = 1.5, color = genColor(), marker = '.')
	plt.show()

	

class MyManager(BaseManager):
	"""
	Blueprint for process object manager
	"""
	pass

class MyFigure:
	"""
	Blueprint for objects managed between processes
	"""
	figure = plt.figure()
	axis = Axes3D(figure)

	def __init__(self, ysize):
		self.axis.set_xlim([0, ysize])
		self.axis.set_ylim([0, ysize])
		self.axis.set_zlim([0, ysize])
		self.axis.view_init(elev = 40, azim = 50)

	def plot(self,coords, start = 0, count = "ALL"):
		num = 0

		if len(coords) == 4:	
			    x,y,z,v = coords;# print("sizes: x =%d y =%d z = %d" % (len(x[0]),len(y[0]), len(z[0])))
			    self.axis.scatter(x, y, z, s = 100*v, marker="o", facecolor = genColor(), linewidths = 1.2, edgecolor = genColor())
			    num += 1

		elif type(count) is str and count.upper() == "ALL": 
				for point in coords:
					    x,y,z,v = point;# print("sizes: x =%d y =%d z = %d" % (len(x[0]),len(y[0]), len(z[0])))
					    self.axis.scatter(x, y, z, s = 100*v, marker="o", facecolor = genColor(), linewidths = 1.2, edgecolor = genColor())
					    num += 1
		else:
				for i in range(start,count):
					    i = 0;
					    x,y,z,v = coords[i]
					    self.axis.scatter(x, y, z, s = 100*v, marker="o", facecolor = genColor(), linewidths = 1.2, edgecolor = genColor())
					    num += 1

		self.figure.savefig("figure%d.png" % num)
		print("Done... from plot in MyFigure class")

"""
MAX = 30
with open("fullT10.dat","rb") as file:
        count = 0
        while True:
                data = file.read(1080)
                if not data:
                        break
                data  = struct.unpack(getfmt(134),data)
                print(data)
                if count > MAX:
                        break
                count += 1
"""
"""
    Example:

        filename = "full0.txt"
        data = parseData(filename)
        parsed = np.array(data[0])
        print(parsed)

        filename = "full0.txt"
        data = parseData(filename)

        pdata = unpack(data[0])

        print(pdata.count(.998))
        parsed = np.array(unpack(data[0]))
        print(parsed)
"""














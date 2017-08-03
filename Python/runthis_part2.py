import subprocess, sys
from itertools import combinations

# the firs element is the binary to run. Assumed to be in the same dir as this script.
# the second element  is the first arg to passed to the binary. In this case, it is the simulation time'
# the third element is the second arg  and it is the time to start contraining the lowerbound;
# this time should always be at least half the simulation time but always less than the simulation time
cmd = ["./allcombinations","15","10"]

# array of indices used to generate combinations to be run
params = [0,1,2,3,4,5,6,7]

# sample size
NUM = len(params)

def buildCmd(args):
	"""
	builds the command used to run the algorithm
	to ensure that each command is performed
	
	args: at tuple of combinations built from params

	command: the new command generated using cmd and args
	"""
	global cmd
	command = cmd[:]
	for i in args:
		command.append(str(i))

	return command

if __name__ == "__main__":
	for i in xrange(1, NUM + 1):
		obj = combinations(params, i)

		if( i < 4):
			cmd[1] = "2"
			cmd[2] = "1"

		cmd.append(str(i)); 
		# Each combination is run in a different process
		while True:
			try:
				arg = obj.next()
				new_command = buildCmd(arg)
				subprocess.call(new_command)
				#print("Finished running this command: %s" % str(new_command))
				
			except Exception as e:
				print e
				break			
		
		cmd.pop()

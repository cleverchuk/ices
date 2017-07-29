import subprocess, os
from multiprocessing import pool

def job(num):
	cmd = ["./a.out","15","0",str(num)]
	subprocess.call(cmd)

if __name__ == "__main__":
	pool = pool.Pool()
	pool.map(job,[ i for i in range(1,8)])
	cmd = "python email_sender.py"
	os.system(cmd)

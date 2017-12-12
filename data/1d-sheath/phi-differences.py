import numpy as np
from scipy import interpolate

exactdata = np.loadtxt("1d-test-files/unperturbed.dat", usecols=(0,1))
exact = interpolate.interp1d(exactdata[:,0], exactdata[:,1])

simulation = np.loadtxt("raw/width10-50x500.dat", usecols=(1,3))
differences = []
for i in range(5,500):
  differences.append([simulation[i,0], simulation[i,1]-exact(simulation[i,0]-0.09)])
np.savetxt("differences/width10-50x500diff.dat",differences)

simulation = np.loadtxt("raw/width20-50x1000.dat", usecols=(1,3))
differences = []
for i in range(5,1000):
  differences.append([simulation[i,0], simulation[i,1]-exact(simulation[i,0]-0.09)])
np.savetxt("differences/width20-50x1000diff.dat",differences)

simulation = np.loadtxt("raw/width30-50x1500.dat", usecols=(1,3))
differences = []
for i in range(5,1500):
  differences.append([simulation[i,0], simulation[i,1]-exact(simulation[i,0]-0.09)])
np.savetxt("differences/width30-50x1500diff.dat",differences)

simulation = np.loadtxt("raw/width40-50x2000.dat", usecols=(1,3))
differences = []
for i in range(5,2000):
  differences.append([simulation[i,0], simulation[i,1]-exact(simulation[i,0]-0.09)])
np.savetxt("differences/width40-50x2000diff.dat",differences)

simulation = np.loadtxt("raw/width20-12x240.dat", usecols=(1,3))
differences = []
for i in range(1,240):
  differences.append([simulation[i,0], simulation[i,1]-exact(simulation[i,0]-0.0416667)])
np.savetxt("differences/width20-12x240diff.dat",differences)

simulation = np.loadtxt("raw/width20-25x500.dat", usecols=(1,3))
differences = []
for i in range(2,500):
  differences.append([simulation[i,0], simulation[i,1]-exact(simulation[i,0]-0.1)])
np.savetxt("differences/width20-25x500diff.dat",differences)

simulation = np.loadtxt("raw/width20-100x2000.dat", usecols=(1,3))
differences = []
for i in range(10,2000):
  differences.append([simulation[i,0], simulation[i,1]-exact(simulation[i,0]-0.095)])
np.savetxt("differences/width20-100x2000diff.dat",differences)

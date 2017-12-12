import numpy as np

nulist = [0,0.1,0.2,0.5,1,2,3,4,5]
surffield = []
for nu in nulist:
  field = np.loadtxt("nu"+str(nu)+".dat", usecols=(5,))
  if nu == 0:
    zerofield = min(field)
  maxfield = min(field)
  surffield.append([nu, maxfield, maxfield/zerofield])

np.savetxt("enhancement-factors.dat",surffield)

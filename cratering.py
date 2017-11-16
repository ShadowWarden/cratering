# cratering/cratering.py
#
# Omkar H. Ramachandran
# omkar.ramachandran@colorado.edu
#
# Cratering project for ASTR 3750 - Fall 17
#

import numpy as np
import matplotlib.pyplot as plt
import pylab

# Create matplotlib subplots
fig, ax = plt.subplots()

# Simulation parameters. Number of asteroid impacts
Nt = 200

# Upper limit of the simulated area
XMax = 500
YMax = 500

# Crater size. Keep constant at 10 for now. Eventually, make it 
# a positive gaussian
R0 = 10.
rr = np.linspace(10,50,1001)
dr = (50-10)/1001.
sigma_R = 10.
P_r = np.exp(-(rr-R0)**2/2/sigma_R**2)/np.sqrt(2*np.pi*sigma_R)*dr
# Normalize to make probabilities sum to 1
P_r /= np.sum(P_r)

# Running count of number of craters as a function of t
n = np.zeros([Nt])
# Keep track of all the centers. Initialize and append as we go along
X = np.zeros([2]) 
R = np.zeros([1])

# f(R)/f(R0) = c*r**(-2)/c*R0**(-2)

# p(x) = (alpha-1)/R0*(x/xmin)**(-alpha)

for i in range(Nt):
    n[i] = n[i-1]
    x = np.random.uniform()*XMax
    y = np.random.uniform()*YMax
    r = np.random.choice(rr,p=P_r)
    # Find creaters within the wipe condition
    if(n[i] >= 1):
        d = np.sqrt((X[:,0] - x)**2 + (X[:,1] - y)**2) 
    else:
        d = np.sqrt((X[0] - x)**2 + (X[1]-y)**2)
        d = np.array([d])
    ii = np.where(d < r+R)
    flag = 0
    for j in range(len(ii[0])-1):
#        c = 2*R0**2*np.arccos(d[ii[0][j]]*0.5/R0)-d[ii[0][j]]*0.5*np.sqrt(4*R0**2-d[ii[0][j]]**2)
        jj = ii[0][j]
        c =  R[jj-flag]**2*np.arccos((d[jj]**2-r**2+R[jj-flag]**2)/(2*d[jj]*R[jj-flag])) + r**2*np.arccos((d[jj]**2+r**2-R[jj-flag]**2)/(2*d[jj]*r))+1/2.*np.sqrt((-d[jj]+r+R[jj-flag])*(d[jj]+r-R[jj-flag])*(d[jj]-r+R[jj-flag])*(d[jj]+r+R[jj-flag]))
        if( c >= 0.5*np.pi*R[jj-flag]**2):
            # We have an obliterated crater. Get rid of row in X
            R = R[~np.logical_and(X[:,0] == X[jj-flag,0],X[:,1] == X[jj-flag,1])]
            X = X[~np.logical_and(X[:,0] == X[jj-flag,0],X[:,1] == X[jj-flag,1])]
            flag += 1

    if(n[i] == 0):
        X[0] = x
        X[1] = y
        R[0] = r
    else:
        X = np.append(X,[x,y])
        R = np.append(R,r)
    n[i] = len(X)/2
    X = X.reshape([int(n[i]),2])


for i in range(len(X)):
    circle = pylab.Circle((X[i,0],X[i,1]),radius=R[i],edgecolor='r',linestyle='--',facecolor='none')
    ax.add_artist(circle)
    

plt.xlim([0,XMax])
plt.ylim([0,YMax])
plt.xlabel(r"$x/km$")
plt.ylabel(r"$y/km$")
plt.show()


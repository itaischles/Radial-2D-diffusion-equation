import numpy as np
import matplotlib.pyplot as plt

###################################################################
# This script uses Finite Difference Method (implicit) to numerically 
# solve the partial differential equation:
# dc/dt = D(d^2c/dr^2 + 1/r dc/dr) - kc
# For details see "Summary.docx" in this folder.
###################################################################

###################################################################
# parameters
###################################################################
n_r = 100 # number of elements in radial vector
r_max = 3.0 # maximum radial distance (um)
n_t = 100 # number of elements in time vector
t_max = 7000.0 # maximum time (ps)
D = 0.1 # diffusion coefficient (cm^2/s)
k = 1/100 # decay rate (ps^-1)
sigma_r = 0.5 # initial width of Gaussian concentration profile (um)
###################################################################

# convert diffusion coefficient from cm^2/s to um^2/ps
D = D*1e-4

# set up space and time vectors, meshgrid, and step sizes
r = np.linspace(0.0,r_max,n_r)
t = np.logspace(-1,np.log10(t_max),n_t)
rgrid,tgrid = np.meshgrid(r,t) # for plotting if needed
dr = r[1]-r[0] # (constant) step size in radial position

# this is the index used in the .docx summary. The elements of r run from 0..n
# so that the number of elements in r is actually greater than n by 1.
n = len(r)-1

# initialize concentration vector to be zeros
conc = np.zeros((len(r),len(t)))
# set initial concentration profile (t=0) to a radial Gaussian
conc[:,0] = np.exp(-0.5 * r**2/sigma_r**2)
# conc[:,0] = np.ones(r.shape)
conc[1,0] = conc[0,0] # set derivative to zero at r=0
# normalize initial concentration profile by radially integrating
conc[:,0] = conc[:,0] / np.sum([i*conc[i,0] for i in range(len(r))])

###################################################################
# Using Finite Difference Method (implicit) to solve the partial differential equation.
# For details see "Summary.docx" in this folder.
###################################################################

# solve for the concentration at later times
for j in range(1,len(t)):
    
    dt = t[j]-t[j-1] # time step (changes with time)

    # define shorthand constants (also used in the .docx file)
    s = D*dt/(dr)**2
    kappa = k*dt
    
    # integrate radially to find total concentration at previous time step
    last_total_conc = np.sum([i*conc[i,j-1] for i in range(0,n)])
    
    # setting up the matrix m
    m = np.zeros((n-1, n-1))
    m[0,0:2] = [1+3*s/2+kappa, -3*s/2] # first row of m
    m[-1,-2:] = [-s+s/(2*(n-1)), 1+2*s+kappa] # last row of m
    for i in range(n-1): # add the boundary condition to the last row of m
        m[-1,i] = m[-1,i] + (i+1)/n * (s+s/(2*(n-1)))
    for i in range(1,n-2): # all other rows of m (= excluding first and last row)
        m[i,i-1:i+2] = [-s+s/(2*(i+1)), 1+2*s+kappa, -s-s/(2*(i+1))]
    
    # setting up the vector b
    b = np.zeros((n-1,1))
    b[-1] = (-s-s/(2*(n-1))) * 1/n * last_total_conc/(1+kappa)
    
    # solve for the concentrations in the boundary interior
    conc[1:-1,j] = np.matmul(np.linalg.inv(m),np.expand_dims(conc[1:-1,j-1],1)-b)[:,0]
    # solve for the boundary edges
    conc[0,j] = conc[1,j] # this is the "dc(0,t)/dr=0" boundary condition
    conc[-1,j] = 1/n * (last_total_conc/(1+kappa) - np.sum([i*conc[i,j] for i in range(0,n)]))
    
# renormalize the entire concentration surface so that the concentration at r=0,t=0 is 1
conc = conc/conc[0,0]

# integrate the concentration under the probe shape (assumed Gaussian too) to get ~deltaA(t)
deltaA = np.array([np.trapz(r*conc[:,0]*conc[:,j], r, axis=0) for j in range(len(t))])

#############################
# plot the result
#############################

fig1 = plt.figure(figsize=(15,5))

ax1 = fig1.add_subplot(131)
ax1.pcolor(rgrid,tgrid,conc.T,shading='auto',cmap=plt.cm.turbo)
ax1.contour(rgrid,tgrid,conc.T,levels=40,colors='black',alpha=0.2)
ax1.set_xlabel(r'Radial position [$\mu$m]')
ax1.set_ylabel('Time [ps]')
ax1.set_yscale('log')

ax2 = fig1.add_subplot(132)
times = [0,1,5,10,50,100,500,1000,5000]
colors = plt.cm.turbo_r(np.linspace(0.1,0.9,len(times)))
for j,time in enumerate(times):
    i = np.argmin(np.abs(t-time))
    ax2.plot(r, conc[:,i], label='t='+str(np.round(t[i])), color=colors[j])
ax2.set_xlabel(r'Radial position [$\mu$m]')
ax2.set_ylabel('Conc.')
ax2.legend()

ax3 = fig1.add_subplot(133)
radii = np.linspace(0,r_max,5)
colors = plt.cm.turbo_r(np.linspace(0.1,0.9,len(radii)))
for j,radius in enumerate(radii):
    i = np.argmin(np.abs(r-radius))
    ax3.plot(t, conc[i,:], label='r='+str(r[i]), color=colors[j])
ax3.plot(t, deltaA/max(deltaA), ls='dashed', color='black', label=r'Norm. $\Delta$A')
ax3.set_xlabel('Time [ps]')
ax3.set_ylabel('Conc.')
ax3.set_xscale('log')
ax3.legend()

plt.tight_layout()
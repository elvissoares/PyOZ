import numpy as np
import sys
sys.path.insert(0, '../src/')
from oz import rdf,uhs
import matplotlib.pyplot as plt

plt.style.use(['science'])

MCdata = np.loadtxt('data/radialdistribution_MCdata.dat',skiprows=1)
rMC, g20MC, g50MC, g90MC = MCdata[:,0],MCdata[:,1],MCdata[:,2],MCdata[:,3]

plt.scatter(rMC,g20MC,marker='o',edgecolors='C0',facecolors='none',label='MC')
plt.scatter(rMC,g50MC+1,marker='o',edgecolors='C1',facecolors='none')
plt.scatter(rMC,g90MC+2,marker='o',edgecolors='C2',facecolors='none')

dr = 0.01
r = np.arange(0.0,10.0,dr)+0.5*dr

# Testing PY
g20PY = rdf(r,0.2,u=uhs,params=1.0,model='PY')
g50PY = rdf(r,0.5,u=uhs,params=1.0,model='PY')
g90PY = rdf(r,0.9,u=uhs,params=1.0,model='PY')

plt.plot(r,g20PY,'--k',label='PY')
plt.plot(r,g50PY+1,'--k')
plt.plot(r,g90PY+2,'--k')

# Testing PY
g20HNC = rdf(r,0.2,u=uhs,params=1.0,model='HNC')
g50HNC = rdf(r,0.5,u=uhs,params=1.0,model='HNC')
g90HNC = rdf(r,0.9,u=uhs,params=1.0,model='HNC')

plt.plot(r,g20HNC,'-k',label='HNC')
plt.plot(r,g50HNC+1,'-k')
plt.plot(r,g90HNC+2,'-k')

plt.xlim(1,3)
plt.ylim(0,12)
plt.legend(loc='best')
plt.xlabel('$r/\sigma$')
plt.ylabel('$g(r)$')
plt.text(2.5,3.3,r'$\rho \sigma^3 = 0.9$')
plt.text(2.5,2.2,r'$\rho \sigma^3 = 0.5$')
plt.text(2.5,1.1,r'$\rho \sigma^3 = 0.2$')
plt.savefig('radialdistributionfunction-hardspheres.png',dpi=200)
plt.show()
plt.close()

####################################################################
# Contact Value
MCdata2 = np.loadtxt('data/radialdistribution_sigma_MCdata.dat',skiprows=1)
rhoMC, gsigmaMC = MCdata2[:,0],MCdata2[:,1]

rho = np.arange(0.01,1.0,0.05)

gsigmaPY = np.zeros_like(rho)
gsigmaHNC = np.zeros_like(rho)

for i in range(rho.size):
    gg = rdf(r,rho[i],u=uhs,params=1.0,model='PY')
    gsigmaPY[i] = gg[100]
    gg = rdf(r,rho[i],u=uhs,params=1.0,model='HNC')
    gsigmaHNC[i] = gg[100]

plt.scatter(rhoMC,gsigmaMC,marker='o',edgecolors='C0',facecolors='none',label='MC')
plt.plot(rho,gsigmaPY,'--k',label='PY')
plt.plot(rho,gsigmaHNC,'-k',label='HNC')
plt.legend(loc='best')
plt.xlabel(r'$\rho \sigma^3$')
plt.ylabel(r'$g(\sigma)$')
plt.xlim(0.0,1.0)
plt.ylim(1,6.0)
plt.savefig('contactvalue-rdf-hardspheres.png',dpi=200)
plt.show()
plt.close()

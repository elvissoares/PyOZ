import numpy as np
import sys
sys.path.insert(0, '../src/')
from oz import rdf
import matplotlib.pyplot as plt

plt.style.use(['science'])

def ulj(r,params):
    [sigma,eps] = params
    return 4*eps*((sigma/r)**12-(sigma/r)**6)

dr = 0.01
r = np.arange(0.0,10.0,dr)+0.5*dr

# plt.plot(r,ulj(r,[1.0,1.0]))
# plt.xlim(0.0,10)
# plt.ylim(-1,2)
# plt.xlabel('$r/\sigma$')
# plt.ylabel('$u(r)/\epsilon$')
# plt.text(3,0.8,r'$u(r) = 4\epsilon\left[\left(\frac{\sigma}{r}\right)^{12} -\left(\frac{\sigma}{r}\right)^{6}\right]$')
# plt.savefig('lennardjones-potential.png',dpi=200)
# plt.show()
# plt.close()

MCdata = np.loadtxt('data/radialdistribution-argon.dat',skiprows=1)
rMC, gMC = MCdata[:,0],MCdata[:,1]
plt.scatter(rMC/3.405,gMC,marker='o',edgecolors='C0',facecolors='none',label='MD')

# Testing PY
gPY = rdf(r,0.84,0.71,u=ulj,params=[1.0,1.0],model='PY')
plt.plot(r,gPY,'--k',label='PY')

# Testing PY
gHNC = rdf(r,0.84,0.71,u=ulj,params=[1.0,1.0],model='HNC')
plt.plot(r,gHNC,'-k',label='HNC')

plt.xlim(0,8)
plt.ylim(0,3.2)
plt.legend(loc='best')
plt.xlabel('$r/\sigma$')
plt.ylabel('$g(r)$')
plt.text(2.5,2.2,r'$\rho \sigma^3 = 0.84$')
plt.text(2.3,1.9,r'$k_B T/\epsilon = 0.71$')
plt.savefig('radialdistributionfunction-argon-lennardjones.png',dpi=200)
plt.show()
plt.close()
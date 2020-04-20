from pyscf import gto, dft
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import numpy as np


# basis sets to use
basis = ['STO-3g', 'cc-pVDZ', 'cc-pVTZ', 'cc-pVQZ', 'cc-pV5Z']
sets = len(basis)

# calculate energy per basis set
e = []
for i in range(len(basis)):
    print('[{}/{}]'.format(i+1, sets))
    tmp = gto.Mole()
    tmp.atom = 'O 0, 0, 0; H 0, 1, 0; H 0, 0, 1'
    tmp.basis = basis[i]
    tmp.build()
    e.append(dft.RKS(tmp).run(xc='lda,pw').e_tot)

# approximate the basis set limit y
f = lambda x, y, a, b: y + a * np.exp(-b * x)
p, _ = curve_fit(f, np.arange(1, sets+1), e)
print('Basis set limit: {} Ha'.format(p[0]))

# plot the calculated energies and the fitted curve
plt.figure()
plt.xlabel('basis functions per atomic orbital')
plt.ylabel('energies [Ha]')
plt.plot(np.arange(1, sets+1), e)
plt.plot(np.arange(1, sets, 0.1), f(np.arange(1, sets, 0.1), *p))
plt.show()

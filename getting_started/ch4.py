from pyscf import gto, dft
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import time


ch4 = gto.Mole()
ch4.atom = 'C 0.5, 0.5, 0.5; H 0, 0, 0; H 0, 1, 1; H 1, 1, 0; H 1, 0, 1'
ch4.basis = 'sto-3g'
ch4.build()

# get the calculation time per grid level
mf = dft.RKS(ch4)
mf.xc = 'lda,pw'
e_tot = []
times = []
for i in range(10):
    mf.grids.level = i
    start = time.time()
    mf.kernel()
    end = time.time()
    e_tot.append(mf.e_tot)
    times.append(end - start)

# plot energy and time per grid level
fig, ax1 = plt.subplots()
ax1.yaxis.set_major_formatter(FormatStrFormatter('%.4f'))
ax1.set_xlabel('grid level')
ax1.set_ylabel('energies [Ha]', color='blue')
ax1.plot(e_tot, color='blue')
ax2 = ax1.twinx()
ax2.set_ylabel('time [s]', color='orange')
ax2.plot(times, color='orange')
plt.show()

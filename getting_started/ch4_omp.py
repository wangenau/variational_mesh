from pyscf import lib, gto, dft
import matplotlib.pyplot as plt
import time


ch4 = gto.Mole()
ch4.atom = 'C 0.5, 0.5, 0.5; H 0, 0, 0; H 0, 1, 1; H 1, 1, 0; H 1, 0, 1'
ch4.basis = 'sto-3g'
ch4.build()

# get the calculation time per number of omp threads
mf = dft.RKS(ch4)
mf.xc = 'lda,pw'
mf.grids.level = 9
threads = range(1, 9)
times = []
for i in threads:
    lib.num_threads(i)
    start = time.time()
    mf.kernel()
    end = time.time()
    times.append(end - start)

# plot time per omp threads
plt.figure()
plt.xlabel('threads')
plt.ylabel('time [s]')
plt.plot(threads, times)
plt.show()

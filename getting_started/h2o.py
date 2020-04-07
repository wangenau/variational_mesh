from pyscf import gto, dft


h2o = gto.Mole()
h2o.atom = 'O 0 0 0; H 0 1 0; H 0 0 1'
h2o.basis = 'sto-3g'
#h2o.basis = 'dzp'
h2o.build()
# same as:
# h2o = gto.Mole()
# h2o.build(atom='O 0 0 0; H 0 1 0; H 0 0 1', basis='sto-3g')
# or:
# h2o = gto.M(atom='O 0 0 0; H 0 1 0; H 0 0 1', basis='sto-3g')

mf = dft.RKS(h2o)
mf.xc = 'lda'
mf.verbose = 4
mf.kernel()
# same as:
# mf = dft.RKS(h2o).set(xc='lda',verbose=4)
# mf.kernel()
# or:
# mf = dft.RKS(h2o).run(xc='lda',verbose=4)

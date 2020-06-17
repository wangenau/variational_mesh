from pyscf import gto
from variational_mesh.gen_mesh import *

mol = gto.Mole()
mol.atom = [
    ['O', ( 0    , 0, 0.11779)],
    ['H', ( 0, 0.75545, -0.47116)],
    ['H', ( 0, -0.75545, -0.47116)]]
mol.build()
mol.verbose = 5
mesh = gen_mesh(mol=mol, error=1e-2)
print(mesh.coords.shape)

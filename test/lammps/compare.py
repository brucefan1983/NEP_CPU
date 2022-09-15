import os
import numpy as np
from ase.io import read, write
from pynep.calculate import NEP

atoms = read('POSCAR')
atoms = atoms * (5, 5, 5)
write('C.data', atoms, format='lammps-data')
calc = NEP()
e1 = calc.get_potential_energy(atoms)
f1 = calc.get_forces(atoms)

os.system('./lmp_serial < in.lammps > out')
atoms = read('out.dump', format='lammps-dump-text')
f2 = atoms.calc.results['forces']

df = f1 - f2
print(np.max(np.abs(df)))

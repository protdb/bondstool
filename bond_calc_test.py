from bond_calc import find_bonds
from pymol import cmd

cmd.reinitialize()
cmd.fetch('1edi')
find_bonds('all', '', False, True)
cmd.save('1edi.pse')
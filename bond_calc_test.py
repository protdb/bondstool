from bond_calc import find_bonds
from pymol import cmd

cmd.reinitialize()
cmd.fetch('4hxj')
cmd.remove('not polymer')
cmd.h_add()
find_bonds('chain A and resi 217-279', '', True, True)
cmd.save('1edi.pse')
import sys
from collections import defaultdict
import __main__
__main__.pymol_argv = ['pymol', '-cq']
import pymol
pymol.finish_launching()
from pymol import cmd

input_sdf, pdb_out, mol2_out = sys.argv[1:4]
resn_name = 'PHT'

cmd.load(input_sdf, 'mol')
cmd.alter('mol', f'resn="{resn_name}"')

model = cmd.get_model('mol')
element_counter = defaultdict(int)

for atom in model.atom:
    symbol = atom.symbol.capitalize()
    element_counter[symbol] += 1
    new_name = f"{symbol}{element_counter[symbol]}"
    atom.name = new_name

cmd.delete('mol')
cmd.load_model(model, 'mol')

cmd.save(pdb_out, 'mol')
# cmd.h_add('mol')
cmd.save(mol2_out, 'mol')
cmd.quit()

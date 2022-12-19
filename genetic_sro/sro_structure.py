import numpy as np
from genetic_sro.calculators import sro_calculator
from genetic_sro.mcdft import MCDFT

def get_sro_structures(atoms,calc,n_structure,Temp,N):
    e=calc.calculate_energy(mc_step=0)
    mc_run=MCDFT(atoms,calc,e,Temp,N,traj=None)
    mc_run.build_traj()
    structures=mc_run.structures
    energy=mc_run.energy
    idx=np.argsort(energy)
    accepted_structures=[structures[i] for i in idx]
    return accepted_structures[:n_structure]
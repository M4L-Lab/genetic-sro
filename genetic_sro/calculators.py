from ase.calculators.vasp import Vasp
import os
import numpy as np
from genetic_sro.Structure import Structure
from genetic_sro.sro import get_neighbor_count, get_neighbor_count_matrix
from abc import ABC, abstractclassmethod


class calculator(ABC):
    @abstractclassmethod
    def calculate_energy(self):
        pass

    @abstractclassmethod
    def calculate_dE(self):
        pass


# class vasp_calculator(calculator):
#     def __init__(self, atoms, dir, command):
#         self.atoms = atoms
#         self.dir = dir
#         self.command = command

#     @property
#     def structure(self):
#         return self.atoms

#     @structure.setter
#     def structure(self, atoms):
#         self.atoms = atoms

#     def calculate_energy(self, mc_step):
#         calc_dir = self.dir + f"/{mc_step}"
#         calc = Vasp(
#             directory=calc_dir,
#             command=self.command,
#             istart=0,
#             ibrion=-1,
#             isif=2,
#             ialgo=48,
#             nsw=0,
#             ismear=0,
#             sigma=0.1,
#             ediff=0.1e-04,
#             prec="Normal",
#             xc="PBE",
#             lreal="Auto",
#             npar=4,
#         )
#         self.atoms.pbc = True
#         self.atoms.calc = calc
#         self.atoms.get_potential_energy()
#         return self.atoms.get_potential_energy()

#     def calculate_dE(self, Ei, Ef):
#         dE = Ef - Ei
#         return dE


class sro_calculator(calculator):
    def __init__(self, structure, target_sro, cutoffs):
        self.structure = structure
        self.cutoffs = cutoffs
        self.target_sro = target_sro

    # @property
    # def structure(self):
    #     return self.atoms

    # @structure.setter
    # def structure(self, atoms):
    #     self.atoms = atoms

    def calculate_energy(self, mc_step):
        ndatas = get_neighbor_count(self.structure,self.cutoffs)
        ndatas_nmatrix=get_neighbor_count_matrix(ndatas)
    
        e = self.get_sro_diff(ndatas_nmatrix[0])

        return e

    def get_sro_diff(self, sro):

        # return np.abs(np.abs(sro[0][0]-self.target_sro[0][0]))/np.sum(self.target_sro[0])
        return np.sum(np.abs(sro - self.target_sro)) / np.sum(self.target_sro)

    def calculate_dE(self, Ei, Ef):
        dE = Ef - Ei
        return dE

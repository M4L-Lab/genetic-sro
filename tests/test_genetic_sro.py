from collections import Counter
import pytest
import numpy as np
from genetic_sro.Structure import Structure
from genetic_sro.sro import get_neighbor_count, get_neighbor_count_matrix
from genetic_sro.sro_structure import get_sro_structures
from genetic_sro.calculators import sro_calculator
from ase.io import read


def neighbour_count_matcher(ndatas_1,ndatas_2):
    for ndata_1,ndata_2 in zip(ndatas_1,ndatas_2):
            for key in ndata_1.keys():
                for k in ndata_1[key].keys():
                    assert ndata_1[key][k]==ndata_2[key][k]



class TestAlloy:
    SA1=read('tests/CONTCAR-SA1-T100')
    SA1_alloy=Structure(SA1)
    SA1_ndata=get_neighbor_count(SA1_alloy,cutoffs=[2.872,3.00])
    target_ndata=[{
                'W': Counter({'W': 54, 'Ta': 9, 'Cr': 0}),
                'Cr': Counter({'Cr': 20, 'W': 0, 'Ta': 0}), 
                'Ta': Counter({'W': 9, 'Cr': 0})
                }]
    calc=sro_calculator(SA1_alloy,target_ndata,cutoffs=[2.872,3.00])
    SA1_target_sro_alloy=get_sro_structures(SA1_alloy,calc,n_structure=5,Temp=100,N=1000)
    SA1_target_sro_ndata=get_neighbor_count(SA1_target_sro_alloy,cutoffs=[2.872,3.00])

    @pytest.mark.parametrize('alloy, expected_ratio',
    [
        (SA1_alloy,{'W':86.0/128.0,'Cr':34.0/128.0,'Ta':8.0/128.0})
    ])
    def test_atomic_ratio(self,alloy, expected_ratio):
        for key in expected_ratio.keys():
            assert alloy.atom_ratio[key]==expected_ratio[key]

    @pytest.mark.parametrize('alloy_ndatas, expected_ndatas',
    [
        (SA1_ndata,
            [
                {
                'W': Counter({'W': 54, 'Ta': 9, 'Cr': 0}),
                'Cr': Counter({'Cr': 20, 'W': 0, 'Ta': 0}), 
                'Ta': Counter({'W': 9, 'Cr': 0})
                }, 
                
                {
                'W': Counter({'W': 410, 'Cr': 243, 'Ta': 35}), 
                'Cr': Counter({'W': 243, 'Ta': 29}), 
                'Ta': Counter({'W': 35, 'Cr': 29})
                }
            ]
        )
    ])
    def test_ndata_matrix(self,alloy_ndatas,expected_ndatas):
        neighbour_count_matcher(alloy_ndatas,expected_ndatas)
    
    @pytest.mark.parametrize('ndatas, expected_ndatas',
    [
        (SA1_ndata,
        np.array([
            [[54,0,9],
            [0,20,0],
            [9,0,0]],
            
            [[410,243,35],
            [243,0,29],
            [35,29,0]]
        ])
        )
    ])
    def test_ndata_matrix(self,ndatas,expected_ndatas):
        calc_ndatas=get_neighbor_count_matrix(ndatas)
        assert np.sum(calc_ndatas-expected_ndatas)==0

    

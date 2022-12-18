import pytest
from genetic_sro.Structure import Structure
from ase.io import read

class TestAlloy:
    SA1=read('tests/CONTCAR-SA1-T100')
    SA1_alloy=Structure(SA1)

    @pytest.mark.parametrize('alloy_name, expected_ratio',
    [
        (SA1_alloy,{'W':86.0/128.0,'Cr':34.0/128.0,'Ta':8.0/128.0})
    ])
    def test_atomic_ratio(self,alloy_name, expected_ratio):
        for key in expected_ratio.keys():
            assert alloy_name.atom_ratio[key]==expected_ratio[key]

from collections import Counter



class Structure:
    """Not implemented"""

    def __init__(self,atoms):
        self.atoms=atoms
        self.elements = list(self.atoms.symbols)
        self.atom_dict = {k: v for k, v in enumerate(self.elements)}
        self.atom_ratio = self.get_atoms_ratio()

    def get_atoms_ratio(self):
        ratio={
            key: Counter(self.elements)[key] / sum(Counter(self.elements).values())
            for key in Counter(self.elements).keys()
        }
        return ratio
'''
This is the "Material" class.
It is used to define the material properties for the components.
The properties are:
1) "name": a string of characters acting as a label
2) "sigma_u": Ultimate tensile strength expressed in [MPa]
3) "sigma_y": Yielding strength expressed in [MPa]
4) "sigma_Dm1": Material specific fatigue limit expressed in [MPa]
5) "HB": Hardness of the material
6) "sigma_Dm1C": Component specific fatigue limit expressed in [MPa]
'''
class Material:

    # Constructor
    def __init__(self, name="", sigma_u=0, sigma_y=0, sigma_Dm1=0, HB=0):
        self.name = name
        self.sigma_u = sigma_u
        self.sigma_y = sigma_y
        self.sigma_Dm1 = sigma_Dm1
        self.HB = HB
    
    # Component Fatigue Limit
    def ComponentFatigueLimit(self, section=None):
        self.sigma_Dm1C = self.sigma_Dm1*section.FLCF.C
        return self
'''
This is the "Oil" class.
It is used to define the oil properties of the oil being used.
The properties are:
1) "name": a string of characters acting as a label
2) "temperature": the oil operating temperature expressed in [degC]
3) "v1": Rated viscosity of the oil expressed in [mm^2/s]
4) "v": viscosity of the oil expressed in [mm^2/s]
5) "k": viscosity ratio of the oil
'''
class Oil:

    # Constructor
    def __init__(self, name="", temp=0, v1=0, v=0):
        self.name = name
        self.temperature = temp
        self.v1 = v1
        self.v = v
        self.l = self.v / self.v1
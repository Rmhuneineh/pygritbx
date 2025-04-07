'''
This is the "Component" class.
It's a parent class that defines general properties that are common among different components.
Other classes will inherit these properties instead of having to define every single time.
'''
class Component:

    # Constructor
    def __init__(self, name, material, axis, loc, F_tot, T_tot, omega):
        self.name = name
        self.material = material
        self.axis = axis
        self.loc = loc
        self.F_tot = F_tot
        self.T_tot = T_tot
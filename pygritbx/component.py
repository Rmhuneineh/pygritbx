'''
This is the "Component" class.
It's a parent class that defines general properties that are common among different components.
Other classes will inherit these properties instead of having to define every single time.
Properties are:
1) "name": a string of characters acting as a label
2) "material": a "Material" object of the material class defining the material properties of the component
3) "axis": a 3-element vector representing the axis along which the component is rotating with respect to a defined reference frame
4) "loc": a 3-element vector representing the location of the component with respect to a defined reference frame
5) "F_tot": the total force acting on the component expressed in [N]
6) "T_tot": the total torque acting on the component expressed in [Nm]
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
        self.omega = omega
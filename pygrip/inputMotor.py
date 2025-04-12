'''
This is the "InputMotor" class. It inherits from "Component" class.
The input motor class defines the properties of an electric motor based on:
I) Given properties:
--> 1) "name": a string of characters acting as a label
--> 2) "material": a "Material" object of the material class defining the material properties of the component
--> 3) "axis": a 3-element vector representing the axis along which the motor is rotating with respect to a defined reference frame
--> 4) "loc": a 3-element vector representing the location of the motor with respect to a defined reference frame
--> 5) "power": scalar value representing the power of the motor expressed in [W]
--> 6) "n": scalar value representing the speed expressed in [rpm]
II) Calculated properties
--> 1) "omega": a 3-vector element representing the rotational velocity of the vector expressed in [rad/s]
--> 2) "F_tot": the total force acting on the motor expressed in [N] (initialized to 0)
--> 3) "T_tot": the total torque acting on the motor expressed in [Nm] (based on power and speed)
'''
import numpy as np
from .component import Component
from .force import Force
from .torque import Torque
from math import pi

class InputMotor(Component):
    
    # Constructor
    def __init__(self, name, power, n, axis, loc):
        # Given properties
        super().__init__(name=name, material=None, axis=axis, loc=loc, F_tot=None, T_tot=None, omega=None)
        self.power = power
        self.n = n
        # Calculated properties
        self.omega = self.n * pi / 30 * self.axis
        self.T_tot = Torque(np.array([0 if o == 0 else self.power/o for o in self.omega]), self.loc)
        self.F_tot = Force(np.array([0, 0, 0]), self.loc)
        

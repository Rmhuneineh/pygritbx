'''
This is the "Torque" class.
It defines a torque vector based on two simple properties:
1) "torque": a 3-element torque vector representing the torque expressed in [Nm]
2) "loc": a 3-element vector representing the point of application of the torque expressed in [mm]
'''
import numpy as np
class Torque:

    # Constructor
    def __init__(self, torque=np.zeros(3), loc=0):
        self.torque = torque
        self.loc = loc
    
    # Overload Addition
    def __add__(obj1, obj2):
        return obj1.torque + obj2.torque
    
    # Overload Subtraction
    def __sub__(obj1, obj2):
        return obj1.torque - obj2.torque
    
    # Overload Negative
    def __neg__(obj):
        return -obj.torque
    
    # Overload Equal
    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return all(self.torque == other.torque) and all(self.loc == other.loc)
        return False
    
    # Overload Call
    def __call__(self):
        print(f"Torque: {self.torque}\nLoc: {self.loc}")

    # Calculate Magnitude
    def mag(self):
        return np.sqrt(np.sum(self.torque * self.torque))


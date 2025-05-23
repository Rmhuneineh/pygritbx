'''
This is the "Force" class.

It defines a force vector based on two simple properties:
--> 1) "force": a 3-element force vector representing the force expressed in [N]
--> 2) "loc": a scalar or 3-element vector representing the point of application of the force expressed in [mm]
'''
import numpy as np
class Force:

    # Constructor
    def __init__(self, force=np.zeros(3), loc=0):
        self.force = force
        self.loc = loc
    
    # Overload Addition
    def __add__(obj1, obj2):
        return obj1.force + obj2.force
    
    # Overload Subtraction
    def __sub__(obj1, obj2):
        return obj1.force - obj2.force
    
    # Overload Negative
    def __neg__(obj):
        return -obj.force
    
    # Overload Equal
    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return all(self.force == other.force) and all(self.loc == other.loc)
        return False
    
    # Overload Call
    def __call__(self):
        print(f"Force: {self.force}\nLoc: {self.loc}")

    # Calculate Magnitude
    def mag(self):
        return np.sqrt(np.sum(self.force * self.force))
    
    # Calculate moment at location around axis
    def moment(self, location=np.array([]), axis=np.array([]), projection=np.array([])):
        if projection.size == 0:
            if axis.size == 0:
                return np.cross(self.force, (self.loc - location) * 1e-3)
            else:
                return np.cross(self.force, (self.loc - location) * 1e-3) * np.abs(axis)
        else:
            if location.size == 0:
                return np.sum(np.cross(self.force * projection, self.loc * axis * 1e-3))
            else:
                return np.sum(np.cross(self.force * projection, location * axis * 1e-3))
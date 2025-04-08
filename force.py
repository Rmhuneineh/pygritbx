'''
This is the "Force" class.
It defines a force vector based on two simple properties:
1) "force": a 3-element force vector representing the force expressed in [N]
2) "loc": a 3-element vector representin the point of application of the force expressed in [mm]
'''
class Force:

    # Constructor
    def __init__(self, force, loc):
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
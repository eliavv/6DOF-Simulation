from numpy import *
from DynamicModel import DynamicModel
from Types import Vector

class PhysicalObject(DynamicModel):
    """description of class"""
    def __init__(self, initialPosition=Vector(0,0,0), initialVelocity=Vector(0,0,0), initialAngularVelocity=Vector(0,0,0), initialOrientation=Vector(0,0,0)):
        super(PhysicalObject, self).__init__()
        self.Position = Vector(initialPosition)
        self.Velocity = Vector(initialVelocity)
        self.Orientation = Vector(initialOrientation)
        self.AngularVelocity = Vector(initialAngularVelocity)

    
    


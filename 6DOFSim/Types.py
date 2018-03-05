from numpy import *

class Vector(object):
    """description of class"""
    INERTIAL = "INERTIAL"
    ECEF = "ECEF"
    NED = "NED"
    BODY = "BODY"

    def __init__(self, X=0.0, Y=0.0, Z=0.0, otherVector=None, coordinateFrame="INERTIAL", relativeTo='INERTIAL'):
        if (otherVector is not None):
            if (type(otherVector)==Vector):
                self.X = otherVector.X
                self.Y = otherVector.Y
                self.Z = otherVector.Z
            else:
                self.X = otherVector.flatten()[0]
                self.Y = otherVector.flatten()[1]
                self.Z = otherVector.flatten()[2]
        else:
            self.X = X
            self.Y = Y
            self.Z = Z
        self.CoordinateFrame = coordinateFrame
        self.RelativeTo = relativeTo
        self.iterIndex = 0

    def __iter__(self):
        return self

    def __next__(self):
        if self.iterIndex == 3 :
            raise StopIteration
        else:
            if (self.iterIndex == 0):
                self.iterIndex += 1
                return self.X
            if (self.iterIndex == 1):
                self.iterIndex += 1
                return self.Y
            if (self.iterIndex == 2):
                self.iterIndex += 1
                return self.Z


    def MultiplyByMatrix(self,matrix):
        if (matrix is array):
            if (matrix.shape[1] == 3):
                temp = Vector()
                temp.X = matrix[0,0]*self.X + matrix[0,1]*self.Y + matrix[0,2]*self.Z
                temp.Y = matrix[1,0]*self.X + matrix[1,1]*self.Y + matrix[1,2]*self.Z
                temp.Z = matrix[2,0]*self.X + matrix[2,1]*self.Y + matrix[2,2]*self.Z
            else:
                raise Exception('Given matrix is not 3X3!')
        raise Exception('Given matrix is not valid!')

    def __str__(self):
        return '[' + str(self.X) + ' ' + str(self.Y) + ' ' + str(self.Z) + ']\''

    def __repr__(self):
        return 'Vector[' + str(self.X) + ' ' + str(self.Y) + ' ' + str(self.Z) + ']\''

    def AddToVector(self,vec):
        if (type(vec) is Vector):           
            temp = Vector()
            temp.X = vec.X + self.X
            temp.Y = vec.Y + self.Y
            temp.Z = vec.Z + self.Z
            return temp
        else:
            raise Exception('Given parameters is not a vector!')

    def To_ndarray(self, type='column'):
        import numpy as np
        if (type=='column'):
            return np.array([[self.X], [self.Y], [self.Z]])
        if (type=='row'):
            return np.array([self.X, self.Y, self.Z])


class Angle(object):
    def __init__(self, angle_deg = None, angle_rad = None):
        if ((angle_deg==None) and (angle_rad is not None)):
            self.Radian = angle_rad
        else:
            if ((angle_deg is not None) and (angle_rad is None)):
                self.Degree = angle_deg
            else:
                self.Radian = 0.0;

    def __setattr__(self, name, value):
        if name=="Degree":
            super(Angle, self).__setattr__(name, value)
            super(Angle, self).__setattr__("Radian", deg2rad(self.Degree))
        if name=="Radian":
            super(Angle, self).__setattr__(name, value)
            super(Angle, self).__setattr__("Degree", rad2deg(self.Radian))
                
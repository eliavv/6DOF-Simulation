from Types import Vector, Angle
import math
import copy

class BaseAerodynamics(object):
    """
    Airframe class represents the aerodynamic model of a missile. 
    It consists, by default, a linear aerodynamic model but can be easilly expand to any aerosynamic model.
    
    Written by: Eliav Vaknin
    Date: 3.8.2017
    """

    def __init__(self):
        self.AeroForce = Vector()
        self.AeroForce.CoordinateFrame = self.AeroForce.BODY
        self.AeroMoment = Vector()
        self.AeroMoment.CoordinateFrame = self.AeroMoment.BODY
        self.AOA = Angle(0.0)
        self.Beta = Angle(0.0)
        self.Mach = 0.0
        self.Delta_r = Angle(0.0)
        self.Delta_p = Angle(0.0)
        self.Delta_y = Angle(0.0)
        self.Altitude_m = 0.0
        self.Reynolds = 0.0
        self.AirDensity_kg_m3 = 0.0
        self.TotalVelocity = 0.0
        self.RefLength = 0.0
        self.RefArea = 0.0

    def GetAeroForceAndMomentd(self):
        """
        This methods returns a tuple consists both the aerodynamic forces and moments vector (in this order).
        """
        forcesAndMoments = (self.AeroForce, self.AeroMoment)
        return forcesAndMoments

    def SetAeroParameters(self, **kwargs):
        '''
        This method sets the aerodynamic parameters. Its input should be a dictionary conssists the parameter name as a key and its value as the dictionary value.
        '''
        for name, value in kwargs.items():
            if name=="AOA_deg":
                self.AOA.Degree = value
            if name=="Beta_deg":
                self.Beta.Degree = value
            if name=="Mach":
                self.Mach = value
            if name=="Altitude_m":
                self.Altitude_m = value
            if name=="Reynolds":
                self.Reynolds = value


    def __calculateCoefficients(self):
        #Forces
        self.CN_alpha = 5    # 1/rad
        self.CN_delta_p = -2    # 1/rad
        self.CY_beta = self.CN  # 1/rad
        self.CY_delta_y = self.CN_delta_p  # 1/rad
        self.CA_alpha = 2
        self.CA_delta_p = self.CA_delta_y = 0.2
        
        #Moments
        self.Cm_alpha = -3    # 1/rad
        self.Cm_delta_p = 8    # 1/rad
        self.Cn_beta = 3    # 1/rad
        self.Cn_delta_y = 8    # 1/rad
        self.Cl_alphaTotal = 0    # 1/rad

    def CalculateForcesAndMoments(self):
        self.__calculateCoefficients(self)
        dynamicPressure = 0.5*self.AirDensity_kg_m3*self.TotalVelocity**2  #Calculates the dynamic pressure
        #Force coefficients calculations
        self.AeroForce.X = self.CA_alpha*self.AOA.Radian + self.CA_delta_p*self.Delta_p.Radian + self.CA_delta_y*self.Delta_y.Radian
        self.AeroForce.Y = self.CY_beta*self.Beta.Radian + self.CY_delta_y*self.Delta_y.Radian
        self.AeroForce.Z = -(self.CN_alpha*self.AOA.Radian + self.CN_delta_p*self.Delta_p.Radian)
        #Moments coefficients calculations
        self.AeroMoment.X = self.Cl_alphaTotal*math.sqrt(self.AOA.Radian**2 + self.Beta.Radian**2)
        self.AeroMoment.Y = self.Cm_alpha*self.AOA.Radian + self.Cm_delta_p*self.Delta_p.Radian
        self.AeroMoment.Z = self.Cn_beta*self.Beta.Radian + self.Cn_delta_y*self.Delta_y.Radian
        #Aerodynamic forces calculations
        self.AeroForce.X *= dynamicPressure*self.RefArea
        self.AeroForce.Y *= dynamicPressure*self.RefArea
        self.AeroForce.Z *= dynamicPressure*self.RefArea
        #Aerodynamic moments calculations
        self.AeroMoment.X *= dynamicPressure*self.RefArea*self.RefLength
        self.AeroMoment.Y *= dynamicPressure*self.RefArea*self.RefLength
        self.AeroMoment.Z *= dynamicPressure*self.RefArea*self.RefLength


from PhysicalObject import PhysicalObject
from Types import *
from skaero.atmosphere import coesa as isa
import navpy

class MissileSimplified(PhysicalObject):
    from enum import Enum
    class DynamicsTypes(Enum):
        LINEAR = 1
        NON_LINEAR = 2
        SIMPLIFIED = 3

    def __init__(self, initialPosition = Vector(0,0,0), initialVelocity = Vector(0,0,0), initialAngularVelocity = Vector(0,0,0), initialOrientation = Vector(0,0,0), 
                 FinalTime = 100.0, Mass = 300.0, Ixx = 25, Iyy = 200, RefLength = 0.3):
        super(MissileSimplified, self).__init__(initialPosition, initialVelocity, initialAngularVelocity, initialOrientation)
        self.FinalTime = FinalTime
        self.DYNAMICS_TYPE = self.DynamicsTypes.LINEAR
        self.GravityVector = Vector(0, 0, 9.81)
        # Initializes mass properties
        self.Mass = Mass
        self.Ixx = Ixx
        self.Iyy = self.Izz = Iyy
        #Initializes thrust
        self.Thrust = 0.0
        #Initializes reference sizes
        self.RefLength = RefLength
        self.RefArea = pi*RefLength*RefLength/4
        #Initializes forces and moments
        self._aeroForce = Vector(0,0,0)
        self._aeroMoment = Vector(0,0,0)
        self._thrustForce = Vector(0,0,0)
        self._thrustMoment = Vector(0,0,0)
        self.ExtForce = Vector(0,0,0)
        self.ExtMoment = Vector(0,0,0)
        self.Delta_p = 0.0
        self.Delta_y = 0.0
        # Adds the position, velocity, angular velocity and orientation to the states list
        posDict = {key : val for key, val in zip(['Position_x', 'Position_y', 'Position_z'],initialPosition)}
        velDict = {key : val for key, val in zip(['Velocity_x', 'Velocity_y', 'Velocity_z'],initialVelocity)}
        angularVelDict = {key : val for key, val in zip(['AngularVelocity_x', 'AngularVelocity_y', 'AngularVelocity_z'],initialAngularVelocity)}
        orientationDict = {key : val for key, val in zip(['Orientation_x', 'Orientation_y', 'Orientation_z'],initialOrientation)}
        self.InitialConditions = dict(posDict)
        self.InitialConditions.update(velDict)
        self.InitialConditions.update(angularVelDict)
        self.InitialConditions.update(orientationDict)
        # Creates states list
        self.CreatStateList(self.InitialConditions)
        
    # ------------ Probably can be deleted -----------------------
    def Update(self):   
        from numpy import arcsin, arccos, arctan2

        # Assigns local variables
        Vx = self.States['Velocity_x']
        Vy = self.States['Velocity_y']
        Vz = self.States['Velocity_z']
        h = -self.States['Position_z']      # Assumption : Z Position_z is the altitude
        
        self.Speed = sqrt(Vx*Vx +Vy*Vy +Vz*Vz)
        self.Alpha = arctan2(Vz, Vx)
        self.Beta = arctan2(Vy, self.Speed)

        h, T, p, self.AirDensity = isa.table(h)
        c = (331.5 + 0.6*(T - 273.15))
        self.Mach = self.Speed/c

        return super(MissileSimplified, self).Update()
    # ------------------------------------------------------------
        
    def CalcDerivatives(self):
        from numpy import sin, cos, tan
        super(MissileSimplified, self).CalcDerivatives()
        #region ----------- Assign local variables ------------------
        Px = self.States['Position_x']
        Py = self.States['Position_y']
        Pz = self.States['Position_z']
        Vx = self.States['Velocity_x']
        Vy = self.States['Velocity_y']
        Vz = self.States['Velocity_z']
        p = self.States['AngularVelocity_x']
        q = self.States['AngularVelocity_y']
        r = self.States['AngularVelocity_z']
        phi = self.States['Orientation_x']
        theta = self.States['Orientation_y']
        psi = self.States['Orientation_z']
        #endregion --------------------------------------------------
        if (self.DYNAMICS_TYPE == self.DynamicsTypes.LINEAR):
            from numpy import array
            inertialVel = self.DCM_Inertial2Body.transpose().dot(array([Vx, Vy, Vz]))
            self.Derivatives['Position_x'] = inertialVel[0]
            self.Derivatives['Position_y'] = inertialVel[1]
            self.Derivatives['Position_z'] = inertialVel[2]
            self.Derivatives['Velocity_x'] = self.ExtForce.X/self.Mass - (q*Vz - r*Vy)
            self.Derivatives['Velocity_y'] = self.ExtForce.Y/self.Mass - (r*Vx - p*Vz)
            self.Derivatives['Velocity_z'] = self.ExtForce.Z/self.Mass - (p*Vy - q*Vx)
            self.Derivatives['AngularVelocity_x'] = (self.ExtMoment.X - q*r*(self.Izz - self.Iyy))/self.Ixx
            self.Derivatives['AngularVelocity_y'] = (self.ExtMoment.Y - p*r*(self.Ixx - self.Izz))/self.Iyy
            self.Derivatives['AngularVelocity_z'] = (self.ExtMoment.Z - q*p*(self.Iyy - self.Ixx))/self.Izz
            self.Derivatives['Orientation_x'] = p + (q*sin(phi) + r*cos(phi))*tan(theta)
            self.Derivatives['Orientation_y'] = q*cos(phi) - r*sin(phi)
            self.Derivatives['Orientation_z'] = (q*sin(phi) + r*cos(phi))/cos(theta)

        if (self.DYNAMICS_TYPE == self.DynamicsTypes.NON_LINEAR):
            pass
        if (self.DYNAMICS_TYPE == self.DynamicsTypes.SIMPLIFIED):
            pass

    def Start(self):
        while (self.Time <= self.FinalTime):
            #region -------------- Update Inputs ----------------------
            self.Update()
            self.UpdateMassProperties()
            self.UpdateThrust()
            self.UpdateAero()
            self.DCM_Inertial2Body = navpy.angle2dcm(self.States['Orientation_z'], self.States['Orientation_y'], self.States['Orientation_x'])
            gravityVec = self.GravityVector.To_ndarray()
            gravityForceInBody = self.DCM_Inertial2Body.dot(gravityVec)*self.Mass
            gravityForceInBody = Vector(otherVector=gravityForceInBody)
            self.ExtForce = self._aeroForce.AddToVector(self._thrustForce).AddToVector(gravityForceInBody)
            self.ExtMoment = self._aeroMoment.AddToVector(self._thrustMoment)
            #endregion ------------------------------------------------
            #region -------------- Step -------------------------------
            self.Step()
            #endregion ------------------------------------------------
            #region -------------- Record -------------------------------
            self.Record()
            if (self.CheckTerminalConditions()):
                break
            #endregion ------------------------------------------------
    def UpdateMassProperties(self):
        time = self.Time

        if (0.0 <= time <=2.0):
            self.Mass = 300.0 + -50.0*time      #kg
            self.Ixx = 25.0 -2.5*time       #kg*m^2
            self.Iyy = self.Izz = 200.0 - 30.0*time       #kg*m^2
        if (2.0 < time <= 8.0):
            self.Mass = 100.0 - 16.666*(time - 8)   #kg
            self.Ixx = 20.0 -2.5*(time - 8)       #kg*m^2
            self.Iyy = self.Izz = 140.0 - 30.0*(time - 8)       #kg*m^2
        if (8.0 < time or time < 0.0):
            self.Mass = 100.0
            self.Ixx = 20.0
            self.Iyy = self.Izz = 140.0

    def UpdateThrust(self):
        time = self.Time
        if (0 < time <=2.0):
            self.Thrust = 40000.0   #Newtons
        if (2.0 < time <= 8.0):
            self.Thrust = 5000.0   #Newtons
        if (8.0 < time or time < 0.0):
            self.Thrust = 0.0
        self._thrustForce.X = self.Thrust

    def CheckTerminalConditions(self):
        if (self.States['Position_z'] > 2): #Hit-the-ground condition
            return True

    def UpdateAero(self):
        Alpha = self.Alpha
        Beta = self.Beta
        Delta_p = self.Delta_p
        Delta_y = self.Delta_y
        dynamicPressure = 0.5*self.AirDensity*self.Speed*self.Speed

        CA_alpha = 0.5
        CN_alpha = 2.0
        CN_delta = 0.5
        Cm_alpha = -3.0
        Cm_delta = -10.0

        CA = CA_alpha*Alpha
        CY = -CN_alpha*Beta - CN_delta*Delta_y
        CN = CN_alpha*Alpha + CN_delta*Delta_p
        
        Cl = 0.0
        Cm = Cm_alpha*Alpha + Cm_delta*Delta_p
        Cn = -Cm_alpha*Beta - Cm_delta*Delta_y

        self._aeroForce.X = -dynamicPressure*self.RefArea*CA
        self._aeroForce.Y = dynamicPressure*self.RefArea*CY
        self._aeroForce.Z = -dynamicPressure*self.RefArea*CN
        self._aeroMoment.X = dynamicPressure*self.RefArea*self.RefLength*Cl
        self._aeroMoment.Y = dynamicPressure*self.RefArea*self.RefLength*Cm
        self._aeroMoment.Z = dynamicPressure*self.RefArea*self.RefLength*Cn


if __name__ == '__main__':
    Mis = MissileSimplified(initialOrientation= Vector(0, deg2rad(60), 0))
    Mis.Start()
    table = array(Mis.RecordedTable)
    a = table[1:,:].astype(dtype='float')
    fid = open("MisTest.ono",'wb')
    savetxt(fid, a, delimiter=',', fmt='%.4f', header=str(table[0,:])[1:-1])
    fid.close()

    #region ---------------------------------- Post-processing --------------------------------------------------------
    #Imports relevant libraries
    import numpy as np
    #import matplotlib as mpl
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    #Loads the recorded data
    data = np.loadtxt("MisTest.ono", delimiter=',', skiprows=3)
    timeVec = data[:,0]

    #Plots the inertial position
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.plot(data[:,1] , data[:,2] , -data[:,3])
    ax.grid()

    # X Vs. Altitude
    fig, ax = plt.subplots(1,1)
    ax.plot(data[:,1], -data[:,3])
    ax.grid()

    # total speed Vs. time
    totalSpeed = np.sqrt(data[:,4]*data[:,4] + data[:,5]*data[:,5] + data[:,6]*data[:,6])
    fig, ax = plt.subplots(1,1)
    ax.plot(timeVec, totalSpeed)
    ax.grid()
    # alpha and beta Vs. time
    alphaVec = np.arctan2(data[:,6], data[:,4])
    betaVec = np.arcsin(data[:,5]/totalSpeed)
    fig, ax = plt.subplots(2,1)
    ax[0].plot(timeVec, alphaVec)
    ax[1].plot(timeVec, betaVec)
    ax[0].grid()
    ax[1].grid()

    plt.show()
    #endregion ---------------------------------- Post-processing --------------------------------------------------------

from numpy import *
from abc import ABCMeta, abstractmethod

class DynamicModel(object):
    """This class represents a dynamic module (module that its behavior can be modeled by differntial equations))"""
    #Enums
    EULER = "Euler"
    RUNGEKUTTA4 = "RungeKutta4"
    __metaclass__ = ABCMeta

    def __init__(self):
        """Instantiate DynamicModel class

        Keyword arguments:
        initialConditions -- Dictionary specifies the states name and initial values
        time -- Initial time in [sec] (default 0.0)
        deltaTime -- Time step size in [sec] (default 0.01)
        solver -- Specifies the solver's type [string]
        """
        self.__stateListCreated = False
        self.__recordTableInitialized = False

    def CreatStateList(self, initialConditions, initialTime=0.0, deltaTime=0.01, solver=EULER):
        """
        This method initializes the states of the module

        Keyword arguments:
        initialConditions -- Dictionary specifies the states name and initial values
        time -- Initial time in [sec] (default 0.0)
        deltaTime -- Time step size in [sec] (default 0.01)
        solver -- Specifies the solver's type [string]
        """
        self.Time = initialTime
        self.States = initialConditions.copy() #Must be dictionary type
        self.Derivatives = initialConditions.copy() #Must be dictionary type
        self.InitialConditions = initialConditions
        self.__InitalizeDerivatives()
        self.Solver = solver
        self.DeltaTime = deltaTime
        self.__stateListCreated = True

    #region ------------ Abstract Methods ----------------------

    @abstractmethod
    def CalcDerivatives(self):
        pass
    @abstractmethod
    def Update(self):
        pass

    #endregion -------------------------------------------------

    def Step(self):
        if (self.__stateListCreated==False):
            raise BaseException("State list was not created for this module")
        self.CalcDerivatives()
        if (self.Solver==self.EULER):
            self.__IncrementStatesByEuler()
        elif (self.Solver==self.RUNGEKUTTA4):
            self.IncrementStatesByRK4()
        self.Time += self.DeltaTime

    def __IncrementStatesByEuler(self):
        for stateName, value in self.States.items():
            self.States[stateName] = value + self.Derivatives[stateName]*self.DeltaTime    

    def __InitalizeDerivatives(self):
        for key, val in self.Derivatives.items():
            self.Derivatives[key] = 0.0

    def SetDerivativesByMatrices(self, A, B=None, U=None):
        """
        This method set the derivatives of the current dynamic module according to the following equation:
        X_dot = A*X + B*U
        
        where:
        X - State vector (Nx1)
        X_Dot = State vector derivatives (Nx1)
        A -- Square matrix of type ndarray (NxN)
        B -- Control matrix (NxM)
        U -- Controllers vector (Nx1)

        """
        X = array([[val] for k,val in self.States.items()])
        if (B is None):
            B = zeros(X.shape)
        if (U is None):
            U = zeros((X.shape[0],1))
        X_Dot = A*X + B*U
        i = 0
        for stateName in self.States.keys():
            self.Derivatives[stateName] = X_Dot.item(i)
            i += 1

    def Record(self):
        if (not self.__recordTableInitialized):
            self.RecordedTable = list()
            header = list(self.States.keys())
            header.insert(0,"Time[sec]")
            self.RecordedTable.append(header)
            self.__recordTableInitialized = True

        line = list(self.States.values())
        line.insert(0,self.Time)
        self.RecordedTable.append(line)
from matplotlib import pyplot as plt
import numpy as np 
import sympy as sp

class Fin_1D_Model():
    """
        This class models a Fin with a 1D temperature distribution.
        Pay attention, this class doesn't check if the error of a 1D
        aproximation is low enougth to make the results valid.
    """
    def __init__(self, physics_data, geometry_data):

        """
            Parameters
            ----------
            physics_data: dictionary
                A dictionary that must contain the keys below:
                
                - "k" : Fin's thermal conductivity [W/mK];
                - "h" : Convection heat transfer coefficient [W/Km²];
                - "T_base": Fins base temperature [K or C];
                - "T_env": Temperature of the fluid (enviroment) around the Fin [K or C];
                
                If any of the values above is a function of x, define its value as
                a function inside a string just as:
                    >> physics_data = {..., "k":"x**2+7"}

            physics_data: dictionary
                A dictionary that must contain the keys below:
                
                - "cross_area" : Fin's cross section area [m²];
                - "perimeter" : The Fin's cross section perimeter [m];
                - "length": Fins length [m];
                
                If any of the values above is a function of x, define its value as
                a function inside a string just as:
                    >> geometry_data = {..., "cross_area":"-x**3+9"}
        """

        self.x = sp.symbols('x', positive=True, real=True)
        self.T = sp.Function('T')(self.x)
        
        self.k, self.h, self.T_base, self.T_inf = self._check_data(physics_data)
        self.A, self.p, self.L = self._check_data(geometry_data)
        
        self.Lc, self.T_sol = (0,)*2

    def _check_data(self, data):
        for key in data:
            if isinstance(data[key], str):
                data[key] = sympify(data[key])
        
        return data.values()
    

    def solve(self, boundary_condition='tip_convection', tip_temperature=0):
        """
            This function solves the Fin differential equation
            and save the result in the self.T_sol propertie.

            Parameter
            ---------
            boundary_condition: string
                Defines the set of boundary conditions to be used
                by the solver. It can be:
                - 'infinitely_long_fin': 
                    The Fin is long enougth to be consider 
                    infinitely long as it's tip has the same 
                    temperature of the enviroment.

                - 'tip_convection': 
                    The Fin's tip transfers heat by convection.
                
                - 'adiabatic_tip':
                    The Fin's tip is insulated.
                
                - 'tip_defined_temperature':
                    The fin's tip temperature is known.
            
            tip_temperature: float
                Is the fin's tip temperature and must be informed 
                if the boundary_condition is "tip_defined_temperature".

        """
        fin_equation = (self.k*self.A*self.T.diff(self.x,2) - 
                            self.h*self.p*self.T + self.h*self.p*self.T_inf)

        sol_family = sp.dsolve(fin_equation)

        T_sol = sol_family.args[1]

        system_eq = []

        if boundary_condition == 'tip_convection':
            system_eq = [
                (self.k*T_sol.diff(self.x).subs({self.x:self.L}) + self.h*T_sol.subs({self.x:self.L}) 
                    - self.h*self.T_inf), 
                T_sol.subs({self.x:0})-self.T_base
            ]

        elif boundary_condition == 'infinitely_long_fin':
            T_sol = T_sol.subs({sp.symbols('C2'):0})
            system_eq = [ T_sol.subs({self.x:0}) - self.T_base]
        
        elif boundary_condition == 'adiabatic_tip':
            system_eq = [ 
                T_sol.subs({self.x:0}) - self.T_base,
                T_sol.diff(self.x).subs({self.x:self.L})
                ]
        elif boundary_condition == 'tip_defined_temperature':
            system_eq = [ 
                T_sol.subs({self.x:0}) - self.T_base,
                T_sol.subs({self.x:self.L}) - tip_temperature
                ]
        else:
            print("This class doesn't supports this set of boundary conditions")

        constants = sp.solve(system_eq)

        self.T_sol = T_sol.subs(constants)

        print('Model solved.')
    
    def get_heat_transport(self, position=0):
        """
            This function calculates the heat transfer in any point
            of the fin.

            Parameters
            ----------
            position: float
                Is the position where one wants to calculate the heat transfer.

            returns
            -------
                The heat tranfer in W.
        """

        if position > self.L:
            print("There is no Fin here. The Fin has a length of {}m".format(self.length))
        elif position < 0:
            print("The model doesn't supports negative positions.")
        else:
            return self.k*self.A*self.T_sol.diff(self.x).subs({self.x:position})

    def get_temperature(self, position=0):
        """
            This function calculates the temperature in any point
            of the fin.

            Parameters
            ----------
            position: float
                Is the position where one wants to calculate the temperature.

            returns
            -------
                The temeprature tranfer in K or C, its depends on the temperatures
                unit choosed in the physics_data dictionary.
        """
        if position > self.L:
            print("There is no Fin here. The Fin has a length of {}m".format(self.length))
        elif position < 0:
            print("The model doesn't supports negative positions.")
        else:
            return self.T_sol.subs({self.x:position})

    def get_temperature_profile(self):
        """
            Plots the Fin's temperature distribution along the x axis
        """
        fig = plt.figure()
        ax = fig.add_subplot(111)

        position = np.arange(0, self.L, self.L/1000)
        y = sp.lambdify(self.x, self.T_sol, 'numpy')

        ax.plot(position, y(position), label='Temperature Profile')

        plt.draw()
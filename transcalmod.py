from matplotlib.pyplot import draw, figure
from numpy import pi, log, arange
from sympy import symbols, Function, dsolve, solve, integrate, lambdify

def thermal_resistance_1D(resistance_dic, suppress=False):
    """     
        Calculates each thermal resistance described in the resistance_dic

        Parameters
        ----------
        resistance_dic: dictionary
            A dictionary that must contain as key the resistance id given by the user
            and as values lists of 5 elements. See the example below:
            
            Example: 
                >> example_resistance_dic = {
                    'R_id' = [mechanism, direction, area, thermal_coefficient, data],
                    'R_id2' = [mechanism, direction, area, thermal_coefficient, data],
                    ...
                }

            Values description:
            - mechanism: string
                The heat transport mechanism, this module supports:
                'contact', 'convection', 'radiation' and 'conduction'
            
            - direction: string
                If the mechanism is 'conduction', it can be:
                'axial', 'cylinder_radial' or 'sphere_radial'
                if the mechanism is other then 'conduction' it 
                can be filled with ''.

            - thermal_coefficient: float
                The thermal coefficient (k, h, h_rad or h_c) in SI
            
            - data: dictionary
                A dictionary that contains the specific informations 
                of each 'conduction' mechanism configuration. If the mechanism
                is other then 'condution', it can be a empty dictionary {}.

                For: 'direction'                'data'
                     'axial'                    {'length'}
                     'cylinder_radial'          {'angle','length','R_inner', 'R_outter'}
                     'sphere_radial'            {'R_inner', 'R_outter'}
            
        suppress: boolean, Deafult: False
            If True the function wont print the resistances values to the console
            when called.
        
        ------
        return
            A dictionary containing the resistances ids as keys and its values [K/W].
        """
    resistances = {}
    for resistance in resistance_dic:
        resistance_value = _thermal_resistance_1D_calculations(*resistance_dic[resistance])
        resistances[resistance] = resistance_value
        if not suppress:
            print('{} = {:0.4e} K/W'.format(resistance, resistance_value))
    return resistances

def _thermal_resistance_1D_calculations(mechanism, direction, area, thermal_coefficient, data):
    if mechanism == 'contact' or mechanism == 'convection' or mechanism == 'radiation':
        return 1/(area*thermal_coefficient)
    elif mechanism == 'conduction':
        if direction == 'axial':
            return data['length']/(area*thermal_coefficient)
        elif direction == 'cylinder_radial':
            coeff = 1/(thermal_coefficient*data['angle']*data['length'])
            return coeff*log(data['R_outter']/data['R_inner']) 
        elif direction == 'sphere_radial':
            radius_term = (data['R_outter']-data['R_inner'])/(data['R_outter']*data['R_inner'])
            return radius_term/(4*pi*thermal_coefficient)
        else:
            print("This module doesn't supports this direction.")
    else:
        print("This module doesn't supports this mechanism.")

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

            geometry_data: dictionary
                A dictionary that must contain the keys below:
                
                - "cross_area" : Fin's cross section area [m²];
                - "perimeter" : The Fin's cross section perimeter [m];
                - "length": Fins length [m];
                
                If any of the values above is a function of x, define its value as
                a function inside a string just as:
                    >> geometry_data = {..., "cross_area":"-x**3+9"}
        """

        self.x = symbols('x', positive=True, real=True)
        self.T = Function('T')(self.x)
        
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

        sol_family = dsolve(fin_equation)

        T_sol = sol_family.args[1]

        system_eq = []

        if boundary_condition == 'tip_convection':
            system_eq = [
                (self.k*T_sol.diff(self.x).subs({self.x:self.L}) + self.h*T_sol.subs({self.x:self.L}) 
                    - self.h*self.T_inf), 
                T_sol.subs({self.x:0})-self.T_base
            ]

        elif boundary_condition == 'infinitely_long_fin':
            T_sol = T_sol.subs({symbols('C2'):0})
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

        constants = solve(system_eq)

        self.T_sol = T_sol.subs(constants)

        print('Model solved.')
    
    def get_heat_transfer(self, position=0):
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

        if position > self.L or position < 0:
            print("Position out of the Fin")
        elif position == 0:
            return self.k*self.A*self.T_sol.diff(self.x).subs({self.x:position})
        else:
            heat_transfer = self.h*self.p*(self.T_inf-self.T_sol)
            return integrate(heat_transfer, (self.x, 0, position))
            

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
        fig = figure()
        ax = fig.add_subplot(111)

        position = arange(0, self.L, self.L/1000)
        y = lambdify(self.x, self.T_sol, 'numpy')

        ax.plot(position, y(position), label='Temperature Profile')

        draw()
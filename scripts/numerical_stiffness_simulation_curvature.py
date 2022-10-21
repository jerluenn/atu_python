import sys
sys.path.insert(0, '../../utils')
import os
import shutil

from acados_template import AcadosModel
from acados_template import AcadosOcp, AcadosOcpSolver, AcadosSimSolver, AcadosSim
from casadi import *

import numpy as np

import time
import tetherunit, rod_parameterbuilder
from matplotlib import pyplot as plt

class stiffness_simulator:

    def __init__(self, robot_dict, solution, delta_range, initConditions, no_of_simulations): 

        self.robot_dict = robot_dict
        self.boundary_length = robot_dict['tether_length'] # used for plotting
        self.integration_steps = robot_dict['integration_steps']
        self.sol = solution
        self.N = no_of_simulations
        self.initConditions = initConditions
        self.use_Fz = False
        self.use_My = True
        self.data = np.zeros((3,self.N))

        self.create_range_around_solution(delta_range)
        self.run_main()

    def set_and_integrate(self): 

        self.tetherObject._Integrator.set('x', self.initConditions)
        self.tetherObject._Integrator.solve()
        self.distalConditions = self.tetherObject._Integrator.get('x')
        # print(self.tetherObject._Integrator.get('S_forw')[7:13])

        return self.distalConditions

    def plot(self):

        pass 

    def create_range_around_solution(self, delta_range): 

        self.data_range_Fz = np.linspace(self.sol[9] - delta_range, self.sol[9] + delta_range, self.N)
        self.data_range_My = np.linspace(self.sol[9] - delta_range, self.sol[9] + delta_range, self.N)

    def save_data(self): 

        pass 

    def run_main(self): 

        if self.use_Fz:

            for i in range(self.N): 

                self.initConditions[9] = self.data_range_Fz[i]
                self.set_and_integrate()
                self.save_data()

        if self.use_My:

            for i in range(self.N): 

                self.initConditions[11] = self.data_range_My[i]
                self.set_and_integrate()
                self.save_data()
                
                 







if __name__ == "__main__":

    robot_dict = {}
    robot_dict['type'] = 'hollow_rod'
    robot_dict['outer_radius'] = 0.002
    robot_dict['inner_radius'] = 0.0006
    robot_dict['elastic_modulus'] = 1.80e9
    robot_dict['mass_distribution'] = 0.035
    robot_dict['tether_length'] = 5.1
    robot_dict['shear_modulus'] = 0.75e9
    robot_dict['integration_steps'] = 100
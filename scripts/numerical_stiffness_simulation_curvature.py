import sys
sys.path.insert(0, '../../utils')
import os
import shutil

from acados_template import AcadosModel
from acados_template import AcadosOcp, AcadosOcpSolver, AcadosSimSolver, AcadosSim
from casadi import *

import numpy as np
import pandas as pd

import time

sys.path.insert(0, '../src')

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
        self.data = np.zeros((self.N, 3))

        self.create_range_around_solution(delta_range)

        builder = rod_parameterbuilder.Rod_Parameter_Builder()
        builder.createHollowRod(robot_dict)
        self.tetherObject = tetherunit.TetherUnit(builder)

        self.run_main()

    def set_and_integrate(self): 

        self.tetherObject._Integrator.set('x', self.initConditions)
        self.tetherObject._Integrator.solve()
        self.distalConditions = self.tetherObject._Integrator.get('x')
        # print(self.tetherObject._Integrator.get('S_forw')[7:13])

        return self.distalConditions

    def plot(self):

        plt.plot(self.data_range_My, self.data[:,0])
        plt.show()
        plt.plot(self.data_range_My, self.data[:,1])
        plt.show()
        plt.plot(self.data_range_My, self.data[:,2])
        plt.show()

    def create_range_around_solution(self, delta_range): 

        self.data_range_Fz_lower = np.linspace(self.sol[9] - delta_range, self.sol[9], int(self.N/2))
        self.data_range_Fz_upper = np.linspace(self.sol[9], self.sol[9] + delta_range, int(self.N/2))

        self.data_range_Fz = np.concatenate((self.data_range_Fz_lower, self.data_range_Fz_upper))

        self.data_range_My_lower = np.linspace(self.sol[11] - delta_range, self.sol[11], int(self.N/2))
        self.data_range_My_upper = np.linspace(self.sol[11], self.sol[11] + delta_range, int(self.N/2))

        self.data_range_My = np.concatenate((self.data_range_My_lower, self.data_range_My_upper))

    def create_csv(self): 

        elastic_mod_to_weight_ratio = (self.robot_dict['elastic_modulus'])/(self.robot_dict['mass_distribution']*self.robot_dict['tether_length'])

        df = pd.DataFrame(self.data_range_Fz)
        df.to_csv('data_range_Fz_' + str(int(elastic_mod_to_weight_ratio)) + '.csv') 

        df = pd.DataFrame(self.data_range_My)
        df.to_csv('data_range_My_' + str(int(elastic_mod_to_weight_ratio)) + '.csv') 

        df = pd.DataFrame(self.data)
        df.to_csv('data_Curv_Norm_Z_' + str(int(elastic_mod_to_weight_ratio)) + '.csv')

    def run_main(self): 

        if self.use_Fz:

            for i in range(self.N): 

                self.initConditions[9] = self.data_range_Fz[i]
                self.set_and_integrate()
                self.data[i, 0] = self.distalConditions[13]
                self.data[i, 1] = np.linalg.norm([self.distalConditions[7:13]])
                self.data[i, 2] = self.distalConditions[2]

        if self.use_My:

            for i in range(self.N): 

                self.initConditions[11] = self.data_range_My[i]
                self.set_and_integrate()
                self.data[i, 0] = self.distalConditions[-1]
                self.data[i, 1] = np.linalg.norm([self.distalConditions[7:13]])
                self.data[i, 2] = self.distalConditions[2]
            
        self.plot()
        self.create_csv()


if __name__ == "__main__":

    robot_dict = {}
    robot_dict['type'] = 'hollow_rod'
    robot_dict['outer_radius'] = 0.002
    robot_dict['inner_radius'] = 0.0006
    robot_dict['elastic_modulus'] = 1.80e9
    robot_dict['mass_distribution'] = 0.035
    robot_dict['tether_length'] = 3.1
    robot_dict['shear_modulus'] = 0.75e9
    robot_dict['integration_steps'] = 50

    sol = np.array([0, 0, 0, 1, 0, 0, 0, robot_dict['tether_length']*robot_dict['mass_distribution']*9.81, -7.21548500e-26, -3.62844316e-33, 4.22730307e-26,
   0.21544033, -1.91589977e-24, 0])
    no_of_simulations = 10000
   
    stiffness_simulator(robot_dict, sol, 0.5e-1, np.copy(sol), no_of_simulations)
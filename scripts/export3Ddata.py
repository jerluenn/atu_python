import sys
import os
import shutil

from acados_template import AcadosModel
from acados_template import AcadosOcp, AcadosOcpSolver, AcadosSimSolver, AcadosSim
from casadi import *

import pandas as pd
import numpy as np
from scipy.integrate import odeint
from scipy.optimize._lsq.least_squares import least_squares

import time

from matplotlib import pyplot as plt

import sys
sys.path.insert(0, '../src')
import tetherunit, rod_parameterbuilder
from atu_inverse_solver import atu_solver

class Export3DClass: 

    def __init__(self, robot_dict): 

        self.robot_dict = robot_dict
        self.boundary_length = robot_dict['tether_length'] # used for plotting
        self.integration_steps = robot_dict['integration_steps']
        self.NUMITERATIONS = 1200
        self.initialiseSolver()
        

    def initialiseSolver(self): 

        self.solver_obj = atu_solver(self.robot_dict)
        self.solver, self.integrator = self.solver_obj.createSolver()
        self.ref = np.zeros((14, ))
        self.ref[0:7] = -self.robot_dict['tether_length']*0.8, 0.0, self.robot_dict['tether_length']*0.5, 1, 0, 0, 0

        next_step_sol = np.array([0, 0, 0, 1, 0, 0, 0, 2.80607410e-01, 6.81042090e-02,  3.66954738e-01, -1.23376586e-02, 3.90957500e-02, -8.68201654e-04, -2.23278849e-11])
        self.solver.set(0, 'x', next_step_sol)

        for i in range(robot_dict['integration_steps']): 

            self.integrator.set('x', next_step_sol)
            self.integrator.solve()
            next_step_sol = self.integrator.get('x')
            self.solver.set(i+1, 'x', next_step_sol)   


    def solve_for_initConditions(self, yref):

        yref_full = np.zeros((14, ))
        yref_full[0:7] = yref

        self.solver.set(self.robot_dict['integration_steps'], 'yref', yref_full)

        for i in range(self.NUMITERATIONS):

            self.solver.solve()

        return self.solver.get(0, 'x')

    def plot_and_save_data(self, initialConditions, n): 

        print(initialConditions)
        poseData = np.zeros((self.integration_steps + 1, 7))
        poseData[0, :] = initialConditions[0:7]
        states_i = initialConditions

        for i in range(self.integration_steps):

            self.integrator.set('x', states_i)
            self.integrator.solve()
            states_i = self.integrator.get('x')
            poseData[i+1, :] = states_i[0:7]

        df = pd.DataFrame(poseData)
        df.to_csv('poseData' + str(n+3) + '.csv') 

        ax = plt.axes(projection='3d')
        ax.set_xlim3d([0, self.boundary_length])
        ax.set_xlabel('Z')

        ax.set_ylim3d([-self.boundary_length/2, self.boundary_length/2])
        ax.set_ylabel('Y')

        # ax.set_zlim3d([0, self.boundary_length])
        ax.set_zlabel('X')
        ax.plot3D(poseData[:, 2], poseData[:, 1], -poseData[:, 0])
        plt.show()

             



if __name__ == "__main__":

    robot_dict = {}
    robot_dict['type'] = 'hollow_rod'
    robot_dict['outer_radius'] = 0.002
    robot_dict['inner_radius'] = 0.0006
    robot_dict['elastic_modulus'] = 1.0e9
    robot_dict['mass_distribution'] = 0.035
    robot_dict['tether_length'] = 3.1
    robot_dict['shear_modulus'] = 0.75e9
    robot_dict['integration_steps'] = 50

    exporter = Export3DClass(robot_dict)

    # yrefarray = np.array([[ -0.700340053384834,	-0.0377451533007010,	2.16182509811051,	0.998821053622603, -0.0181238806427636,	-0.0412531096842862,	-0.0175759724366079	],
    # [-0.753289036239369,	-1.04395020794545,	2.61988018784962,	0.998328765526116, -0.00520792646806670,	0.00483585012504650,	-0.0569515774103468	],
    # [-0.683261809964671,	-0.255748909958264,	2.36646104140311, 0.897328040761407,	-0.441012028696954,	-0.000195756832309100,	-0.0173666657900758	],
    # [-0.715080156471735,	-1.07350340558482,	2.04655547642830 , 0.661472907165279, -0.748830723301170,	0.0308061393025759	,0.0280764259660197	],
    # [-0.884783512689611,	-0.935779641720020,	1.68607770211977, 0.999755740869439,	-0.00334585082255210,	0.00385142306866910,	-0.0214673755319294],
    # [-0.728891130395668,	-0.909457144478588,	2.68134673532333,	0.999158642757614, 0.0156529175590650,	-0.0106199189847138,	-0.0357980160758634	]])

    yrefarray = np.array([                    [ -1.560975,	-0.0,	2.439024,	1, 0, 0, 0	],
                    [ -1.65,	-0.0,	1.111667,	1, 0, 0, 0	],
                    [ -1.61667,	-0.0,	2.05,	1, 0, 0, 0	]])

    for n, i in enumerate(yrefarray): 

        initConditions = exporter.solve_for_initConditions(i)
        exporter.plot_and_save_data(initConditions, n)

   



	

    
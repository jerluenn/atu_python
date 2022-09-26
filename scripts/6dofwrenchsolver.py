import queue

from traitlets import Float
import roslib
import rospy
from geometry_msgs.msg import *
from std_msgs.msg import *

import time
import numpy as np
import pyquaternion

import sys
sys.path.insert(0, '../src')

from atu_inverse_solver import atu_solver 


class wrench_solver:

    def __init__(self, num_atu, pos_w_p, pos_centroid_d, robot_dict, friction_coefficient): 

        self.num_atu = num_atu
        self.robot_dict = robot_dict
        self.friction_coefficient = friction_coefficient
        self.pos_w_centroid = np.array([-self.robot_dict['tether_length']*0.8, 0.0, self.robot_dict['tether_length']*0.5, 1, 0, 0, 0])
        self.pos_w_d = np.zeros((7, num_atu))
        self.ref = np.zeros((14, num_atu))
        self.initialise_inverse_solvers(pos_w_p)
        self.initialise_6dof_sensor(pos_centroid_d) 

    def initialise_inverse_solvers(self, pos_w_p): 

        self.inverse_solvers_list = []
        self.integrators_list = []


        for i in range(self.num_atu):
             
            solver_obj = atu_solver(self.robot_dict)
            solver, integrator = solver_obj.createSolver() 
            self.inverse_solvers_list.append(solver)
            self.integrators_list.append(integrator)

            next_step_sol = np.array([0, 0, 0, 1, 0, 0, 0, 1.35202744e-01,  8.59444117e-11,  1.38997104e-01, -3.21851497e-11,  6.09179901e-03, -5.40886376e-12, 0])
            next_step_sol[0:7] = pos_w_p[i]
            self.solver.set(0, 'x', next_step_sol)

            for k in range(robot_dict['integration_steps']): 

                self.integrators_list[i].set('x', next_step_sol)
                self.integrators_list[i].solve()
                next_step_sol = self.integrators_list[i].get('x')
                self.inverse_solvers_list[i].set(k+1, 'x', next_step_sol)   

    def initialise_6dof_sensor(self, pos_centroid_d):

        NUM_ITERATIONS = 1000

        for i in range(self.num_atu):

            self.pos_w_d[0:3, i] = self.pos_w_centroid + pos_centroid_d[0:3, i]
            self.ref[0:7, i] = self.pos_w_d[0:7, i]
            self.inverse_solvers_list[i].set(self.ref[:, i])

            for k in range(NUM_ITERATIONS):

                self.inverse_solvers_list[i].solve()
            

    def initialiseSubscribers(self):

        rospy.Subscriber('/relative_pose', PoseStamped, self.relative_pose_callback)

    def solve_structure_problem(self): 

        pass 

    def relative_pose_callback(self):

        pass 

    
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

    num_atu = 4 
    ground_positions = np.array([[0., -0.15, 0., 1, 0, 0, 0],[0., 0., 0., 1, 0, 0, 0],[0., 0.15, 0., 1, 0, 0, 0],[0., 0.30, 0., 1, 0, 0, 0]])
    centroid_distal_positions = np.array([[0.027, 0.015, -0.013, 1, 0, 0, 0], [-0.027, 0.015, -0.013, 1, 0, 0, 0], [0.027, -0.015, -0.013, 1, 0, 0, 0], [-0.027, -0.015, -0.013, 1, 0, 0, 0]])
    wrench_solver(num_atu, ground_positions, centroid_distal_positions, robot_dict, 0.0569)
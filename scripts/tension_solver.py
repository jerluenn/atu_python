import roslib
import rospy
from geometry_msgs.msg import *
from std_msgs.msg import *

import time
import numpy as np

import sys
sys.path.insert(0, '../src')

from atu_inverse_solver import atu_solver 

class tension_solver: 

    def __init__(self, robot_dict): 

        self.robot_dict = robot_dict
        self.initialiseVariables()
        self.initialiseSolver()
        self.initialiseNode()

    def initialiseVariables(self):

        self.curvature_msg = Float64()
        self.pose_filtered = PoseStamped()
        self.pose_estimated = PoseStamped()
        self.pose_estimated_numpy = np.zeros((7, ))
        self.pose_filtered_numpy = np.zeros((7, ))
        self.endeffector_pose = np.zeros((7, ))
        self.median_pose = np.zeros((7, 7))
        self.median_i = 0

    def initialiseSolver(self): 

        self.solver_obj = atu_solver(self.robot_dict)
        self.solver, self.integrator = self.solver_obj.createSolver()
        self.ref = np.zeros((14, ))
        self.ref[0:7] = -self.robot_dict['tether_length']*0.8, 0.0, self.robot_dict['tether_length']*0.5, 1, 0, 0, 0

        next_step_sol = np.array([0, 0, 0, 1, 0, 0, 0, 1.35202744e-01,  8.59444117e-11,  1.38997104e-01, -3.21851497e-11,  6.09179901e-03, -5.40886376e-12, 0])
        self.solver.set(0, 'x', next_step_sol)

        for i in range(robot_dict['integration_steps']): 

            self.integrator.set('x', next_step_sol)
            self.integrator.solve()
            next_step_sol = self.integrator.get('x')
            self.solver.set(i+1, 'x', next_step_sol)   

    def initialiseNode(self):

        rospy.init_node('tension_solver')

        self.initialiseSubscribers()
        self.initialisePublishers()
        self.main_loop()

        rospy.spin()

    def initialiseSubscribers(self):

        rospy.Subscriber('/loadcell_measurements', Float64MultiArray, self.loadcell_measurements_callback)
        rospy.Subscriber('/relative_pose', PoseStamped, self.relative_pose_callback)

    def initialisePublishers(self): 

        self.estimated_tension_publisher = rospy.Publisher('/estimated_tension', Float64, queue_size=1)
        self.estimated_pose_publisher = rospy.Publisher('/estimated_pose', PoseStamped, queue_size=1)
        self.curvature_publisher = rospy.Publisher('/curvature_publisher', Float64, queue_size=1)
        self.filtered_pose_publisher = rospy.Publisher('/filtered_pose', PoseStamped, queue_size=1)

    def medianFilterStep(self): 

        self.median_pose[self.median_i%7, :] = self.endeffector_pose
        self.median_i += 1
        self.pose_filtered_numpy[0:7] = np.median(self.median_pose[:, 0]), np.median(self.median_pose[:, 1]), np.median(self.median_pose[:, 2]),\
        np.median(self.median_pose[:, 3]), np.median(self.median_pose[:, 4]), np.median(self.median_pose[:, 5]), np.median(self.median_pose[:, 6])
        self.pose_filtered.pose.position.x = self.pose_filtered_numpy[0]
        self.pose_filtered.pose.position.y = self.pose_filtered_numpy[1]
        self.pose_filtered.pose.position.z = self.pose_filtered_numpy[2]
        self.pose_filtered.pose.orientation.w = self.pose_filtered_numpy[3]
        self.pose_filtered.pose.orientation.x = self.pose_filtered_numpy[4]
        self.pose_filtered.pose.orientation.y = self.pose_filtered_numpy[5]
        self.pose_filtered.pose.orientation.z = self.pose_filtered_numpy[6]

    def loadcell_measurements_callback(self, data):

        pass 

    def relative_pose_callback(self, data):

        self.endeffector_pose[0] = data.pose.position.x
        self.endeffector_pose[1] = data.pose.position.y
        self.endeffector_pose[2] = data.pose.position.z
        self.endeffector_pose[3] = np.abs(data.pose.orientation.w)
        self.endeffector_pose[4] = data.pose.orientation.x
        self.endeffector_pose[5] = data.pose.orientation.y
        self.endeffector_pose[6] = data.pose.orientation.z

        self.medianFilterStep()

    def publishData(self): 

        self.estimated_pose_publisher.publish(self.pose_estimated)
        # self.estimated_tension_publisher.publish()
        self.filtered_pose_publisher.publish(self.pose_filtered)
        self.curvature_publisher.publish(self.curvature_msg)

    def updateSolver(self):


        self.ref[0:7] = self.pose_filtered_numpy
        self.solver.set(self.robot_dict['integration_steps'], 'yref', self.ref)
        self.solver.solve()
        self.pose_estimated_numpy[0:7] = self.solver.get(self.robot_dict['integration_steps'], "x")[0:7]
        self.pose_estimated.pose.position.x = self.pose_estimated_numpy[0]
        self.pose_estimated.pose.position.y = self.pose_estimated_numpy[1]
        self.pose_estimated.pose.position.z = self.pose_estimated_numpy[2]
        self.pose_estimated.pose.orientation.w = self.pose_estimated_numpy[3]
        self.pose_estimated.pose.orientation.x = self.pose_estimated_numpy[4]
        self.pose_estimated.pose.orientation.y = self.pose_estimated_numpy[5]
        self.pose_estimated.pose.orientation.z = self.pose_estimated_numpy[6]
        self.curvature_msg.data = self.solver.get(self.robot_dict['integration_steps'], "x")[13]
        
        print(self.solver.get(self.robot_dict['integration_steps'], "x"))

    def main_loop(self): 

        while not rospy.is_shutdown(): 

            t0 = time.time()

            self.updateSolver()
            self.publishData()

            time.sleep(0.0001) 

            print(f"Total time (s): {time.time() - t0}")



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

    tension_solver(robot_dict)
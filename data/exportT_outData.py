import sys
import os
import shutil

from acados_template import AcadosModel
from acados_template import AcadosOcp, AcadosOcpSolver, AcadosSimSolver, AcadosSim
from casadi import *

import roslib
import rospy
from geometry_msgs.msg import *
from std_msgs.msg import *

import pandas as pd
import numpy as np
from scipy.integrate import odeint
from scipy.optimize._lsq.least_squares import least_squares

import time

from matplotlib import pyplot as plt

import sys
sys.path.insert(0, '../src')

class T_out_exporter: 

    def __init__(self): 

        self.initialiseVariables()
        self.initialiseSolver()
        self.initialiseNode()

    def initialiseVariables(self):

        self.T_out = Float64MultiArray() 
        self.friction_coefficient = 0.0569
        self.loadcell_measurements = np.zeros((2, ))
        self.T_out_numpy = np.zeros((2, ))
        self.curvature = 12.25
        self.K_constant = 0.4981

    def initialiseSolver(self):

        pass

    def initialiseNode(self):

        rospy.init_node('tension_solver')

        self.initialiseSubscribers()
        self.initialisePublishers()
        self.main_loop()

        rospy.spin()

    def initialiseSubscribers(self):

        rospy.Subscriber('/loadcell_measurements', Float64MultiArray, self.loadcell_measurements_callback)
        rospy.Subscriber('/curvature', Float64, self.curvature_callback)

    def initialisePublishers(self):

        self.estimated_tension_publisher = rospy.Publisher('/estimated_tension', Float64MultiArray, queue_size=1)

    def loadcell_measurements_callback(self, data): 

        self.loadcell_measurements[0:2] = data.data
        self.compute_estimated_tension()

    def curvature_callback(self, data): 

        self.curvature = data.data

    def compute_estimated_tension(self): 

        if self.loadcell_measurements[1] < 0: 

            self.T_out_numpy[0:2] = 0, 0

        else: 

            self.T_out_numpy[0] = self.loadcell_measurements[1]*self.K_constant
            self.T_out_numpy[1] = self.loadcell_measurements[1]*np.exp(-self.friction_coefficient*self.curvature)
            self.T_out.data = self.T_out_numpy


    def publishData(self): 

        self.estimated_tension_publisher.publish(self.T_out)

    def main_loop(self): 

        while not rospy.is_shutdown(): 

            self.publishData()

            time.sleep(0.001) 

if __name__ == "__main__":

    T_out_exporter()
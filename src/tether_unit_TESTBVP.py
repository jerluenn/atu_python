import sys
sys.path.insert(0, '../../utils')
import os
import shutil

from acados_template import AcadosModel
from acados_template import AcadosOcp, AcadosOcpSolver, AcadosSimSolver, AcadosSim
from casadi import *

import numpy as np
from scipy.integrate import odeint
from scipy.optimize._lsq.least_squares import least_squares

import time
import tetherunit, rod_parameterbuilder
from matplotlib import pyplot as plt


class TetherUnitBoundarySolver: 

    def __init__(self, robot_dict, initConditions, distalPose): 

        self.build_tetherObject(robot_dict)
        self.boundary_length = robot_dict['tether_length'] # used for plotting
        self.integration_steps = robot_dict['integration_steps']
        self.distalPose = distalPose 
        self.initConditions = initConditions
        self.distalConditions = None 

    def build_tetherObject(self, robot_dict): 

        builder = rod_parameterbuilder.Rod_Parameter_Builder()
        builder.createHollowRod(robot_dict)
        try:
            self.tetherObject = tetherunit.TetherUnit(builder, sys.argv[1])
        except: 
            self.tetherObject = tetherunit.TetherUnit(builder)

    def set_and_integrate(self): 

        self.tetherObject._Integrator.set('x', self.initConditions)
        self.tetherObject._Integrator.solve()
        self.distalConditions = self.tetherObject._Integrator.get('x')
        print(self.tetherObject._Integrator.get('S_forw')[7:13])

        return self.distalConditions

    def set_integrate_eigenvalues(self): 

        self.poseData = np.zeros((self.integration_steps + 1, 7))
        self.internalWrenchData = np.zeros((self.integration_steps + 1, 6))
        self.arcData = np.zeros((self.integration_steps + 1, 1))
        self.poseData[0, :] = self.initConditions[0:7]
        self.internalWrenchData[0, :] = self.initConditions[7:13]
        self.eigenvaluesData = np.zeros((self.integration_steps + 1, 17), dtype = 'complex_')
        states_i = self.initConditions
        self.eigenvaluesData[0, :] = np.linalg.eig(self.tetherObject.linearisedEquationsFunction(states_i).full())[0]


        for i in range(self.integration_steps): 
            
            self.tetherObject._stepIntegrator.set('x', states_i)
            self.tetherObject._stepIntegrator.solve()
            states_i = self.tetherObject._stepIntegrator.get('x')
            self.eigenvaluesData[i + 1, :] = np.linalg.eig(self.tetherObject.linearisedEquationsFunction(states_i).full())[0]
            self.poseData[i + 1, :] = states_i[0:7]
            self.arcData[i + 1, :] = states_i[14]
            self.internalWrenchData[i + 1, :] = states_i[7:13]
            


    def solveBVP(self, debug, plot): 

        t_start = time.time()
        # sol = least_squares(self.getResiduals, np.zeros(12), method = 'lm')
        initial_guess = np.array([1.0*9.81*self.tetherObject._tether_length*self.tetherObject._mass_distribution, 0.0001, 0.0001, 0.00001, 0.185, 0.0001])
        lower_bound = np.array([0.9999*9.81*self.tetherObject._tether_length*self.tetherObject._mass_distribution, -0.0001, -0.0001, -0.0001, 0.18, -0.0001])
        upper_bound = np.array([1.0001*9.81*self.tetherObject._tether_length*self.tetherObject._mass_distribution, 0.0001, 0.0001, 0.0001, 0.22, 0.0001])
        sol = least_squares(self.getResiduals, initial_guess, bounds = (lower_bound, upper_bound) ,method = 'trf')
        self.initConditions[7:13] = sol.x[0:6] 
        self.set_and_integrate()

        if debug == True: 

            print(f"Internal forces and moments at proximal end: {sol.x[0:6]}.")
            print(f"External wrench at end: {sol.x[6:12]}")
            print(f"Internal forces and moments at distal end: {self.distalConditions[7:13]}")
            print(f"Pose at distal end: {self.distalConditions[0: 3]}")
            print(f"Tension at proximal end: {self.initConditions[14]}")
            print(f"Tension at distal end: {self.distalConditions[14]}")
            print(f"Total curvature: {self.distalConditions[16]}")
            print(f"Total cost: {sol.cost}")
            print(f"Time taken: {time.time() - t_start} s.")
            

        if plot == True: 

            self.plotData(debug)

    def plotData(self, debug):
 
        self.poseData = np.zeros((self.integration_steps + 1, 7))
        self.internalWrenchData = np.zeros((self.integration_steps + 1, 6))
        self.arcData = np.zeros((self.integration_steps + 1, 1))
        self.poseData[0, :] = self.initConditions[0:7]
        self.internalWrenchData[0, :] = self.initConditions[7:13]
        states_i = self.initConditions

        for i in range(self.integration_steps): 
            
            self.tetherObject._stepIntegrator.set('x', states_i)
            t = time.time()
            self.tetherObject._stepIntegrator.solve()
            states_i = self.tetherObject._stepIntegrator.get('x')
            # print(f"One integration step of a 1 step: {time.time() - t}")
            self.poseData[i + 1, :] = states_i[0:7]
            self.arcData[i + 1, :] = states_i[14]
            self.internalWrenchData[i + 1, :] = states_i[7:13]

            if debug == True: 

                pass
                # print("States at node", i, ": ", states_i)

        if debug == True: 

            for i in range(7): 

                if i == 0: 

                    plt.plot(self.arcData, self.poseData[:, 0:3])
                    plt.show()
                    plt.plot(self.arcData, self.poseData[:, 3:7])
                    plt.show()

                else: 

                    plt.plot(self.arcData, self.internalWrenchData[:, i - 1])
                    plt.show()

            

        ax = plt.axes(projection='3d')
        ax.set_xlim3d([0, self.boundary_length])
        ax.set_xlabel('Z')

        ax.set_ylim3d([-self.boundary_length/2, self.boundary_length/2])
        ax.set_ylabel('Y')

        # ax.set_zlim3d([0, self.boundary_length])
        ax.set_zlabel('X')
        ax.plot3D(self.poseData[:, 2], self.poseData[:, 1], -self.poseData[:, 0])
        plt.show()

        

    def getResiduals(self, guess): 

        external_wrench = np.zeros(6)
        self.initConditions[7:13] = guess[0:6]
        self.set_and_integrate()
        # residual = np.hstack((50*(self.distalConditions[0:3] - self.distalPose[0:3]), 50*(self.distalConditions[3:7] - self.distalPose[3:7])))
        residual = np.hstack(50*(self.distalConditions[7:13] - external_wrench) + 50*self.distalConditions[16])
        # print(norm_2(residual))

        return residual


def test_Function2(testClass, initConditions):

    testClass.initConditions = initConditions

    t = time.time() 

    testClass.set_and_integrate()

    print("Time taken for one integration step: ", time.time() - t)

    testClass.plotData(True)
    plt.show()
    testClass.initConditions = testClass.distalConditions
    print(testClass.distalConditions)

def solveIteratively(testClass, initConditions, numIterations): 

    testClass.initConditions = initConditions 
    rodData = np.zeros((11, 17))
    rodData[0, :] = initConditions

    for i in range(numIterations): 

        data = testClass.set_and_integrate()
        testClass.initConditions = data 
        rodData[i+1, :] = data 

    print("Done")

    ax = plt.axes(projection='3d')
    ax.set_xlabel('Z')
    ax.set_xlim3d([0, 6])

    ax.set_ylabel('Y')
    ax.set_ylim3d([-1.5, 1.5])

    # ax.set_zlim3d([0, self.boundary_length])
    ax.set_zlabel('X')
    ax.plot3D(rodData[:, 2], rodData[:, 1], -rodData[:, 0])
    plt.show()

def getEigenvalues(testClass, initConditions):

    testClass.initConditions = initConditions
    testClass.set_integrate_eigenvalues()   
    
    return testClass.eigenvaluesData

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

    initConditions = np.array([0, 0, 0, 1, 0, 0, 0, robot_dict['tether_length']*robot_dict['mass_distribution']*9.81*10, -7.21548500e-26, -3.62844316e-33, 4.22730307e-26,
   0.21544033, -1.91589977e-24, 5, 0, 0.05, 0]) 
    initConditions = np.array([0, 0, 0, 1, 0, 0, 0, robot_dict['tether_length']*robot_dict['mass_distribution']*9.81, -7.21548500e-26, -3.62844316e-33, 4.22730307e-26,
   0.21544033, -1.91589977e-24, 5, 0, 0.05, 0]) 
    initConditions = np.array([0, 0, 0, 1, 0, 0, 0, robot_dict['tether_length']*robot_dict['mass_distribution']*9.81,  5.394324155e-06,    0.00, -3.025669809e-06,  0.278456471, -3.501944344e-07, 0]) 
    distalPose = np.array([-0.6, 0, 0.485, 1, 0, 0, 0])
    testClass = TetherUnitBoundarySolver(robot_dict, initConditions, distalPose)
    # solveIteratively(testClass, initConditions, 10)
    test_Function2(testClass, testClass.initConditions)
    eig = getEigenvalues(testClass, testClass.initConditions)
    print(testClass.distalConditions[16])
    # print(eig)


    # print(testClass.tetherObject._Kbt)

    # a = testClass.tetherObject.linearisedEquationsFunction(testClass.distalConditions)
    # b = np.linalg.eig(a.full())
    # print(b[0])






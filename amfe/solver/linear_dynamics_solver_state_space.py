#
# Copyright (c) 2018 TECHNICAL UNIVERSITY OF MUNICH, DEPARTMENT OF MECHANICAL ENGINEERING, CHAIR OF APPLIED MECHANICS,
# BOLTZMANNSTRASSE 15, 85748 GARCHING/MUNICH, GERMANY, RIXEN@TUM.DE.
#
# Distributed under 3-Clause BSD license. See LICENSE file for more information.
#

"""
Super class of all state-space linear dynamics solvers.
"""

import numpy as np
from time import time

from .solver import Solver
from ..linalg.linearsolvers import PardisoLinearSolver
__all__ = [
    'LinearDynamicsSolverStateSpace'
]


class LinearDynamicsSolverStateSpace(Solver):
    '''
    General class for solving the linear dynamic problem of the state-space system linearized around zero-state.

    Parameters
    ----------
    mechanical_system : Instance of MechanicalSystemStateSpace
        State-space system to be solved.
    options : Dictionary
        Options for solver.

    References
    ----------
       [1]  M. Géradin and D.J. Rixen (2015): Mechanical vibrations. Theory and application to structural dynamics.
            ISBN 978-1-118-90020-8.
    '''

    def __init__(self, mechanical_system, **options):
        self.mechanical_system = mechanical_system
        self.E = self.mechanical_system.E()
        self.A = self.mechanical_system.A()

        # read options
        if 'linear_solver' in options:
            self.linear_solver = options['linear_solver']
        else:
            # TODO: What is a good choice for constrained state space systems?
            print('Attention: No linear solver object was given, setting linear_solver = PardisoSolver(...).')
            self.linear_solver = PardisoLinearSolver(A=None, mtype='nonsym')

        if ('initial_conditions' in options) and ('x0' in options['initial_conditions']):
            x0 = options['initial_conditions']['x0']
            # TODO: The following section is commented out because this prevents solving reduced mechanical systems
            # TODO: because these systems do not have a ndof property.
            # if len(x0) != 2*self.mechanical_system.dirichlet_class.no_of_constrained_dofs:
            #     raise ValueError('Error: Dimension of x0 is not valid for state-space system.')
        else:
            print('Attention: No initial state was given, setting x0 = 0.')
            x0 = np.zeros(2 * self.mechanical_system.dirichlet_class.no_of_constrained_dofs)
        self.initial_conditions = {'x0': x0}

        if 't0' in options:
            self.t0 = options['t0']
        else:
            print('Attention: No initial time was given, setting t0 = 0.0.')
            self.t0 = 0.0

        if 't_end' in options:
            self.t_end = options['t_end']
        else:
            print('Attention: No end time was given, setting t_end = 1.0.')
            self.t_end = 1.0

        if 'dt' in options:
            self.dt = options['dt']
        else:
            print('Attention: No time step size was, setting dt = 1.0e-4.')
            self.dt = 1.0e-4

        if 'output_frequency' in options:
            self.output_frequency = options['output_frequency']
        else:
            print('Attention: No output frequency was given, setting output_frequency = 1.')
            self.output_frequency = 1

        if 'verbose' in options:
            self.verbose = options['verbose']
        else:
            print('Attention: No verbose was given, setting verbose = False.')
            self.verbose = False
        return

    def overwrite_parameters(self, **options):
        pass

    def effective_stiffness(self):
        pass

    def effective_force(self, x_old, dx_old, t, t_old):
        pass

    def update(self, x, x_old, dx_old):
        pass

    def solve(self):
        '''
        Solves the linear dynamic problem of the state-space system linearized around zero-state.
        '''

        # start time measurement
        t_clock_start = time()

        # initialize variables and set parameters
        self.mechanical_system.clear_timesteps()
        t = self.t0
        x = self.initial_conditions['x0'].copy()

        # evaluate initial derivative
        self.linear_solver.set_A(self.E)
        dx = self.linear_solver.solve(self.A @ x + self.mechanical_system.F_ext(x, t))

        # LU-decompose effective stiffness
        K_eff = self.effective_stiffness()
        self.linear_solver.set_A(K_eff)
        if hasattr(self.linear_solver, 'factorize'):
            self.linear_solver.factorize()

        # write output of initial conditions
        self.mechanical_system.write_timestep(t, x.copy())

        # time step loop
        output_index = 0
        while t < self.t_end:

            # save old variables
            x_old = x.copy()
            dx_old = dx.copy()
            t_old = t

            # solve system
            output_index += 1
            t += self.dt
            F_eff = self.effective_force(x_old, dx_old, t, t_old)

            x = self.linear_solver.solve(F_eff)

            # update variables
            dx = self.update(x, x_old, dx_old)

            # write output
            if output_index == self.output_frequency:
                self.mechanical_system.write_timestep(t, x.copy())
                output_index = 0

            print('Time: {0:3.6f}.'.format(t))

            # end of time step loop

        self.linear_solver.clear()

        # end time measurement
        t_clock_end = time()
        print('Time for solving linear dynamic problem: {0:6.3f} seconds.'.format(t_clock_end - t_clock_start))
        return

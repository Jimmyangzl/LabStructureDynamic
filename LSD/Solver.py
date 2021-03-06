# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 14:04:42 2019

@author: q446161
"""

import numpy as np 

from scipy.sparse import linalg  as LAS
from scipy import linalg  as LA

import time

class Solver:
    def __init__(self,M,K):
        # The solver class consist of a mass and stiffness matrix
        self.M = M
        self.K = K
        
        # Each solution procedure has an own object
        # The classes are defined below
        self.Static = StaticSolver(self.K)
        self.Eigen  = Eigensolver(self.M, self.K)
        self.Transient  = TransientSolver(self.M, self.K)        
        self.Harmonic   = HarmonicSolver(self.M, self.K)
        
        
class StaticSolver:
    def __init__(self,K):
        # The StaticSolver only depends on the stiffness matrix
        self.K = K

         
    def Inverse(self, load = 'all', N = 1):  
        # load: right hand side vector of the static problem
        # N: number of repeated evaluations for time measurement
        if load == 'all':
            Load = np.eye(np.size(self.K,1))   

        else:
            Load = load
            
        K = self.K
        
        # Compute the full inverse of the stiffness matrix using LAS (scipy.sparse.linalg).
        # The function your need is called inv().
        # Use a for loop to compute the inverse multiple times 
        # for a reliable time-measurement
        start = time.time()
        #++++++++++++ Something is missing here +++++++++++++++
        for i in range(N):
            K_inv = LAS.inv(K)
        # Calculation of inverse is the 'offline' time
        end = time.time()
        tOffline = (end-start)/N    
        

        # Compute the static response by multiplying the inverse matrix on the
        # load vector.
        # Use a for loop to perform the multiplication multiple times 
        # for a reliable time-measurement
        start = time.time()
        #++++++++++++ Something is missing here +++++++++++++++
        for i in range(N):
            StaticResponse = K_inv @ Load;
        # The evaluation of the matrix multiplication is the 'online' time
        end = time.time()
        tOnline = (end-start)/N    
    
        # Return the static response, the offline and online time
        
        return StaticResponse, tOffline, tOnline
    

    
    def Factorization(self, load = 'all', N = 1):
        # load: right hand side vector of the static problem
        # N: number of repreated evaluations for time measurement
               
        if load == 'all':
            Load = np.eye(np.size(self.K,1))   

        else:
            Load = load
        K = self.K
        # Compute the factorization of the stiffness matrix using LAS (scipy.sparse.linalg)
        # The function your need is called splu(). This function returns and object.
        # Use a for loop to compute the inverse multiple times 
        # for a reliable time-measurement
       
        start = time.time()
        #++++++++++++ Something is missing here +++++++++++++++
        for i in range(N):
            K_splu = LAS.splu(K) # returns an object, which has a solve method.
        # Calculation of inverse is the 'offline' time
        end = time.time()
        tOffline = (end-start)/N    

        # Compute the static response by solving the factorized problem.
        # Use the solve() command of the SPLU object.
        # Use a for loop to perform the solution multiple times 
        # for a reliable time-measurement
        start = time.time()
        #++++++++++++ Something is missing here +++++++++++++++
        for i in range(N):
            StaticResponse = K_splu.solve(Load)
     # The evaluation of the matrix factorization is the 'online' time
        end = time.time()
        tOnline = (end-start)/N    
    
        # Return the static response, the offline and online time
        return StaticResponse, tOffline, tOnline


class TransientSolver:
    def __init__(self,M,K):
        # The transient solver requires the stiffness and mass matrix
        self.M = M
        self.K = K
        # The maximum frequency is required to calculate the stability limits
        Eig_solver = Eigensolver(M,K)
        self.omegaMax = Eig_solver.MaxOmega()
               
    def stabilityLimits(self):    
        
        return {"PurelyExplicit": 0,
        "CentralDifference": 2 / self.omegaMax,
        "Fox_N_Goodwin": 2.45 / self.omegaMax,
        "LinearAcceleration": 3.46 / self.omegaMax,
        "AverageConstantAcceleration": "unlimited"}
           
    def PurelyExplicit(self,x0,xd0,p,t,alpha=0,beta=0):
        # Interface to the purely explicit Newmark integration
        # Constants 0, 0 
        x, xd, xdd = self.LinearNewmarkIntegration(x0,xd0,p,t,0,0,alpha,beta) 
        return x, xd, xdd
    
    def CentralDifference(self,x0,xd0,p,t,alpha=0,beta=0):
        # Interface to the central difference Newmark integration
        # Constants 1/2, 0 
        x, xd, xdd = self.LinearNewmarkIntegration(x0,xd0,p,t,0.5,0,alpha,beta) 
        return x, xd, xdd  
      
    def Fox_N_Goodwin(self,x0,xd0,p,t,alpha=0,beta=0):
        # Interface to the Fox and Goodwin Newmark integration
        # Constants 1/2, 1/12       
        x, xd, xdd = self.LinearNewmarkIntegration(x0,xd0,p,t,0.5,1/12,alpha,beta) 
        return x, xd, xdd

    def LinearAcceleration(self,x0,xd0,p,t,alpha=0,beta=0):    
        # Interface to the linear acceleration Newmark integration
        # Constants 1/2, 1/6       
        x, xd, xdd = self.LinearNewmarkIntegration(x0,xd0,p,t,0.5,1/6,alpha,beta) 
        return x, xd, xdd
     
    def AverageConstantAcceleration(self,x0,xd0,p,t,alpha=0,beta=0):
        # Interface to the average constant acceleration Newmark integration
        # Constants 1/2, 1/4         
        x, xd, xdd = self.LinearNewmarkIntegration(x0,xd0,p,t,0.5,1/4,alpha,beta) 
        return x, xd, xdd  
   
    def LinearNewmarkIntegration(self,x0,xd0,p,t,y,b,alpha,beta):  
        # x0: initial state
        # sd0: derivative of x0
        # p: load generated by createLoad
        # t: time generated by createLoad
        # y: gamma in approximation formulas for the Newmark method
        # b: beta in approximation formulas for the Newmark method
        # alpha&beta: damping coefficients for damping matrix
        
        # General Newmark integration scheme
        
        # M & K
        M = self.M
        K = self.K
        
        # Initialize variables 
        nT      = np.size(t)
        n       = np.size(M,1)
        
        # Damping matrix for structural damping coefficients alpha and beta
        C       = alpha*K+beta*M

        # Displacement, velocity and acceleration
        x       = np.zeros((n,nT))
        xd      = np.zeros((n,nT))
        xdd     = np.zeros((n,nT))

        # Set the initial displacement and velocity
        x[:,0]      = x0
        xd[:,0]     = xd0
        
        
        # Calculate the initial acceleration from the given initial 
        # velocity and displacement. Recall: M xdd+C xd+K x = f
        #++++++++++++ Something is missing here +++++++++++++++
        xdd[:,0] = LAS.inv(M) @ (p[:,0] - K @ x[:,0] - C @ xd[:,0])

        # Retrieve the time step. It is supposed to be constant.
        h = t[1]- t[0]
        
        # A constant time step allows the precalculation of the matrix S
        # Use the LU factorization of the scipy package
        #++++++++++++ Something is missing here +++++++++++++++
        S_LU       = LAS.splu(M + y*h*C + b*h**2*K)


        for i in range(0,np.size(t)-1):
            
            # Prediciton step, index i
            
            # Extract the current displacement, velocity and acceleration
            #++++++++++++ Something is missing here +++++++++++++++
            xCurr = x[:,i]
            xdCurr = xd[:,i]
            xddCurr = xdd[:,i]
            
            # Calculcate the intermediate displacement and velocity
            #++++++++++++ Something is missing here +++++++++++++++
            xCorrP  = xCurr + h*xdCurr + h**2*(0.5-b)*xddCurr
            xdCorrP = xdCurr + (1-y)*h*xddCurr
        
            
            # Solve the LU-decomposet problem to obtain the acceleration
            pP      = p[:,i+1]
            #++++++++++++ Something is missing here +++++++++++++++
            xddP    = S_LU.solve(p[:,i+1] - C@xdCorrP - K@xCorrP)
            
            #Correction step
            
            # Correct the displacement and velocity, store the data in the 
            # next time step i+1
            #++++++++++++ Something is missing here +++++++++++++++
            x[:,i+1] = xCorrP + h**2*b*xddP
            xd[:,i+1] = xdCorrP + y*h*xddP
            xdd[:,i+1] = xddP
          
        
        return x, xd, xdd   

    def createLoad(self,T,LoadDof, Mag=1,Type = 'Rect',deltaT=0.0001):
        # This funtion can be used to generate time series for forcing terms
        
        n       = np.size(self.M,1)
        Load = np.zeros((n , int(T/deltaT)))
        Time = np.linspace(0, T, int(T/deltaT))
 
        if Type == 'Rect':
            Load[LoadDof,10:] = Mag 
            
        elif Type == 'Impulse':
 
             Load[LoadDof,10:15] = Mag 

        
        elif Type == 'Jump':
             print('Not implemented')
           
        elif Type == 'AM':
            freq=25
            xThresh= 0.38
            xAM =np.linspace(-xThresh,xThresh,int(T/deltaT))
            E=0.1+(np.sin(xAM)**2)+0.35*np.exp(-0.5*(xAM/0.4)**2);
            H=np.cos(25*xAM);
            G=np.exp(-0.5*((xAM-0.15)/0.1)**2);
            Load[LoadDof,:] =E*H*(1-0.5*G)+G/6;

            Time = np.linspace(0, T, int(T/deltaT))
        elif Type == 'AMlong':
            freq=25
            xThresh= 3
            xAM =np.linspace(-xThresh,xThresh,int(T/deltaT))
            E=0.1+(np.sin(xAM)**2) +0.35*np.exp(-0.5*(xAM/0.4)**2);
            H=np.cos(25*xAM);
            G=np.exp(-0.5*((xAM-0.15)/0.1)**2);
            Load[LoadDof,:] =E*H*(1-0.5*G)+G/6;

            Time = np.linspace(0, T, int(T/deltaT))        
        else :
            print('Not implemented')
        
        return  Load,Time   
        

class Eigensolver:
    
    def __init__(self,M,K):
        # Eigensolver requires mass and stiffness information 
        self.K = K
        self.M = M
     
    def MaxOmega(self):
        # Calcuate the maximum frequency of the system

        #++++++++++++ Something is missing here +++++++++++++++
        K_inv = LAS.inv(self.K)
        M_inv = LAS.inv(self.M)
        SolverFix_inv = Eigensolver(M_inv,K_inv)
        omega,_ = SolverFix_inv.ScipySparse(n = 1)
#         omega,_ = self.ScipySparse(n = self.K.shape[0]-2)

        maxOmega = 1 / omega[0]
        
        return np.real(maxOmega)
        
    def ScipySparse(self,n=1,shift = 0, epsilon = 1e-4, nit = 100):
        # Solve the structural eigenproblem by using the sparse eigensolver of 
        # the scipy package. 
   
        # Standard Eigenproblem Kx = w2 Mx with which='LM' meaning largest magnitude eigensolutions 
        # and sigma=shift**2 applying the shift factor. If shift=0, the lowest eigenfrequencies 
        # are found.
        
        [eVal, eVec] = LAS.eigs(self.K,n,self.M, which='LM', sigma=shift**2)
        
        # The eigenvalues and eigenvectors must be resorted.
        
        # The numpy array omega2 can be resorted by calling sort()
        eVal.sort()
        # Sort the eigenvectors column-wise. Hint: use the argsort() function
        # of the numpy package to get the required indices.
        eVec = eVec[:,np.argsort(eVal)] 
        
        
        # It is beneficial to mass-normalize the eigenmodes
        eVec = eVec/np.sqrt(np.diag(eVec.T @ self.M @ eVec))

        # Take the real value and avoid NaN by removing negative solutions
        omega2 = np.abs(np.real(eVal))
        # Take the sqrt to obtain the eigenfrequencies
        omega = np.sqrt(omega2)

        # Return the eigenfrequencies ([rad/s]) and eigenmodes 
        return omega, eVec
    
    def PowerIteration(self,n=1,shift = 0, epsilon = 1e-4, nit = 100):
        # Solve the structural eigenproblem by using the power iteration.
        
        # Initialize eigenvectors, -values and the 
        # LU-factorization of the stiffness matrix, regard the shift
        
        eVec = np.zeros([np.size(self.K,1),n])
        eVal = np.zeros([n])
        #++++++++++++ Something is missing here +++++++++++++++
        K = self.K
        M = self.M
        sys_dim = M.shape[0] # dimension of the system
        K_shift = K - (shift**2) * M
        K_lu_shift = LAS.splu(K_shift)

        # For loop: iteration over all to-be-calculated eigenmodes
        # from 0 to n-1
        P = np.eye(sys_dim)
        for it in range(0,n):
            
            # Initialize variables which are required in the while loop
            count = 0
            eval_     = 100
            eval_prev = 10
            
            # There are several possibilities to initialize the current eigenvector
            #++++++++++++ Something is missing here +++++++++++++++
            z = np.reshape(np.ones(sys_dim),(sys_dim,1)) # reshape z to size [sys_dim,1]
            
            # There are several termination criteria, one possibility is to use 
            # the relative change in the calculated eigenfrequency.
            
            while np.abs(eval_ - eval_prev) > epsilon: #++++++++++++ Something is missing here +++++++++++++++
                
                # Calculate the updated eigenvector
                #++++++++++++ Something is missing here +++++++++++++++
                z = K_lu_shift.solve(M @ z)
                
                # Deflation
                # If there are already calculated eigenvectors, deflate them.
                # This can be implemented in one line of code.
                # Alternatively, a for loop iteration over all
                # previously calculated eigenvectors is requried.
                
                if it > 0:
                    #++++++++++++ Something is missing here +++++++++++++++
                    z_i = np.reshape(eVec[:,it-1],(sys_dim,1))
                    P = P @ (np.eye(sys_dim) - z_i @ z_i.T @ M)
                    z = P @ z


                # Mass normalize the current eigenvector
                #++++++++++++ Something is missing here +++++++++++++++
                z = z / np.sqrt(z.T @ M @ z)
                    
                # Set the previous and current values for termination-checking
                eval_prev = eval_
                
                # The currently estimated eigenvalues can be calculated as 
                # ration of mass and stiffness matrix, each projected on the 
                # current eigenvector
                #++++++++++++ Something is missing here +++++++++++++++
                eval_ = z.T @ K_shift @ z
                z = np.reshape(z,(sys_dim,))
                
                # Increase counter, terminate when limit is reached.
                count = count+1
                if count == nit:
                    break 
             
            # Store frequency and vector     
            eVal[it] = eval_ 
            eVec[:,it] = z
            
           
        omega2 = eVal
        eVec = eVec/np.sqrt(np.diag(eVec.T @ self.M @ eVec))
        
        eVec = eVec[:,np.argsort(omega2)]    
        omega2.sort()

        # Resolve shift and problems with square root
        omega2 = omega2+(shift**2)
        omega_nan = np.isnan(np.sqrt(omega2))
        omega = np.sqrt(omega2)
        omega[omega_nan] = omega2[omega_nan]
        return omega,eVec
          
        
        
    def KrylovSubspace(self,n=1,shift=0.0,epsilon= 1e-4, nit = 1000):        
        # Solve the structural eigenproblem by using the krylov subspace method

        # Initialize the shifted stiffness matrix and 
        # the total number of buffer vectors
        K = self.K
        M = self.M
        K_shifted = K - (shift**2) * M
        n_total = M.shape[0]
        # Factorize the stiffness matrix
        K_lu_shift = LAS.splu(K_shifted)
        # counter
        p = 0

        # Define initial guesses for eigenvectors and values
        # Get unitary shapes. Maximize the angles. ??
        
        z = np.reshape(np.ones(n_total),(n_total,1))
        m = M @ z
        gamma = np.sqrt(z.T @ m)
        z = z / gamma
        m = m / gamma
        z_prev = z
 
        Z = z
        M_z = m
        
        
        # Initialize the to-be-calculated eigenmodes
        X = np.zeros((n_total,n))
        T = np.eye(1)
        # Remember: alternative initilization are available


        # Initialize iteration variables
        omegaSqPrev = np.zeros([n]);
        convergence = 0
#         counter = 0 
        
        
        while (convergence == 0) and (p<nit):
            
            # Update the basis of the eigenvectors. Use the factirized, shifted
            # stiffness matrix.
            
#             X = K_lu_shift.solve(M @ z) #++++++++++++ Something is missing here +++++++++++++++

            if p > 0:
                T = np.insert(T, T.shape[0], values = 0, axis = 0)
                T = np.insert(T, T.shape[1], values = 0, axis = 1)
                T[p-1, p] = gamma
                T[p, p-1] = gamma
            z_prev = z
            z = K_lu_shift.solve(m)

            # Project mass and stiffness matrix on the current subspace
            # which results in reduced stiffness and mass matrix
            #++++++++++++ Something is missing here +++++++++++++++
            for i in range(0,p+1):
                a = np.reshape(M_z[:,i],(1,n_total)) @ z
                z = z - a * np.reshape(Z[:,i],(n_total,1))
                if i == p:
                    T[p, p] = a
            m = M @ z
            gamma = np.sqrt(z.T @ m)
            z = z / gamma
            m = m / gamma
            Z = np.column_stack((Z,z))  
            M_z = np.column_stack((M_z,m))
            
            # Solve the reduced eigenproblem. Note: this is not a sparse problem.
            # Use the standard linalg module from ths scipy package.
            if p+1 >= n:
                [omegaSq,Y] = LA.eig(T) #++++++++++++ Something is missing here +++++++++++++++
                omegaSq = np.abs(np.real(omegaSq[0:n]))  
                omegaSq = 1.0 / omegaSq   
                Y = Y[:, 0:n]
                
            
            # Keep the eigenfrequencies and -vectors real.
#                 omegaSq = np.real(omegaSq[n-1])          
            # Project the reduced eigenvector out of the subspace
            # ++++++++++++ Something is missing here +++++++++++++++
#             X = np.real(Z@Xred)    
            
            # Scaling is required. For example, scale the eigenvectors to one.
#             X = X/np.sqrt(np.diag(X.T @ self.M @ X))
            # Convergence and termination criteria
            
            # Calculate the difference of all TARGETED eigenfreqeuncies
            #++++++++++++ Something is missing here +++++++++++++++
                domega = np.absolute(np.sqrt(omegaSq) - np.sqrt(omegaSqPrev)) 
            
            # Increase counter and check convergence
             
                if np.all(domega < epsilon):
                    convergence = 1
            # Update the eigenfrequencies 
                omegaSqPrev = omegaSq
            p = p + 1
         
        # Sort and shift back, avoid problems with square root
        X = Z[:,0:-1]@Y
        omegaSq = np.reshape(omegaSq,(n,))
#         for i in range(0,X.shape[1]):    
#             X[:, i] = X[:, i] / np.linalg.norm(X[:,i])
#         for i in range(0,X.shape[1]):    
#             X[:, i] = X[:, i] / np.amax(X[:,i])
        X = X[:,np.argsort(omegaSq)]    
        omegaSq.sort()

        omega2 = np.real(omegaSq[:n])
        omega2 = omega2+(shift**2)
#         omega_nan = np.isnan(np.sqrt(omega2))
        omega = np.sqrt(omega2)
#         omega[omega_nan] = omega2[omega_nan]
    


        return omega, np.real(X)
        

class HarmonicSolver:
    def __init__(self,M,K):
        # Harmonic solver requires mass and stiffness information
        self.K = K;
        self.M = M;
           
    def DirectHarmonics(self,omegas = np.linspace(0,10,1000)):
        # Direct method: 
        # z = (K - omega^2 * M)^(-1) @ s
        # where s is the load, z is the output(amplitude)
        # For a certain omega, s can be seen as an identity matrix because each input is a unit input
        # in the from of [0,0,...,1,...0]^T, the ith element of the ith s is 1 and the other
        # elements will be 0, all these loads build up an identity matrix.
        # The out put will be (K - omega^2 * M)^(-1)
        
        # Initialize the harmonic response as 3D empty numpy array.
        #++++++++++++ Something is missing here +++++++++++++++    
        K = self.K
        M = self.M
        sys_dim = M.shape[0]
        # The 1st dimension is the input dof while the 2nd the output dof
        # input dof: the node on which the load(amplitude 1) is applied
        # output: contain all the nodes
        HarmonicResponse = np.zeros((sys_dim, sys_dim, omegas.shape[0]))
        for it,omega in enumerate(omegas): 
            # Calculate the entire inverse for the current frequency.
            # The scipy sparse linalg module can be used.
            #++++++++++++ Something is missing here +++++++++++++++
            currHarmonicResponse = LAS.inv(K - (omega**2) * M) 
            
            # Note that a dense and sparse array are different. The todense()
            # function must be applied to assign the harmonic response.
            currHarmonicResponseD = currHarmonicResponse.todense()
            HarmonicResponse[:,:,it]=currHarmonicResponseD
            
        return HarmonicResponse, omegas
    
    def ModalSuperposition(self,omegas = np.linspace(0,10,1000),n = 10):
        
        K = self.K
        M = self.M
        sys_dim = M.shape[0]
        # Modal superposition requries eigenmodes and frequencies 
        # Use an eigensolver of your choise 
        #++++++++++++ Something is missing here +++++++++++++++
        omega_eig, modes = Eigensolver.ScipySparse(self,n=n)
        omega2 = omega_eig**2
        # square the input frequencies
        omegas2 = omegas**2
        # Initialize the harmonic response as 3D empty numpy array.
        HarmonicResponse = np.empty([np.size(self.K,0),np.size(self.K,1),np.size(omegas)])
       
        for it,omega in enumerate(omegas): 
            # Calculate the harmonic response for all modes for the current
            # frequency.
            #++++++++++++ Something is missing here +++++++++++++++
            # loop for superposition, use the first n modes
            currHarmonicResponse = np.empty([np.size(self.K,0),np.size(self.K,1)])
            for i in range(0,n):
                mode = np.real(np.reshape(modes[:,i],(sys_dim,1)))
                m = mode.T @ M @ mode
                currHarmonicResponse += mode @ mode.T / (m * (omega2[i] - omegas2[it]))
            HarmonicResponse[:,:,it]=currHarmonicResponse
        return HarmonicResponse, omegas
        

    def ModeAcceleration(self,omegas = np.linspace(0,10,1000),n = 10):
        # Modal superposition requries eigenmodes and frequencies 
        # Use an eigensolver of your choise       
        # ++++++++++++ Something is missing here +++++++++++++++
        omega_eig, modes = Eigensolver.ScipySparse(self,n=n) 
        omega2 = omega_eig**2
        K = self.K
        M = self.M
        sys_dim = M.shape[0]
        # square the input frequencies
        omegas2 = omegas**2
        # Additionally, the inverse is required. In this case, a factorization has no benefits
        #++++++++++++ Something is missing here +++++++++++++++
        K_inv = LAS.inv(K)
        
        # Initialize the harmonic response as 3D empty numpy array.
        HarmonicResponse = np.empty([np.size(self.K,0),np.size(self.K,1),np.size(omegas)])
        
        for it,omega in enumerate(omegas): 
            # Calculate the harmonic response for all modes for the current
            # frequency. Add the mode acceleration term. 
            #++++++++++++ Something is missing here +++++++++++++++s
            currHarmonicResponse = K_inv
            for i in range(0,n):
                mode = np.real(np.reshape(modes[:,i],(sys_dim,1)))
                m = mode.T @ M @ mode
                currHarmonicResponse += mode @ mode.T / (m * (omega2[i] - omegas2[it])) - mode @ mode.T / (m * omega2[i]) 
            HarmonicResponse[:,:,it]=currHarmonicResponse
            
        
        return HarmonicResponse, omegas
            
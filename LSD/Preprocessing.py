# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 14:02:44 2019

@author: q446161
"""

import numpy as np 
import copy 
from matplotlib import pyplot
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
import scipy.sparse as SP

from amfe.mesh import Mesh
from amfe.component import StructuralComponent
from amfe.material import BeamMaterial, KirchhoffMaterial

class Mapping:
    def __init__(self):
        self.Map  = None

        
    def fromAMfe(self,AMfeMapping,LSDLineGeometry):   
        # Create a map based on a given data from the AMfe toolbox
        N2G = AMfeMapping.nodal2global.values
        NdId = np.arange(np.size(N2G,0)).repeat(6,0).T
        DirId = np.arange(6)
        DirMat= np.repeat(DirId[:,None],np.size(N2G,0),axis=1).T
        Dir = np.reshape(DirMat,np.size(N2G))
        #N2G.sort(axis=1)
        N2G = np.concatenate((N2G[:,3:6],N2G[:,0:3]),axis=1)
        N2G = N2G.reshape((np.size(N2G)))
        
        Coords = LSDLineGeometry.Gridpoints.repeat(6,0)
        
        self.Map  =  np.concatenate((NdId[:,None],N2G[:,None],Dir[:,None],Coords),axis= 1) 
    

    def fromSubMaps(self,SubMaps):  
        # Create a map for an assembly of several substructures
        MapTuple = ()
        for subMap in SubMaps:
            MapTuple = MapTuple+(subMap.Map,)

        self.Map = np.vstack(MapTuple)
        return self
        
            
    def DofNumber(self):
        return np.size(self.Map,0)        
        
    def getByCoords(self,XYZcoords, direction = 'all',tol = 10e-1):
        # Entry ID of the degrees of freedom of a specified point by its xyz coords
        
        foundEntry = np.empty(0,dtype='int')
        for XYZcoord in XYZcoords:
            diff=self.Map[:,3:6] - XYZcoord
            diffnorm = np.linalg.norm(diff,axis=1)
            foundEntry = np.append(foundEntry, np.argwhere(diffnorm < tol)) 

        foundDof = self.Map[foundEntry,1]

        if direction =='all':
            foundDof  = foundDof
        elif direction == 'x':
            foundDof = foundDof[0::6]
        elif direction == 'y':  
            foundDof = foundDof[1::6]
        elif direction == 'z':
            foundDof = foundDof[2::6]
        elif direction == 'rx': 
            foundDof = foundDof[3::6]
        elif direction == 'ry':
            foundDof = foundDof[4::6]
        elif direction == 'rz':     
            foundDof = foundDof[5::6]
            
        allDoF  = np.arange(self.DofNumber())
        unfoundDof = allDoF[np.logical_not(np.isin(allDoF,foundDof))]

        return foundDof.astype(int), unfoundDof.astype(int)
    
    def truncate(self,truncatedDof):
        # If dof of a structure are truncated, the mapping can 
        # be aligned with this function
        
        allDoF  = self.Map[:,1]  
        keptEntries = np.logical_not(np.isin(allDoF,truncatedDof))

        newMap = Mapping()
        newMap.Map    =  self.Map[keptEntries,:]   
        newMap.Map[:,1] = np.arange(np.size(newMap.Map ,0))
        
        return newMap 
    
    def addInteriorDof(self,noInterior):
        # Function to add interior degrees of freedom when craig-bampton 
        # reduction is used
        mapWildcard = np.full((noInterior,6), np.nan)
        mapWildcard[:,1] = np.arange(noInterior)+self.DofNumber() 
                
        newMap = copy.deepcopy(self)
        newMap.Map  =  np.concatenate((newMap.Map,mapWildcard),axis= 0) 
        
        return newMap    
        
    def complementDof(self,truncatedDof):
        # Get all dof which are not specified by truncatedDof
        return self.Map[np.logical_not(np.isin(self.Map[:,1]  ,truncatedDof)),1].astype(int)
    
class MK:
    def __init__(self,M,K,Map,R=None):
        self.M = M
        self.K = K
        self.Map = Map
        self.R = R
            
    def truncate(self,truncatedDof,retainedDof):  
        # Truncate matrices for given degrees of freedom.
        
        # Two steps are required, otherwise the diagonal is extracted
        # Truncate mass matrix columns
        # Truncate mass matrix rows
        # Truncate stiffness matrix columns
        # Truncate stiffness matrix rows
        
        #++++++++++++ Something is missing here +++++++++++++++
        
        M = self.M.toarray()
        M = np.delete(M,truncatedDof,0)
        M = np.delete(M,truncatedDof,1)
        M = SP.csr_matrix(M)
        K = self.K.toarray()
        K = np.delete(K,truncatedDof,0)
        K = np.delete(K,truncatedDof,1)
        K = SP.csr_matrix(K)

            
            
        #++++++++++++ Something is missing here +++++++++++++++
        
        Map = self.Map.truncate(truncatedDof) # Truncate the map
        
        R_ = SP.eye(np.size(truncatedDof)+np.size(retainedDof),dtype=np.int8, format='csr')
        R = R_[:,retainedDof]
        
     
        return MK(M,K,Map,R)

        
class LineGeometry:
    def __init__(self, KeypointFile, LineFile):
        
        self.Keypoints      = np.genfromtxt(KeypointFile, delimiter=',')
        self.Gridpoints     = np.genfromtxt(KeypointFile, delimiter=',')
        LinesArray     = np.genfromtxt(LineFile,dtype= 'int')
        
        # Set Lines as Objects 
        self.Lines = []
        self.LineGroups = []
        self.LineGroupTags = []
              
        for LineData in LinesArray:
            newLine = Line(LineData,self.Keypoints[LineData[1],:],self.Keypoints[LineData[2],:])
            self.Lines.append(newLine)
            if newLine.GRU in self.LineGroupTags:
                for i, LineGroupTag in enumerate(self.LineGroupTags):
                    if newLine.GRU == LineGroupTag:
                        self.LineGroups[i].append(newLine) 
            else:
                newLineGroup = [newLine]
                self.LineGroupTags.append(newLine.GRU)
                self.LineGroups.append(newLineGroup)
            
        self.LineNodes      = 2
    
    def discretizeLineGeometry(self):   
        self.LineElements = [];
        lineCounter = 0;
        
        for LineGroup, NodeNumber in zip(self.LineGroups, self.LineNodes):
            for Line_ in LineGroup:
                
                noNewP = NodeNumber-2;
                IntP = np.linspace(Line_.KP0,Line_.KP1,NodeNumber);
                NewP = IntP[1:-1,:]; 
                pointOffset = np.size(self.Gridpoints,0)
                self.Gridpoints = np.vstack((self.Gridpoints,NewP)); # Unitl here, it's correct
                newPoitsId = np.arange(pointOffset,pointOffset+noNewP,dtype='int');
                NewLinesChained = np.hstack((Line_.KP0id, newPoitsId, Line_.KP1id));
                
                for x,y in zip(NewLinesChained[0:-1], NewLinesChained[1:]):
                    self.LineElements.append(Line([lineCounter,x,y, Line_.TRI,Line_.GRU],self.Gridpoints[x,:],self.Gridpoints[y,:]))
                    lineCounter = lineCounter+1 
                    #self.LineElements = np.vstack((self.LineElements,[x, y]));
    
    def assembleLineElements(self,D = 0.2, d= 0.195,E = 2.1e11,nu = 0.3,rho = 7500 ):# The return MKfree is an MK class
        
        G = E/(2*(1+nu))
        A   = np.pi*(D**2 - d**2)/4
        I_y = np.pi*(D**4 - d**4)/64
        I_z = np.pi*(D**4 - d**4)/64
        J_x = np.pi*(D**4 - d**4)/32
        J_x = np.pi*(D**4 - d**4)/64
        
        X3 = np.array([1E21, 1E21, 1E21]) 
        X3 = np.array([30, 15,15]) 

        mesh = Mesh(3)
        
        for i, Keypoint in enumerate(self.Gridpoints):
            mesh.add_node(Keypoint, i)
            
        s = 'straight_line'
        for i,LineElement in enumerate(self.LineElements):
            mesh.add_element(s, np.array([LineElement.KP0id,LineElement.KP1id]), i)
            
        component = StructuralComponent(mesh)
        matBeam = BeamMaterial(E, G, rho, A, I_y, I_z, J_x, X3)

        component.assign_material(matBeam, s, 'S', 'shape')
        mapping = component._mapping
        
        t0 = 0
        x0 = np.zeros(mapping.no_of_dofs)
        dx0 = x0.copy()
        Mfree = component.M(x0, dx0, t0)
        Kfree = component.K(x0, dx0, t0)
        
        MAPfree = Mapping()
        MAPfree.fromAMfe(mapping,self)
        
        OrigDofOrder = MAPfree.Map[:,1]
        OrigDofOrder_ = np.zeros(MAPfree.DofNumber(),dtype='int')



        for dofId in np.arange(MAPfree.DofNumber()):
            OrigDofOrder_[dofId] =np.argwhere(OrigDofOrder == dofId)
            
        Kfree_ = Kfree[OrigDofOrder,:]
        Kfree = Kfree_[:,OrigDofOrder]
        Mfree_  = Mfree[OrigDofOrder,:]
        Mfree = Mfree_ [:,OrigDofOrder]
        
        
        MAPfree.DofId    = np.arange(MAPfree.DofNumber(),dtype='int')
        MAPfree.Map[:,1] = np.arange(MAPfree.DofNumber(),dtype='int')

        MKfree = MK(Mfree ,Kfree ,MAPfree)
        return MKfree
              

    def plotKeypoints(self, axis = None):
        
        if axis == None:
            fig = pyplot.figure()
            ax  = Axes3D(fig) 
        else:
            ax = axis
        ax.scatter(self.Keypoints[:,0], self.Keypoints[:,1], self.Keypoints[:,2])
        
    def plotLines(self, axis = None, color = 'black'):
        
        if axis == None:
            fig = pyplot.figure()
            ax  = Axes3D(fig) 
        else:
            ax = axis
            
        for Line in self.Lines:
            x_ = [Line.KP0[0], Line.KP1[0]]
            y_ = [Line.KP0[1], Line.KP1[1]]
            z_ = [Line.KP0[2], Line.KP1[2]]
            ax.plot(x_,y_,z_,color , marker = 'o')           

       
    def plotLineElements(self, axis = None,color = 'black'):
        colors = ['olivedrab','crimson','black','steelblue','darkorange']
        if axis == None:
            fig = pyplot.figure()
            ax  = Axes3D(fig) 
        else:
            ax = axis
            
        for Line in self.LineElements:
            x_ = [Line.KP0[0], Line.KP1[0]]
            y_ = [Line.KP0[1], Line.KP1[1]]
            z_ = [Line.KP0[2], Line.KP1[2]]
            ax.plot(x_,y_,z_, color, marker = 'o')           
      
    
    def plotDisplacement(self,dX,  axis = None, amp = 10E2):
        dx = dX[0::6]*amp
        dy = dX[1::6]*amp
        dz = dX[2::6]*amp
        if axis == None:
            fig = pyplot.figure()
            ax  = Axes3D(fig) 
        else:
            ax = axis
        ax.scatter(self.Gridpoints[:,0], self.Gridpoints[:,1], self.Gridpoints[:,2])
        ax.scatter(self.Gridpoints[:,0]+dx, self.Gridpoints[:,1]+dy, self.Gridpoints[:,2]+dz,c='black')
        
    def plotDisplacementLine(self,dX, axis = None, amp = 10E2):
        dx = dX[0::6]*amp
        dy = dX[1::6]*amp
        dz = dX[2::6]*amp
        
        PlotPoints = self.Gridpoints + np.vstack((dx,dy,dz)).T
                
        if axis == None:
            fig = pyplot.figure()
            ax  = Axes3D(fig) 
        else:
            ax = axis
        for Line in self.LineElements:
            x = [Line.KP0[0], Line.KP1[0]]
            y = [Line.KP0[1], Line.KP1[1]]
            z = [Line.KP0[2], Line.KP1[2]]
            dx= [PlotPoints[Line.KP0id ,0],PlotPoints[Line.KP1id ,0]]
            dy= [PlotPoints[Line.KP0id ,1],PlotPoints[Line.KP1id ,1]]
            dz= [PlotPoints[Line.KP0id ,2],PlotPoints[Line.KP1id ,2]]
            
            ax.plot(x,y,z, color = 'black', marker = 'o')   
            ax.plot(dx,dy,dz, color = 'darkred', marker = 'o')  
            

    def animateDisplacementLine(self,dX, axis = None, amp = 10E2):
        dx = dX[0::6]*amp
        dy = dX[1::6]*amp
        dz = dX[2::6]*amp
        fig = pyplot.figure()

        animLines =[]
        PlotPoints = self.Gridpoints + np.vstack((dx,dy,dz)).T
                
        if axis == None:
            fig = pyplot.figure()
            ax  = Axes3D(fig) 
        else:
            ax = axis
        for Line in self.LineElements:
            x = [Line.KP0[0], Line.KP1[0]]
            y = [Line.KP0[1], Line.KP1[1]]
            z = [Line.KP0[2], Line.KP1[2]]
            dx= [PlotPoints[Line.KP0id ,0],PlotPoints[Line.KP1id ,0]]
            dy= [PlotPoints[Line.KP0id ,1], PlotPoints[Line.KP1id ,1]]
            dz= [PlotPoints[Line.KP0id ,2], PlotPoints[Line.KP1id ,2]]
            
            ax.plot(x,y,z, color = 'black', marker = 'o')   
            #animLine = self.ax.plot(dx,dy,dz, color = 'darkred', marker = 'o')   
            animLine, = ax.plot([],[],[], color = 'darkred', marker = 'o')   
            animLines.append(animLine)
        
        animLines =  animLines
           
        def init():
            for animLine in animLines:
                #animLine.set_data([],[],[])
                animLine.set_data([],[])
                animLine.set_3d_properties([])
            return animLines
            
        def animate(i):
            for animLine,Line in zip(animLines,self.LineElements):
                dx= [PlotPoints[Line.KP0id , 0],PlotPoints[Line.KP1id ,0]]
                dy= [PlotPoints[Line.KP0id ,1], PlotPoints[Line.KP1id ,1]]
                dz= [PlotPoints[Line.KP0id ,2], PlotPoints[Line.KP1id ,2]]
                animLine.set_data(i*dx,i*dy)
                animLine.set_3d_properties(i*dz)
                print(i)
            return animLines
        animation.FuncAnimation(fig, animate, init_func=init)#,frames=100, interval=100, blit=True)
        
class Line: 
    def __init__(self,LineArray,KP0, KP1):
        self.ID     = LineArray[0]
        self.KP0id  = LineArray[1] 
        self.KP1id  = LineArray[2] 
        self.TRI    = LineArray[3]         
        self.GRU    = LineArray[4] 
        self.KP0    = KP0
        self.KP1    = KP1

        


 
    



    







    
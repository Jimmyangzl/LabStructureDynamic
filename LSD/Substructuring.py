import numpy as np

from scipy.sparse import linalg  as LAS
import scipy.sparse as SP
import copy 
import LSDstudent.Preprocessing as Prep
import LSDstudent.Solver as Sol
        

class Assembly:
    def __init__(self, Substructures):
        self.Substructures =  Substructures
        
        
    def assemble(self):
        
        Ktuple = ()
        Mtuple = ()
        MapTuple = ()
        
        for Substruture in self.Substructures:
            Ktuple = Ktuple+(Substruture.K,)
            Mtuple = Mtuple+(Substruture.M,)
            MapTuple = MapTuple+(Substruture.Map.Map,)
            
        self.Kb = SP.block_diag(Ktuple)
        self.Mb = SP.block_diag(Mtuple)
        
        self.Map = Prep.Mapping()
        self.Map.Map = np.vstack(MapTuple)
        self.Map.Map[:,1] = np.arange(self.Map.DofNumber())

        self.L,self.minDof = self.getLfromMap()
        
        self.M = self.L.T @self.Mb @self.L
        self.K = self.L.T @self.Kb @self.L
        
        
        self.MapMin = Prep.Mapping()
        self.MapMin.Map = self.Map.Map[self.minDof,:]
        self.MapMin.Map[:,1] = np.arange(self.MapMin.DofNumber())

        
        return Prep.MK(self.M,self.K,self.MapMin)
        

    def getInterfaceCoords(self):
        smallMap = self.Map.Map[0::6,3:6]
        unq, cnt = np.unique(smallMap, axis=0, return_counts=True)
        return unq[cnt>1,:]

    
    def getLfromMap(self):
        DirCoords = self.Map.Map[:,2:6]
        
        foundFlag = np.full(np.size(DirCoords,0), False)
        col = 0
        CRE = ()
        Row = []
        Col = []
        Data = []
        minDof = []

        for i,dirCoord1 in enumerate(DirCoords):
             
            if not foundFlag[i]:
                CRE = CRE+(1,(i,col))
                Data.append(1)
                Col.append(col)
                Row.append(i)
                minDof.append(i)
                if i < np.size(DirCoords,0)-1:
                    
                    for j in range(i+1, np.size(DirCoords,0)):
                        dirCoord2 = DirCoords[j,:]
                        if np.linalg.norm(dirCoord1[1:]-dirCoord2[1:]) < 10e-3 and dirCoord1[0].astype(int) == dirCoord2[0].astype(int):
                            foundFlag[j] = True
                            CRE = CRE+(1,(j,col))  
                            Data.append(1)
                            Col.append(col)
                            Row.append(j)
                    
                col+=1
             
             
        return SP.csc_matrix((Data,(Row,Col))), minDof
            
            


class Reduction:
    def __init__(self,M,K,Map):
        self.M = M
        self.K = K
        self.Map = Map
        
    def Guyan(self, masterDof,slaveDof):
        newOrder = np.concatenate((masterDof, slaveDof),axis=0)

        M_  = self.M[newOrder,:]
        M   = M_[:,newOrder]
        
        K_  = self.K[newOrder,:]
        K   =  K_[:,newOrder]        
        

        # Calculate the static condensation matrix S, several steps are missing. 



        S = SP.eye(np.size(slaveDof),np.size(masterDof)) # This is only a wildcard
        R = np.concatenate((np.eye(np.size(masterDof)),np.asarray(S.todense())),axis=0)
        
        Kr = R.T @K@R
        Mr = R.T @M@R
        
        
        MapR = Prep.Mapping
        MapR.Map = self.Map.Map[masterDof,:]         
        R[newOrder,:] = R
    
        return Mr,Kr,R, MapR
        
    def CraigBampton(self, masterDof, slaveDof, interiorDof = 10):
        Ms = self.M[slaveDof,:]
        Ks = self.K[slaveDof,:]            
    
        Mss = Ms[:,slaveDof]
        
        Ksm = Ks[:,masterDof]
        Kss = Ks[:,slaveDof]        
        
        
        newOrder = np.concatenate((masterDof, slaveDof),axis=0)

        M_  = self.M[newOrder,:]
        M   = M_[:,newOrder]
        
        K_  = self.K[newOrder,:]
        K   =  K_[:,newOrder]        
        

        S = np.zeros(np.size(slaveDof),np.size(slaveDof))
        _, R22 = np.zeros(np.size(slaveDof),np.size(slaveDof))

        # Your task: build the reduction matrix!
#        R1 =...
#        R21 = ...
#        R2  = ...
#        R =  np.concatenate((R1,R2),axis=1)
        
        R = np.eye(np.size(masterDof)+interiorDof,np.size(slaveDof)) # This is only a wildcard
        
        Kr = R.T @K@R
        Mr = R.T @M@R
        
        MapR = self.Map.truncate(slaveDof)
        MapR = MapR.addInteriorDof(interiorDof)  
        R[newOrder,:] = R
        
        return Mr,Kr,R, MapR
        
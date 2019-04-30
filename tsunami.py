#
# PYTHON for FEM DUMMIES 18-19
# Projet "tsunami"
#
# Canevas de départ
#  Vincent Legat
#
# -------------------------------------------------------------------------
# 

import numpy as np

# -------------------------------------------------------------------------

def readMesh(fileName) :
  with open(fileName,"r") as f :
    nNode = int(f.readline().split()[3]) #nbre de noeuds
    xyz   = np.array(list(list(float(w) for w in f.readline().split()[2:]) for i in range(nNode)))
    nElem = int(f.readline().split()[3]) #nbre de triangles
    elem  = np.array(list(list(int(w)   for w in f.readline().split()[2:]) for i in range(nElem))) #triangles
  X = xyz[:,0] #coordonnées x de chaque triangle
  Y = xyz[:,1] #coordonnées y de chaque triangle
  H = xyz[:,2] #bathymétrie 
  return [nNode,X,Y,H,nElem,elem]

# -------------------------------------------------------------------------

def findEdges(fileName) :
    nElem = readMesh(fileName)[4]
    nEdges = nElem*3 
    elem = mesh[5]
    nBoudary = 0
    edges = [[0 for i in range(4)] for i in range(nEdges)]
    for i in range (nElem) :
      for j in range(3) :
        id = i*3 + j
        edges[id][0] = elem[i][j]
        edges[id][1] = elem[i][(j+1)%3]
        edges[id][2] = i
        edges[id][3] = -1
    
    edges.sort(key = lambda item : -(min(item[0:2])*self.nEdges)-max(item[0:2])) 
    index = 0
    
    for i in range(nEdges) :
        if (edges[i][0:2] != edges[i-1][1::-1]) :
         edges[index] = edges[i]
         index += 1
        else:
         edges[index-1][3] = edges[i][2]
    del edges[index:]
    edges.sort(key = lambda item : item[3])
    nBoundary = 2*index - nEdges
    nEdges = index
         
    return [nEdges,nBoundary,edges]

# -------------------------------------------------------------------------

def readResult(fileBaseName,iter,nElem) :
  fileName = fileBaseName % iter
  with open(fileName,"r") as f :
    nSize = int(f.readline().split()[3])
    if (nElem != nSize) :
      print(" ==== Error : incoherent sizes : %d != %d" % (nElem,nSize))     
    E = np.array(list(list(float(w) for w in f.readline().split()[2:5]) for i in range(nElem)))
    print(" === iteration %6d : reading %s ===" % (iter,fileName))
  return E

# -------------------------------------------------------------------------

def writeResult(fileBaseName,iter,E) :
  fileName = fileBaseName % iter
  nElem = E.shape[0]  
  with open(fileName,"w") as f :
    f.write("Number of elements %d\n" % nElem)
    for i in range(nElem):
      f.write("%6d : %14.7e %14.7e %14.7e\n" % (i,*E[i,:]))
    print(" === iteration %6d : writing %s ===" % (iter,fileName))
 
# -------------------------------------------------------------------------

def initialConditionOkada(x,y) :
  R = 6371220;
  x3d = 4*R*R*x / (4*R*R + x*x + y*y);
  y3d = 4*R*R*y / (4*R*R + x*x + y*y);
  z3d = R*(4*R*R - x*x - y*y) / (4*R*R + x*x + y*y);
  lat = np.arcsin(z3d/R)*180/np.pi;
  lon = np.arctan2(y3d,x3d)*180/np.pi;
  lonMin = 142;
  lonMax = 143.75;
  latMin = 35.9;
  latMax = 39.5;
  olon = (lonMin+lonMax)/2;
  olat = (latMin+latMax)/2;
  angle = -12.95*np.pi/180; 
  lon2 = olon + (lon-olon)*np.cos(angle) + (lat-olat)*np.sin(angle);
  lat2 = olat - (lon-olon)*np.sin(angle) + (lat-olat)*np.cos(angle);
  return np.all([lon2 <= lonMax,lon2 >= lonMin,lat2 >= latMin,lat2 <= latMax],axis= 0).astype(int)
#retourne une matrice de dim nbrTriangles X 3
 
# -------------------------------------------------------------------------

def mapEdge(theEdges,iEdge) :
  myEdge = theEdges.edges[iEdge]
  elementLeft  = myEdge[2]
  nodesLeft    = elem[elementLeft]
  mapEdgeLeft  = [3*elementLeft + np.nonzero(nodesLeft == myEdge[j])[0][0]  for j in range(2)]
  mapLeft      = [3*elementLeft + j                                         for j in range(3)]
  elementRight = myEdge[3]
  if (elementRight != -1) :
    nodesRight   = theMesh.elem[elementRight]
    mapEdgeRight = [3*elementRight + np.nonzero(nodesRight == myEdge[j])[0][0] for j in range(2)]
    mapRight     = [3*elementRight + j                                         for j in range(3)]
  else :
    mapEdgeRight = []
    mapRight     = []
  return [mapEdgeLeft,mapLeft,mapEdgeRight,mapRight]
#faut faire en sorte d avoir acces a edges

# -------------------------------------------------------------------------

def compute(theMeshFile,theResultFiles,U,V,E,dt,nIter,nSave):
  #U V et E de longueur nElem
  [nNode,X,Y,H,nElem,elem] = readMesh(theMeshFile)
  xsi3 = np.array([0.166666666666667,0.666666666666667,0.166666666666667])
  eta3 = np.array([0.166666666666667,0.166666666666667,0.666666666666667])
  weight3 = np.array([0.166666666666667,0.166666666666667,0.166666666666667])
  R = 6371220;
  gamma = 10**-7
  g = 9.81
  omega = 2*np.pi/86400
  
  for i in range nIter:
      u = np.copy(U)
      U = np.zeros((nElem,3))
      v = np.copy(V)
      V = np.zeros((nElem,3))
      e = np.copy(E)
      E = np.zeros((nElem,3))
      for iElem in range nElem:
          #boucle pour initialiser tous les elements necessaires aux calcul
      for i in range(nElem):
          #boucle pour les integrales sur les triangles
          nodes = elem[i]    
          
          x = X[nodes]
          y = Y[nodes]
          h = H[nodes]
          for j in range(3):
               stereo = (4*R*R+x[j]*x[j]+y[j]*y[j])/(4*R*R)
               #E[i] += weight3[j]*(u[i]*dphidx[j] + v[i]*dphidy[j])*h[j]*jac*stereo #terme 1 de E
               #E[i] += weight3[j]*h[j]/(R*R)*(x[j]*u[i]+y[j]*V[i]) #demander comment calculer le phi
          if(i%nSave==0):
              #write dans un fichier
          
          

  #calcul de E (l elevation de l eau pour chaque triangle)
  #besoin de dphidx et dphidy,h,u,v,x,y

  


           
  #for i in range(nEdges):
      #[mapEdgeLeft,mapLeft,mapEdgeRight,mapRight] = mapEdge(theEdges,i)
      #boucle pour les integrales sur les segments
      
      
  #calcul de U (la vitesse horizontale de l eau pour chaque triangle)
  #besoin de phi,f,v,gamma,u,dphidx,eta (elevation du niveau de la mer)du triangle et d a cote,nx,x et y
  #nx et ny
               
  #calcul de V (la vitesse verticale de l eau pour chaque triangle)
  
  return [U,V,E]
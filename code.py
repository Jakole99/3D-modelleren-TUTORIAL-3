import numpy as np
import random as rd
import bpy
import bmesh
import itertools 
import math

def random():
    randomNumber = rd.uniform(1,4)
    return round(randomNumber, 1)


def deleteAll():
    bpy.ops.object.select_all(action='SELECT')
    bpy.ops.object.delete()
    
    for material in bpy.data.materials:
        bpy.data.materials.remove(material)
    
    for meshes in bpy.data.meshes:
        bpy.data.meshes.remove(meshes)
        
       
#1 Compute and load a control polygon.
def controlMesh(n, length):

    range = np.linspace(0, length, num=n)
    inner = np.linspace(1, n-2, num=n-2)

    controlPolygon = [(x,y,random()) if np.isin(i, inner) and np.isin(j, inner) 
    else (x,y,0.0) for i, x in enumerate(range) for j, y in enumerate(range)]

    return controlPolygon
        
        

#2 Display a planar mesh in blender.

#I wanted to make a difference between showing the control mesh and the final mesh,
#so that you could see the difference easier.
#I do this by drawing the faces of the final mesh and not of the control mesh.

#I wanted to give the option to the Grader, to see the final mesh with or without the faces.
def showMesh(A, n, withFaces):
    if(withFaces):
        showMeshFaces(A, n)
    else:
        showMeshEdges(A, n)
        

#This one only shows the lines of the edges between vertices.
def showMeshEdges(A, n):
    
    mesh = bpy.data.meshes.new("Surface")
    obj = bpy.data.objects.new("Surface", mesh)

    scene = bpy.context.scene
    scene.collection.objects.link(obj)

    vertices = A
    
    print("vertices")
    print(vertices)

    #Combine every vertex with its x neighbour and y neighbour
    numbers = list(range(0,n*n))
    edges = [(i,j) for (i,j) in itertools.combinations(numbers, 2) if j - i == n or (i,j)[1] == (i,j)[0] + 1 and (i,j)[1]%n != 0] 
    
    mesh.from_pydata(vertices, edges, [])
    mesh.update()
    
    
#This one shows the object with the surface filled.
def showMeshFaces(A, n):
    mesh = bpy.data.meshes.new("Surface")

    bm = bmesh.new()

    vertices = A
    
    print("vertices")
    print(vertices)
    
    #Assign indices to all vertices according to the nxn, so that i can easily 
    #specify with vertex I want.
    numbers = list(range(0,n))
    indices = [(i,j) for i, x in enumerate(numbers) for j, y in enumerate(numbers)] 
    combinations = [val for val in zip(vertices,indices)]
    
    #This is how I get all the different quads
    quads = []
    for x in range(n-1):
        for y in range(n-1):
            quad = [vert for (vert,(i,j)) in combinations if i >= x and i <= x + 1 and j >= y and j <= y+1]
            quads.append(quad)
            
    for quad in quads:
        #rearange the order, because otherwise the face isn't connected correclty.
        tmp = quad[2]
        quad[2] = quad[3]
        quad[3] = tmp
        
        quad = tuple(bm.verts.new(vert) for vert in quad)
        bm.faces.new(quad)
        
    bm.to_mesh(mesh)
    bm.free()

    obj = bpy.data.objects.new("Surface", mesh)

    scene = bpy.context.scene
    scene.collection.objects.link(obj)


def showMeshes(S,n, withFaces):
    if (any(isinstance(i, list) for i in S)): #check if it is a list of meshes, https://blog.finxter.com/how-to-check-if-a-list-is-nested-in-python/
        for mesh in S:
            showMesh(mesh, n, withFaces)
    else:
        #only controlMesh
        showMesh(S, n, withFaces)



#3 Approximate a surface on the polygon using the De Casteljau algorithm.
def interpolation(A):
    res = []
    A = np.array(A)
    for i in range(len(A)-1):
        res.append(A[i])
        res.append((A[i]+A[i+1])/2)
        
    res.append(A[-1])
    
    res = [tuple(arr) for arr in res]
    return list(res)


def casteljauReduction(A, res):
    mid = math.floor(len(res)/2)
# The reduction logic
#    A                                res
#    [0, 3, 3, 0]               ->    [., ., ., ., ., ., .]                 Start
#    [0, 1.5, 3, 3, 3, 1.5, 0]  ->    [., ., ., ., ., ., .]                 add (sum neighbours / 2) between each element of A
#    [1.5, 3, 1.5]              ->    [0, ., ., ., ., ., 0]                 remove and add the outerpoints from A to res in the correct index and remove each second element from A
#    [1.5, 2.25, 3, 2.25, 1.5]  ->    [0, ., ., ., ., ., 0]                 add (sum neighbours / 2) between each element of A
#    [2.25, 2.25]               ->    [0, 1.5, ., ., ., 1.5, 0]             Remove and add the outerpoints from A to res in the correct index and remove each second element from A
#    [2.25, 2.25]               ->    [0, 1.5, ., ., ., 1.5, 0]             Notice there are only 2 elements left in A
#    [2.25, 2.25, 2.25]         ->    [0, 1.5, 2.25, 2.25, 2.25, 1.5, 0]    End
    if (len(A) > 2):
        #add (sum neighbours / 2) between each element of A
        A = interpolation(A)
        
        B = A[1:-1]
        B = B[::2]
        
        index = len(B)
        
        #Add the outerpoints of A to res in the correct index
        res[mid - index] = A[0]
        res[mid + index] = A[-1]
        
        #Remove outerpoints from A
        A = A[1:-1]
        
        #Remove every second element from A
        A = A[::2]
        
        #Go on
        casteljauReduction(A, res)
        
    else:
        #add (sum neighbours / 2) between each element of A
        A = interpolation(A)

        #Add the rest of A to res.
        res[mid-1] = A[0]
        res[mid]   = A[1]
        res[mid+1] = A[-1]
        
    res = [tuple(list) for list in res]
    return res
    
    
def deCasteljau(A, n):
    
    #x
    numpyA = np.array(A) 
    splitA = np.split(numpyA, n) 
    rowsA = [splitA.tolist() for splitA in splitA]
    
    nextSize = 2 * n - 1
    emptyNextSize = np.empty((nextSize,3))
    
    A = [casteljauReduction(rowsA, emptyNextSize) for rowsA in rowsA]
    Ax = []
    for listX in A:
        Ax.extend(listX) #I just want a single list, not a list of lists, so that I can repeat for y.

    #y
    columnsA = [Ax[i::nextSize] for i in range(nextSize)]
    emptyNextSize = np.empty((nextSize,3))
    
    A = [casteljauReduction(columnsA, emptyNextSize) for columnsA in columnsA]
    
    refinedControlPoints = [listY for column in zip(*A) for listY in column]
    

    #split into 4 parts
    numbers = list(range(0,nextSize))
    indices = [(i,j) for i, x in enumerate(numbers) for j, y in enumerate(numbers)] 
    combinations = [val for val in zip(refinedControlPoints,indices)]
    
    Q1 = [coord for (coord,(i,j)) in combinations if i >= 0 and i <= n-1 and j >= 0 and j <= n-1]
    Q2 = [coord for (coord,(i,j)) in combinations if i >= 0 and i <= n-1 and j >= n-1 and j <= nextSize]
    Q3 = [coord for (coord,(i,j)) in combinations if i >= n-1 and i <= nextSize and j >= 0 and j <= n-1]
    Q4 = [coord for (coord,(i,j)) in combinations if i >= n-1 and i <= nextSize and j >= n-1 and j <= nextSize]
    res = [Q1] + [Q2] + [Q3] + [Q4]

    return res
    
    
def iteratedDeCasteljau(A, n, s):
    res = []
    if (s > 0):
        if (any(isinstance(i, list) for i in A)): #https://blog.finxter.com/how-to-check-if-a-list-is-nested-in-python/
            for sublist in A:
                res += deCasteljau(sublist, n)
        else:
            res += (deCasteljau(A, n))
            
        return iteratedDeCasteljau(res, n, s-1)
    
    else: 
        return A
        


#4.  Implement an efficient intersection testl
def intersectTest(A, p1, p2, centre, r):
    
    direction = p2 - p1
    
    a =  np.sum((direction)**2) 
    b = 2 * np.sum((direction) * (p1 - centre))
    c = np.sum((p1 - centre)**2) - r**2 
        
    discriminant = b**2 - 4*a*c
    return discriminant >= 0


def lineIntersect(A, n, p1, p2, e):
    stack = [A]
    while stack:
        #We test the first element in the stack
        tupleMesh = stack.pop()
        
        #setup to check the bounding sphere
        npMesh = np.array(tupleMesh)
        max = np.max(npMesh, axis=0)
        min = np.min(npMesh, axis=0)
        dia = max-min
        centre = dia/2 + min
        r = np.linalg.norm(dia/2)
        
        #check if it hits the bounding sphere
        if (intersectTest(npMesh, p1, p2, centre, r)):
            
            #if so, and the box is bigger then e, we can subdivide
            if (r > e):
                nextMeshes = deCasteljau(tupleMesh, n)
                
            #else, its close enough, so we can return True
            else:
                return True
        else:
            nextMeshes = []
            
        stack.extend(nextMeshes)
        
    #Stack is empty, nothing has been intersected
    return False


def getPointOnLine(p1,p2,t):
    return p1 + (p2-p1)*t
    
    
def drawLine(p1, p2):
    
    mesh = bpy.data.meshes.new("Line")
    obj = bpy.data.objects.new("Line", mesh)
    
    scene = bpy.context.scene
    scene.collection.objects.link(obj)
    
    bm = bmesh.new()
    
    v1 = bm.verts.new(getPointOnLine(p1,p2,8))
    v2 = bm.verts.new(getPointOnLine(p1,p2,-8))

    bm.edges.new((v1, v2))

    bm.to_mesh(mesh)
    bm.free()



#Main
def main():
    
    deleteAll()

    n = 6
    length = 4
    A = controlMesh(n,length)

    #test controlMesh(A, n)
    showMesh(A, n, withFaces=False)

    #test deCasteljau(A, n)
    #showMeshes(deCasteljau(A, n), n)

    #test IteratedDeCasteljau(A, n, s)
    s = 4
    showMeshes(iteratedDeCasteljau(A, n, s), n, withFaces=True)

    # True test lineIntersect(A, n, p1, p2, e)
    p1 = np.array((1,3,1))
    p2 = np.array((2,2,-2))
    e = 0.01
    print("True test")
    print("-----")
    print(lineIntersect(A, n, p1, p2, e))
    print("-----")

    drawLine(p1, p2)
    
    #False test lineIntersect(A, n, p1, p2, e)
    p3 = np.array((1,-3,1))
    p4 = np.array((2,2,-2))
    print("False test")
    print("-----")
    print(lineIntersect(A, n, p3, p4, e))
    print("-----")

    drawLine(p3, p4)

main()







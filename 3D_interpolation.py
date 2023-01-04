from shapely.geometry import Point, Polygon
import scipy
import numpy as np
import pandas as pd
from stl import mesh
from scipy.spatial import ConvexHull, convex_hull_plot_2d


def calcul_z_moy(data):
    sum = 0
    nb_pts = len(data)
    for i in range(nb_pts):
        sum += data[i][2]
    return sum/nb_pts

def calcul_SST(data):
    nb_pts = len(data)
    z_moy = calcul_z_moy(data)
    sum = 0
    for i in range(nb_pts):
        z = data[i][2]
        sum += (z- z_moy)**2 

    return sum

def calcul_SSR(data,A,C):
    z_moy = calcul_z_moy(data)
    n = len(data)
    vecteur = []
    for index in range(n):
        vecteur.append(data[index][2])
    Z = np.dot(A,C)
    sum = 0
    for index in range(n):
        sum += (z_moy-Z[index])**2
    return sum

def calcul_R2(data,A,C):
    SST = calcul_SST(data)
    SSR = calcul_SSR(data,A,C)
    return SSR/SST

def gaussian_curvature(data,X,Y,C,degree):
    f1x = calcul_f1x(X,Y,C,degree)
    f1y = calcul_f1y(X,Y,C,degree)
    f2xy = calcul_f2xy(X,Y,C,degree)
    f2xx = calcul_f2xx(X,Y,C,degree)
    f2yy = calcul_f2yy(X,Y,C,degree)
    hull = ConvexHull(data[:,0:2])
    n = len(X)
    Z = []
    for i in range(n):
        m = len(X[i])
        list = []
        for j in range(m):
            list.append(0)
        Z.append(list)
    for i in range(n):
        m = len(X[i])
        for j in range(m):
            Z[i][j]  = (f2xx[i][j]*f2yy[i][j]-(f2xy[i][j]**2))/((1+f1x[i][j]**2+f1y[i][j]**2)**2)*(10**6)
    return Z

def calcul_f2xx(X,Y,C,degree):
    matrix = []
    i = 0
    for index in range(degree+1):
        list = []
        for index2 in range(degree-index+1):
            list.append(C[i])
            i += 1
        matrix.append(list)
    f2xx = []
    for index in range(len(X)):
        list = []
        for index2 in range(len(X[index])):
            list.append(0)
        f2xx.append(list)
    n = len(X)
    for i in range(n):
        m = len(X[i])
        for j in range(m):
            sum = 0 
            for index in range(2,degree+1):
                for index2 in range(degree-index+1):
                    sum += (index)*(index-1)*matrix[index][index2]*(X[i][j]**(index-2))*(Y[i][j]**index2)  
            f2xx[i][j] = sum
    return(f2xx)

def calcul_f2yy(X,Y,C,degree):
    matrix = []
    i = 0
    for index in range(degree+1):
        list = []
        for index2 in range(degree-index+1):
            list.append(C[i])
            i += 1
        matrix.append(list)
    f2yy = []
    for index in range(len(X)):
        list = []
        for index2 in range(len(X[index])):
            list.append(0)
        f2yy.append(list)
    n = len(X)
    for i in range(n):
        m = len(X[i])
        for j in range(m):
            sum = 0 
            for index in range(0,degree+1):
                for index2 in range(2,degree-index+1):
                    sum += (index2)*(index2-1)*matrix[index][index2]*(X[i][j]**(index))*(Y[i][j]**(index2-2))  
            f2yy[i][j] = sum
    return(f2yy)

def calcul_f2xy(X,Y,C,degree):
    matrix = []
    i = 0
    for index in range(degree+1):
        list = []
        for index2 in range(degree-index+1):
            list.append(C[i])
            i += 1
        matrix.append(list)
    f2xy = []
    for index in range(len(X)):
        list = []
        for index2 in range(len(X[index])):
            list.append(0)
        f2xy.append(list)
    n = len(X)
    for i in range(n):
        m = len(X[i])
        for j in range(m):
            sum = 0 
            for index in range(1,degree+1):
                for index2 in range(1,degree-index+1):
                    sum += (index)*(index2)*matrix[index][index2]*(X[i][j]**(index-1))*(Y[i][j]**(index2-1))  
            f2xy[i][j] = sum
    return(f2xy)

def calcul_f1x(X,Y,C,degree):
    matrix = []
    i = 0
    for index in range(degree+1):
        list = []
        for index2 in range(degree-index+1):
            list.append(C[i])
            i += 1
        matrix.append(list)
    f1x = []
    for index in range(len(X)):
        list = []
        for index2 in range(len(X[index])):
            list.append(0)
        f1x.append(list)
    n = len(X)
    for i in range(n):
        m = len(X[i])
        for j in range(m):
            sum = 0 
            for index in range(1,degree+1):
                for index2 in range(0,degree-index+1):
                    sum += (index)*matrix[index][index2]*(X[i][j]**(index-1))*(Y[i][j]**(index2))  
            f1x[i][j] = sum
    return(f1x)

def calcul_f1y(X,Y,C,degree):
    matrix = []
    i = 0
    for index in range(degree+1):
        list = []
        for index2 in range(degree-index+1):
            list.append(C[i])
            i += 1
        matrix.append(list)
    f1y = []
    for index in range(len(X)):
        list = []
        for index2 in range(len(X[index])):
            list.append(0)
        f1y.append(list)
    n = len(X)
    for i in range(n):
        m = len(X[i])
        for j in range(m):
            sum = 0 
            for index in range(0,degree+1):
                for index2 in range(1,degree-index+1):
                    sum += (index2)*matrix[index][index2]*(X[i][j]**(index))*(Y[i][j]**(index2-1))  
            f1y[i][j] = sum
    return(f1y)


def interpolate(data,degree,k):
    mn = np.min(data, axis=0)
    mx = np.max(data, axis=0)
    X,Y = np.meshgrid(np.linspace(mn[0], mx[0], k), np.linspace(mn[1], mx[1], k))
    XX = X.flatten()
    YY = Y.flatten()
    list = []
    list_z = []
    for index in range(0,degree+1):
        for index2 in range(0,degree-index+1):
            list.append((data[:,0:1]**index)*(data[:,1:2]**index2))
            list_z.append((XX**index)*(YY**index2))
    A=[]
    m = len(list)
    n = len(data)
    for index in range(n):
        new_list = []
        for i in range(m):
            [element] = list[i][index]
            new_list.append(element)
        A.append(new_list)
    M=[]
    n = len(XX)
    for index in range(n):
        new_list = []
        for i in range(m):
            element = list_z[i][index]
            new_list.append(element)
        M.append(new_list)
    C, residues, rank, s = scipy.linalg.lstsq(A, data[:,2])
    Z = np.dot(M,C).reshape(X.shape)

    return A,C,X,Y,Z

def calcul_kmoy(data,degree,C,k,visualisation):
    hull = ConvexHull(data[:,0:2])
    if visualisation:
        fig1 =  plt.figure()
        ax = plt.axes(projection='3d')
        plt.title("Interpolation")
    n = len(data[hull.vertices,0])
    contours_x = [data[hull.vertices,0][-1]]
    contours_y = [data[hull.vertices,1][-1]]
    for index in range(n):
        contours_x.append(data[hull.vertices,0][index])
        contours_y.append(data[hull.vertices,1][index])
    n = len(contours_x)
    pts = []
    for index in range(n):
        pts.append((contours_x[index],contours_y[index]))
    poly = Polygon(pts)
    mn = np.min(data, axis=0)
    mx = np.max(data, axis=0)
    X,Y = np.meshgrid(np.linspace(mn[0], mx[0], k), np.linspace(mn[1], mx[1], k))
    n = len(X)
    Z = []
    f1x = calcul_f1x(X,Y,C,degree)
    f1y = calcul_f1y(X,Y,C,degree)
    f2xy = calcul_f2xy(X,Y,C,degree)
    f2xx = calcul_f2xx(X,Y,C,degree)
    f2yy = calcul_f2yy(X,Y,C,degree)
    for i in range(n):
        m = len(X[i])
        list = []
        for j in range(m):
            list.append(0)
        Z.append(list)
    for i in range(n):
        m = len(X[i])
        for j in range(m):
            Z[i][j]  = (f2xx[i][j]*f2yy[i][j]-(f2xy[i][j]**2))/((1+f1x[i][j]**2+f1y[i][j]**2)**2)
    Z = np.matrix(Z)
    XX = X.flatten()
    YY = Y.flatten()
    ZZ = Z.flatten()
    ZZ = np.array(ZZ)
    [a] = ZZ
    ZZ = a 
    ZZZ = []
    n = len(XX)
    for index in range(n):
        p = Point((XX[index],YY[index]))
        if poly.contains(p):
            ZZZ.append(ZZ[index])
        else :
            ZZ[index]=0
    n = len(ZZZ)
    sum = 0
    for index in range(n):
        sum += ZZZ[index]
    z_moy = sum/n
    if visualisation:
        plt.plot(contours_x, contours_y, 'r-', lw=2)
        #ax.scatter(data[:,0], data[:,1], data[:,2], c='r', s=5,label='Scanned points')
        ax.scatter(XX, YY, ZZ, c='black', s=0.1,label='Scanned points')
        plt.show()
    return z_moy

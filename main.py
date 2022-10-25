import numpy as np
import matplotlib.pyplot as plt

def readCoordsFromFile(filename):
    Xcoordinates = []
    Ycoordinates = []
    dots = open(filename, 'r')
    while True:
        newdot = dots.readline()
        if newdot == '': break
        else:
            coord = newdot.split(',')
            Xcoordinates.append(int(coord[0]))
            Ycoordinates.append(int(coord[1]))
    return np.array(Xcoordinates), np.array(Ycoordinates)

def generateRandomCoords():
    RNG = np.random.default_rng()
    Cnum = np.random.randint(1, 11)
    Crange = np.arange(-100, 101)
    Xcoords = np.sort(RNG.choice(Crange, size=Cnum, replace=False))
    Ycoords = RNG.choice(Crange, size=Cnum)
    return Xcoords, Ycoords

def computeDegree(X, Y):
    prev_m, degree, flag, straightdots = 0.0, 0, 0, 0
    listflag, listdegree, listm = [flag], [degree], [prev_m]
    for i in range(1, len(Y)):
        m = (Y[i]-Y[i-1])/(X[i]-X[i-1])
        if m == prev_m:
            straightdots += 1
        else:
            if flag == 1 and m <= prev_m:
                flag = -1
                degree += 1
            elif flag == -1 and m >= prev_m:
                flag = 1
                degree += 1
            elif flag == 0:
                degree += (straightdots + 1)
                if degree > 1:
                    flag = 1 if prev_m < m else -1
            straightdots = 0
        prev_m = m
        listdegree.append(degree)
        listflag.append(flag)
        listm.append(prev_m)
    if straightdots > 0 and flag != 0: degree += straightdots
    #print('Verificacion de grados:', listdegree)
    #print('Verificacion de banderas:', listm)
    #print('Verificacion de banderas:', listflag)
    return degree

def buildSystem(X, Y, d):
    Xmatrix, Yvector = [], []
    for i in range(len(Y)):
        Xline = []
        for j in range(d+1):
            Xline.append(X[i]**(d-j))
        Xmatrix.append(Xline)
        Yvector.append(Y[i])
    Xmatrix = np.array(Xmatrix)
    print('Matrix A: \n', Xmatrix)
    Yvector = np.array(Yvector)
    print('Vector b:', Yvector)
    return Xmatrix, Yvector

def printCheck(XM, XMT, XM2, INV, XMY, TOTAL):
    print('Original matrix: \n', XM)
    print('Transpose: \n', XMT)
    print('Product: \n', XM2)
    print('Inverse of product: \n', INV)
    print('Transpose times Y: \n', XMY)
    print('Total: \n', TOTAL)

if __name__ == '__main__':
    option, Xcoords, Ycoords = '', [], []
    while option != 'X':
        option = input('¿Desea leer las coordenadas de un archivo (F), '
                       'generarlas aleatoriamente (R) o salir (X)?: ')
        if option in ('F', 'R'):
            cpairs, result1, result2, result0 = '', '', '', ''
            if option == 'F':
                filename = input('Ingrese el nombre del archivo .txt: ')
                Xcoords, Ycoords = readCoordsFromFile(filename)
            elif option == 'R':
                Xcoords, Ycoords = generateRandomCoords()
            print('Las coordenadas obtenidas son: ')
            for cpair in zip(Xcoords, Ycoords):
                cpairs = cpairs + str(cpair) + ' '
            print(cpairs)
            degree = computeDegree(Xcoords, Ycoords)
            X_M, Y_v = buildSystem(Xcoords, Ycoords, degree)
            quoeflist = np.linalg.lstsq(X_M, Y_v, rcond=None)[0]
            #print('Quoefficients:', quoeflist)
            pnf_list = np.polyfit(Xcoords, Ycoords, degree)
            X_M_T = np.transpose(X_M)
            if np.linalg.det(X_M_T @ X_M) != 0:
                calculist = np.linalg.inv(X_M_T @ X_M) @ (X_M_T @ Y_v)
            else: calculist = np.array([0] * len(Y_v))
            printCheck(X_M, X_M_T, (X_M_T @ X_M), np.linalg.inv(X_M_T @ X_M), (X_M_T @ Y_v), calculist)
            assert len(quoeflist) == len(pnf_list) == len(calculist)
            for q in range(len(quoeflist)):
                if q == degree:
                    result1 += '(' + str(quoeflist[q]) + ')'
                    result2 += '(' + str(pnf_list[q]) + ')'
                    result0 += '(' + str(calculist[q]) + ')'
                else:
                    result1 += '(' + str(quoeflist[q]) + ')X^' + str(degree-q) + ' + '
                    result2 += '(' + str(pnf_list[q]) + ')X^' + str(degree-q) + ' + '
                    result0 += '(' + str(calculist[q]) + ')X^' + str(degree - q) + ' + '
            print('La función polinómica que mejor se ajusta con lstsq es', result1)
            print('La función polinómica que mejor se ajusta con polyfit es', result2)
            print('La función polinómica que mejor se ajusta con analítica es', result0)
            #Rmeasure = np.corrcoef(quoeflist, pnf_list)[0, 1] if len(Y_v) > 1 else 1.0
            #print('La medida de correlación (R) entre ambos métodos es', Rmeasure)
            pnf_space = np.arange(-100, 100, 0.01)
            MDQ_fit = np.polyval(quoeflist, pnf_space)
            Autofit = np.polyval(pnf_list, pnf_space)
            Riskyfit = np.polyval(calculist, pnf_space)
            plt.plot(pnf_space, MDQ_fit, 'r')
            plt.plot(pnf_space, Autofit, 'y')
            plt.plot(pnf_space, Riskyfit, 'g')
            plt.plot(Xcoords, Ycoords, 'bs')
            plt.show()
        else: break


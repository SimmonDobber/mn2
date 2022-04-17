import math
import time


class SetOfEquations:
    def __init__(self, index, N, a1, a2, a3):
        self.matrixSize = N
        self.matrix = self.getMatrix(a1, a2, a3)
        self.vector = self.getVector(index)

    def LUSolve(self):
        startTime = time.time()
        _L, _U = self.LUDecomposition()
        y = self.forwardSubstitution(_L, self.vector)
        x = self.backwardSubstitution(_U, y)
        residuumNorm = self.getResiduumNorm(x)
        elapsedTime = time.time() - startTime;
        return x, residuumNorm, round(elapsedTime, 4)

    def jacobiSolve(self, errorThreshold):
        startTime = time.time()
        x = [0 for i in range(self.matrixSize)]
        x1 = [0 for i in range(self.matrixSize)]
        residuumNorm = self.getResiduumNorm(x)
        initialResiduumNorm = residuumNorm
        iterations = 0
        while residuumNorm > errorThreshold:
            if residuumNorm > 10e6 + initialResiduumNorm:
                print("Dla danego przypadku metoda Jacobiego nie zbiega się.")
                return
            for i in range(0, self.matrixSize):
                s = 0
                for j in range(0, self.matrixSize):
                    if j != i:
                        s += self.matrix[i][j] * x[j]
                x1[i] = (self.vector[i] - s) / self.matrix[i][i]
            for i in range(0, self.matrixSize):
                x[i] = x1[i]
            residuumNorm = self.getResiduumNorm(x)
            iterations += 1
        elapsedTime = time.time() - startTime
        return x, iterations, round(elapsedTime, 4)

    def gaussSeidelSolve(self, errorThreshold):
        startTime = time.time()
        x = [0 for i in range(self.matrixSize)]
        residuumNorm = self.getResiduumNorm(x)
        initialResiduumNorm = residuumNorm
        iterations = 0
        while residuumNorm > errorThreshold:
            if residuumNorm > 10e6 + initialResiduumNorm:
                print("Dla danego przypadku metoda Gaussa-Seidla nie zbiega się.")
                return
            for i in range(0, self.matrixSize):
                s = 0
                for j in range(0, self.matrixSize):
                    if j != i:
                        s += self.matrix[i][j] * x[j]
                x[i] = (self.vector[i] - s) / self.matrix[i][i]
            residuumNorm = self.getResiduumNorm(x)
            iterations += 1
        elapsedTime = time.time() - startTime
        return x, iterations, round(elapsedTime, 4)

    def getMatrix(self, a1, a2, a3):
        matrix = [[0 for x in range(self.matrixSize)] for y in range(self.matrixSize)]
        self.fillDiagonal(matrix, 0, a1)
        self.fillDiagonal(matrix, 1, a2)
        self.fillDiagonal(matrix, 2, a3)
        return matrix

    def fillDiagonal(self, matrix, offset, value):
        for i in range(self.matrixSize - offset):
            matrix[i + offset][i] = value
            matrix[i][i + offset] = value

    def getVector(self, index):
        vector = []
        for i in range(self.matrixSize):
            vector.append(self.getVectorElement(i, index))
        return vector

    def getVectorElement(self, i, index):
        f = (index % 10000) / 1000
        return math.sin(i * (f + 1))

    def LUDecomposition(self):
        _L = [[0 for x in range(self.matrixSize)] for y in range(self.matrixSize)]
        _U = [[0 for x in range(self.matrixSize)] for y in range(self.matrixSize)]
        for i in range(0, self.matrixSize):
            for j in range(0, self.matrixSize):
                if j >= i:
                    _L[j][i] = self.matrix[j][i]
                    for k in range(0, i):
                        _L[j][i] = _L[j][i] - _L[j][k] * _U[k][i]
            for j in range(0, self.matrixSize):
                if j == i:
                    _U[i][j] = 1
                elif j >= i:
                    _U[i][j] = self.matrix[i][j] / _L[i][i]
                    for k in range(0, i):
                        _U[i][j] = _U[i][j] - ((_L[i][k] * _U[k][j]) / _L[i][i])
        return _L, _U

    def forwardSubstitution(self, _L, b):
        y = [0 for i in range(self.matrixSize)]
        for i in range(0, self.matrixSize):
            s = 0
            for j in range(0, i):
                s = s + _L[i][j] * y[j]
            y[i] = (b[i] - s) / _L[i][i]
        return y

    def backwardSubstitution(self, _U, b):
        x = [0 for i in range(self.matrixSize)]
        for i in reversed(range(0, self.matrixSize)):
            s = 0
            for j in range(i + 1, self.matrixSize):
                s = s + _U[i][j] * x[j]
            x[i] = (b[i] - s) / _U[i][i]
        return x

    def getResiduumNorm(self, x):
        residuumNorm = 0
        residuumVector = self.getResiduumVector(x)
        for i in range(0, self.matrixSize):
            residuumNorm = residuumNorm + (residuumVector[i] * residuumVector[i])
        residuumNorm = math.sqrt(residuumNorm)
        return residuumNorm

    def getResiduumVector(self, x):
        residuumVector = [0 for i in range(self.matrixSize)]
        for i in range(0, self.matrixSize):
            s = 0
            for j in range(0, self.matrixSize):
                s = s + self.matrix[i][j] * x[j]
            residuumVector[i] = s - self.vector[i]
        return residuumVector





    


import matplotlib.pyplot as plt
from SetOfEquations import SetOfEquations

index = 184574
a1 = 5 + int((index % 1000) / 100)
a2 = -1
a3 = -1
N = 900 + int((index % 100) / 10) * 10 + int(index % 10)
errorThreshold = 1e-9

setOfEquations = SetOfEquations(index, N, a1, a2, a3)
jacobiResults = setOfEquations.jacobiSolve(errorThreshold)
gaussSeidelResults = setOfEquations.gaussSeidelSolve(errorThreshold)

print("Jacobi method results:")
print("Time:", jacobiResults[2], 's')
print("Iterations:", jacobiResults[1], '\n')

print("Gauss-Seidel method results:")
print("Time:", gaussSeidelResults[2], 's')
print("Iterations:", gaussSeidelResults[1], '\n')

setOfEquations = SetOfEquations(index, N, 3, a2, a3)

jacobiResults = setOfEquations.jacobiSolve(errorThreshold)
gaussSeidelResults = setOfEquations.gaussSeidelSolve(errorThreshold)

LUDecompositionResults = setOfEquations.LUSolve()

print("")
print("LU Decomposition method residuum norm:", LUDecompositionResults[1])

N = [100, 500, 1000, 2000, 3000]
timeLUDecomposition = []
timeJacobi = []
timeGaussSeidel = []

for n in N:
    setOfEquations = SetOfEquations(index, n, a1, a2, a3)
    timeLUDecomposition.append(setOfEquations.LUSolve()[2])
    timeJacobi.append(setOfEquations.jacobiSolve(errorThreshold)[2])
    timeGaussSeidel.append(setOfEquations.gaussSeidelSolve(errorThreshold)[2])

plt.plot(N, timeLUDecomposition, label='LU Decomposition', color='green', linewidth=1)
plt.plot(N, timeJacobi, label='Jacobi', color='red', linewidth=1)
plt.plot(N, timeGaussSeidel, label='Gauss-Seidel', color='blue', linewidth=1)
plt.legend(loc="upper left")

plt.xlabel('matrix size')
plt.ylabel('time [s]')
plt.yscale('log')
plt.show()










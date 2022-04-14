from SetOfEquations import SetOfEquations

index = 184574
a1 = 5 + int((index % 1000) / 100)
a2 = -1
a3 = -1
errorThreshold = 1 / 1e9
setOfEquations = SetOfEquations(index, a1, a2, a3)
print(setOfEquations._LUSolve())



load library.sage

branes = ['A','A','A','A','B','C']
algebraType = "Dn"
k = len(branes)
iMat = constructIntersectionMatrix(branes)

print "Root Search:"
roots = junctionSearch(-2,(0,0),1) # I've verified this box size is big enough
simpleroots = [roots[18],roots[14],roots[12],roots[13]]
roots.append([0,0,0,0,0,0])
cartanMatrix = computeCartanMatrix(simpleroots)
print "Cartan Matrix:"
print cartanMatrix
rootsG = generateWeightGraph(roots,simpleroots) # generate a weight graph
# plot in 3d in sage with rootsG.plot3d(vertex_size=.0125,edge_size=.006)

posroots = getPositiveRoots(roots,simpleroots)
rootlevels = freudenthalRecursionFormula(findHighestWeight(roots,simpleroots),posroots,simpleroots)

# finds the 8v, 8s, and 8c related by triality
v8 = junctionSearch(-1,(1,0),2)
s8 = junctionSearch(-1,(0,1),2)
c8 = junctionSearch(-1,(1,1),2)


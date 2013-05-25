load library.sage

branes = ['A','A','A','A','A','A','A','B','C','C']
k = len(branes)
algebraType = "E8"
iMat = constructIntersectionMatrix(branes)

print "Root Search:"
# roots = smartJunctionSearch(-2,(0,0),2,4,2)
# used above search to get roots. now load in
load E8roots.sage
simpleroots = [roots[183],roots[156],roots[140],roots[130],roots[124],roots[122],roots[120],roots[121]] 
cartanMatrix = computeCartanMatrix(simpleroots)
print "Cartan Matrix:"
print cartanMatrix 
rootsG = generateWeightGraph(roots,simpleroots)
rootsG.plot3d(vertex_size=.0125,edge_size=.004)
roots.append([0,0,0,0,0,0,0,0,0,0])
rootsG = generateWeightGraph(roots,simpleroots)
rootsG.plot3d(vertex_size=.0125,edge_size=.004)
highestroot = findHighestWeight(roots,simpleroots)
posroots = getPositiveRoots(roots,simpleroots)


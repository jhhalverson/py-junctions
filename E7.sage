load library.sage

branes = ['A','A','A','A','A','A','B','C','C']
algebraType = "E7"
k = len(branes)
iMat = constructIntersectionMatrix(branes)

print "Root Search:"
roots = smartJunctionSearch(-2,(0,0),1,3,2)
roots.append([0,0,0,0,0,0,0,0,0])
simpleroots = [roots[99],roots[83],roots[73],roots[67],roots[65],roots[63],roots[64]] 
cartanMatrix = computeCartanMatrix(simpleroots)
print "Cartan Matrix:"
print cartanMatrix 

rootsG = generateWeightGraph(roots,simpleroots)

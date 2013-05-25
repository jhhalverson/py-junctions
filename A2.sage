load library.sage

branes = ['A','A','A'] # set up the vanishing cycles
k = len(branes)
iMat = constructIntersectionMatrix(branes) # construct the intersection matrix

print "Root Search:"
roots = junctionSearch(-2,(0,0),2) # search for roots
simpleroots = [roots[4],roots[3]] # these are good simple roots
cartanMatrix = computeCartanMatrix(simpleroots)
print "Cartan Matrix:"
print cartanMatrix

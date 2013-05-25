load library.sage

algebraType = 'An'

branes = ['A','A','A','A']
k = len(branes)
iMat = constructIntersectionMatrix(branes)

print "Root Search:"
roots = junctionSearch(-2,(0,0),2) # I've verified this box size is big enough
roots.append([0,0,0,0])
simpleroots = [roots[9],roots[7],roots[6]]
posroots = getPositiveRoots(roots,simpleroots)
cartanMatrix = computeCartanMatrix(simpleroots)
print "Cartan Matrix:"
print cartanMatrix

the4 = junctionSearch(-2,(2,0),2)


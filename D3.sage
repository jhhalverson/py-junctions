# D3 and A3 are isomorphic as algebras, but have completely different canonical
# bases following Zwiebach and DeWolfe's notation. Do we get the appropriate structures?

load library.sage
algebraType = "Dn" # makes computing self intersections much more efficient

branes = ['A','A','A','B','C']
k = len(branes)
iMat = constructIntersectionMatrix(branes)

# Standard initialization type stuff
print "Root Search:"
roots = junctionSearch(-2,(0,0),1) # I've verified this box size is big enough
roots.append([0,0,0,0,0])
simpleroots = [[1,-1,0,0,0],[0,1,-1,0,0],[0,1,1,-1,-1]]
posroots = getPositiveRoots(roots,simpleroots)
cartanMatrix = computeCartanMatrix(simpleroots)
print "Cartan Matrix:"
print cartanMatrix

# this junction search for AAAA gives the four, not the six = antisymmetric
# interestingly, here it gives the six
theSix = junctionSearch(-1,(1,0),3)
theSixHW = findHighestWeight(theSix,simpleroots)
theSixFreud = freudenthalRecursionFormula(theSixHW,posroots,simpleroots)
print "Data for the d=6 rep"
printFreudenthalData(theSixFreud)


# since the six has Dynkin labels (0,1,0) and the d=20 rep is (0,2,0), mult six HW by two
the20HW = multV(2,theSixHW)
the20Freud = freudenthalRecursionFormula(the20HW,posroots,simpleroots)
print "\n\n\nThis rep. should be the d=20 box rep"
printFreudenthalData(the20Freud)
print "dim is ", countFreudenthal(the20Freud)

# indeed, we found the twenty dimensional representation


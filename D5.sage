load library.sage
algebraType = "Dn" # makes computing self intersections much more efficient

branes = ['A','A','A','A','A','B','C']
k = len(branes)
iMat = constructIntersectionMatrix(branes)

print "Root Search:"
roots = junctionSearch(-2,(0,0),1) # I've verified this box size is big enough
roots.append([0,0,0,0,0,0,0])
simpleroots = [roots[32],roots[26],roots[22],roots[20],roots[21]]
posroots = getPositiveRoots(roots,simpleroots)
highestroot = findHighestWeight(roots,simpleroots)
cartanMatrix = computeCartanMatrix(simpleroots)
print "Cartan Matrix:"
print cartanMatrix

# finds the 16 dimensional spinor representation
print "16 search:"
spinor = junctionSearch(-1,(1,1),2) 
spinorG = generateWeightGraph(spinor,simpleroots) # can plot using spinorG.plot3d()
spinorhw = findHighestWeight(spinor,simpleroots)
print "spinorhw is ", spinorhw
spinorlevels = freudenthalRecursionFormula(spinorhw,posroots,simpleroots)
count = 0
for lev in spinorlevels:
    for w in lev:
        count = count + w[1]
print "Freudenthal found spinor has dim ", count

# knowing that the highest weight of the 126 is twice the highest weight of the
# spinor, find the 126 rep in terms of junctions
onetwosixlevels = freudenthalRecursionFormula(multV(2,spinorhw),posroots,simpleroots)
print "Search for 126? has dim ", countRepDim(onetwosixlevels)


# knowing the the rep with highest weight given by 3 times the highest root has
# dimension 7644, check Freudenthal's formula
seven664 = freudenthalRecursionFormula(multV(3,highestroot),posroots,simpleroots)
count = 0
for lev in seven664:
    for w in lev:
        count = count + w[1]
print "Search for 7644? has dim ", count

load library.sage

branes = ['A','A','A','A','A','B','C','C']
algebraType = "E6"
k = len(branes)
iMat = constructIntersectionMatrix(branes)

print "BEGIN ROOT SEARCH:"
roots = smartJunctionSearch(-2,(0,0),1,2,1)
roots.append([0,0,0,0,0,0,0,0])

simpleroots = [roots[56],roots[46],roots[40],roots[38],roots[36],roots[37]] 
posroots = getPositiveRoots(roots,simpleroots)
cartanMatrix = computeCartanMatrix(simpleroots)
print "Cartan Matrix:"
print cartanMatrix 
rootsG = generateWeightGraph(roots,simpleroots)

# Study the 27, find highest weight, apply Freudenthal to rep with hw = 3*hw27
twoseven = smartJunctionSearch(-1,(1,0),2,3,3)
twosevenG = generateWeightGraph(twoseven,simpleroots)
hw27 = findHighestWeight(twoseven,simpleroots)

# The rep with hw = 3*hw27 has d = 3003. check this
threethousandrep = freudenthalRecursionFormula(multV(3,hw27),posroots,simpleroots)
print "rep with hw = 3*hw27 has dim ", countRepDim(threethousandrep)

# The rep with hw = 3*hw78 has d = 43758 . check this
hw78 = findHighestWeight(roots,simpleroots)
bigrep= freudenthalRecursionFormula(multV(3,hw78),posroots,simpleroots)
print "rep with hw = 3*hw78 has dim ", countRepDim(bigrep)

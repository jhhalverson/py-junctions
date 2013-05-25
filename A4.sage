load library.sage

branes = ['A','A','A','A','A']
k = len(branes)
iMat = constructIntersectionMatrix(branes)

print "Root Search:"
roots = junctionSearch(-2,(0,0),1) # I've verified this box size is big enough
roots.append([0,0,0,0,0])
simpleroots = [roots[16],roots[13],roots[11],roots[10]]

# does a junction search and finds the five and highest weight
five = junctionSearch(-1,(1,0),1,1)
fivehw = findHighestWeight(five,simpleroots)

# does a junction search and finds the five bar and highest weight
fivebar = junctionSearch(-1,(-1,0),1,1)
fivebarhw = findHighestWeight(fivebar,simpleroots)

# does a junction search and finds the ten and highest weight
ten = junctionSearch(-2,(2,0),2,2)
tenhw = findHighestWeight(ten,simpleroots) 

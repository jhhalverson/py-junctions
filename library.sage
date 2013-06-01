#####################################################
# library.sage: backbone for junction computations ##
#####################################################

# This file is the library file which contains all of the basic
# function which perform computations in the paper 
# "Matter from Geometry without Resolution"
# By A. Grassi, J. Halverson, and J. Shaneson
# It is divided up into a few basic sections.

# This library is meant to be loaded and used within
# an "algebra initialization file", e.g. D5.sage.
# See D5.sage for a well commented example.

# A few notes on choices made:
# 1) For speed, vectors (and therefore junctions) are
# representated as lists instead of matrix objects.




#####################################################
# BASIC INITIALIZATION                             ##
#####################################################
import networkx as nx # import graph theory library for weight diagram plotting
import matplotlib.pyplot as plt 

algebraType = None # can label specific algebras to speed up computations

# defines standard notions of A-brane, B-brane, and C-branes
# in terms of the vanishing cycles above the corresponding
# components of the discriminant..
piA = (1,0)
piB = (1,-1)
piC = (1,1)
pi01 = (0,1)
pinA = (-1,0)
pinB = (-1,1)
pinC = (-1,-1)
pin01 = (0,-1)

cycle = {'A': piA, 'B': piB, 'C': piC, '01': pi01, 'Y':pi01,'nA': pinA, 'nB': pinB, 'nC': pinC, 'n01': pin01, 'nY':pin01} # create a useful dictionary



#####################################################
# BASIC FUNCTIONS                                  ##
#####################################################

# Computes the product of two one-cycles pi1 and pi2 on a torus
def oneCycleProduct(pi1,pi2):
    return pi1[0]*pi2[1] - pi1[1]*pi2[0]

# Computes the intersection matrix for a given set of branes, listed in terms
# of labels in an initilization file (e.g. E8.sage). If a wider variety of (p,q)-branes
# is desired, just edit the "cycle" dictionary above accordingly.
def constructIntersectionMatrix(branes):
    k = len(branes)
    entries = [[0 for i in range(k)] for j in range(k)]
    for i in range(k): entries[i][i] = -1
    for i in range(k):
        for j in range(i+1,k): 
            entries[j][i] = oneCycleProduct(cycle[branes[i]],cycle[branes[j]])/2
            entries[i][j] = oneCycleProduct(cycle[branes[i]],cycle[branes[j]])/2
    return entries

# Computes the self-intersection of the junction j. 
# If the algebraType is defined, uses a speeded up algorithm.
# If not, computes using simple matrix multiplication.
def selfIntersect(j):
    if algebraType == "Dn": # Faster computations for Dn. Assumes Dn = A^n BC.
        sum = 0
        nbranes = len(branes)
        numAs = nbranes - 2 
        for i in range(numAs):
            sum = sum -j[i]^2
            sum = sum + j[i]*(j[nbranes-1]-j[nbranes-2])
        sum = sum - (j[nbranes-2]-j[nbranes-1])^2
        return sum    
    elif algebraType == "E6": # Faster computations for E6. Assumes E6 = A^5 BC^2.
        sum = 0
        numAs = 5
        nb = 8
        for i in range(numAs):
            sum = sum -j[i]^2
            sum = sum + j[i]*(j[nb-1]+j[nb-2]-j[nb-3])
        sum = sum - j[nb-1]^2 - j[nb-2]^2 - j[nb-3]^2 + 2*j[nb - 3]*(j[nb-1]+j[nb-2])
        return sum    
    elif algebraType == "E7":  # Similar
        sum = 0
        numAs = 6
        nb = 9
        for i in range(numAs):
            sum = sum -j[i]^2
            sum = sum + j[i]*(j[nb-1]+j[nb-2]-j[nb-3])
        sum = sum - j[nb-1]^2 - j[nb-2]^2 - j[nb-3]^2 + 2*j[nb - 3]*(j[nb-1]+j[nb-2])
        return sum    
    elif algebraType == "E8":  # Similar2
        sum = 0
        numAs = 7
        nb = 10
        for i in range(numAs):
            sum = sum -j[i]^2
            sum = sum + j[i]*(j[nb-1]+j[nb-2]-j[nb-3])
        sum = sum - j[nb-1]^2 - j[nb-2]^2 - j[nb-3]^2 + 2*j[nb - 3]*(j[nb-1]+j[nb-2])
        return sum    
    else: return intersect(j,j)


# Computes the intersection of two junctions j1 and j2. 
# Assumes iMat is appropriately defined elsewhere, e.g. algebra initialization file.
def intersect(j1,j2):
    sum = 0
    for i in range(len(j1)):
        for j in range(len(j2)):
            sum = sum + j1[i]*j2[j]*iMat[i][j]

    return sum
    
# Computes the asymptotic charge of a junction J. 
# Assumes variable branes is appropriately defined.    
def asymptotic(J):
    if type(branes[0]) == str: # specified branes in A,B,C form
        p,q = 0,0
        for i in range(len(J)):
            j = J[i]
            p = p + j*cycle[branes[i]][0]
            q = q + j*cycle[branes[i]][1]
    return (p,q)
    #return (J[0]+J[1]+J[2]+J[3]+J[4]+J[5]+J[6],-J[5]+J[6])

# Computes and returns the Cartan Matrix, given the simple roots.
def computeCartanMatrix(simpleroots):
    cm = [[intersect(j1,j2) for j1 in simpleroots] for j2 in simpleroots]
    return matrix(QQ,cm)

# Searches for junctions such that A-brane, B-brane, and C-brane charge each get their own box.
# Options include a progress meter and the ability to compute rank.
def smartJunctionSearch(selfint,asymp,abox,bbox,cbox,progress=True,computerank=False):
    junctions = []
    junctionboxa, junctionboxb, junctionboxc = [],[],[]
    generateJunctionBox([],branes.count('A'),abox,junctionboxa)
    generateJunctionBox([],branes.count('B'),bbox,junctionboxb)
    generateJunctionBox([],branes.count('C'),cbox,junctionboxc)
    
    # put candidate junctions together in one box
    junctionbox = []
    for ja in junctionboxa:
        for jb in junctionboxb:
            for jc in junctionboxc:
                junctionbox.append(ja + jb + jc)
 
    print "\nSearch: (J,J)=", selfint,", (p,q) = ", asymp, " ", len(junctionbox), " candidate junctions"
    l = len(junctionbox)
    mod  = None
    if l > 1000000: mod = 1000000
    elif l > 100000: mod = 100000
    else: mod=10000
    for i in range(l):
        J = junctionbox[i]
        if progress and i % mod == 0: print i, " of ", l
        if selfIntersect(J) == selfint and asymptotic(J) == asymp:
            junctions.append(J)
    print "There are ", len(junctions), " such junctions."
    if computerank: print "They span a rank ", computeRank(junctions), " lattice."
    return junctions

# Searching for junctions of a give self-intersection and asymptotic charge with the given box size.
# Options include a progress meter and the ability to compute rank.
def junctionSearch(selfint,asymp,box,progress=True,computerank=False,output=True):
    junctions = []
    junctionbox = []
    generateJunctionBox([],k,box,junctionbox)
    if output: print "\nSearch: (J,J)=", selfint,", (p,q) = ", asymp, " ", len(junctionbox), " candidate junctions"
    l = len(junctionbox)
    if l > 1000000: mod = 1000000
    elif l > 100000: mod = 100000
    else: mod=10000
    for i in range(l):
        J = junctionbox[i]
        if progress and output and i % mod == 0: print i, " of ", l
        if selfIntersect(J) == selfint and asymptotic(J) == asymp:
            junctions.append(J)
    if output: print "There are ", len(junctions), " such junctions."
    if len(junctions) != 0 and computerank and output: print "They span a rank ", computeRank(junctions), " lattice."
    return junctions

# Method used to generate junction boxes. See usage in other methods.    
def generateJunctionBox(myList,numbranes,size,biglist,onlypos = False):
    if len(myList) == numbranes:
        biglist.append(myList)
        return

    if not onlypos:
        for i in range(-size,size+1):
            my = myList[0:]
            my.append(i)
            generateJunctionBox(my,numbranes,size,biglist,onlypos=onlypos)
    else:
        for i in range(0,size+1):
            my = myList[0:]
            my.append(i)
            generateJunctionBox(my,numbranes,size,biglist,onlypos=onlypos)
    

# Method to generate boxes.
def generateBox(myList,length,size,biglist,onlyPos=False):
    generateJunctionBox(myList,length,size,biglist,onlypos=onlyPos)

# Puts the given set of weights in a matrix and computes its rank.
def computeRank(weights):
    if type(weights[0]) == list:
        w2 = [matrix(QQ,weight) for weight in weights]
        weights = w2
    tList = []
    for weight in weights: tList.append(weight.list())
    mat = matrix(QQ,tList)
    return mat.rank()

# Searches across asymptotic (p,q) charges for junctions with the given self-intersection.
def pqSearch(maxpq,boxsize,progress=True,computeRank=False):
    print "\n\n\nBEGIN SEARCH OVER MANY (p,q) VALUES"
    for p in range(-maxpq,maxpq+1):
        for q in range(-maxpq,maxpq+1):
            junctionSearch(-1,(p,q),boxsize,progress,computeRank)

# Searching across self-intersections for junctions with the given asymptotic (p,q) charge.
def jjSearch(minjj,maxjj,pq,boxsize,progress=True,computeRank=False,output=True):
    print "\n\n\nBEGIN SEARCH OVER MANY (J,J) VALUES"
    toret = []
    for jj in range(minjj,maxjj+1):
        toret.append(junctionSearch(jj,pq,boxsize,progress,computeRank,output))
    return toret

# Finds the highest weight given a set of weights and the simple roots
def findHighestWeight(weights,simpleroots):
    for weight in weights:
        hw = true
        for root in simpleroots:
            s = [root[i]+weight[i] for i in range(len(root))] 
            if s in weights:
                hw = false
                break
        if hw==true: return weight

# A simple method for getting the level diagram of a set of weights, given the simple roots.
def getLevels(weights,simpleroots):
    hw = findHighestWeight(weights,simpleroots)
    levels, curLevel, foundWeight =[[hw]], 0, True
    while foundWeight == True:
        foundWeight = False
        newlevel = []
        for weight in levels[curLevel]:
             for root in simpleroots:
                 newweight = [weight[k] - root[k] for k in range(len(root))]
                 if newweight in weights and newweight not in newlevel:
                     newlevel.append(newweight)
                     foundWeight=True
        if foundWeight: levels.append(newlevel)
        curLevel = curLevel + 1

    return levels

# Takes a set of weights and the set of simpleroots, uses the getLevels() method, and prints it out for a LaTeX tabular environment.
def printLevelsForLatex(weights,simpleroots):
    levels = getLevels(weights,simpleroots)
    for lev in levels:
        levelstring = str(levels.index(lev)) + "&" + str(len(lev)) + "&" 
        for weight in lev:
            levelstring = levelstring + str(tuple(weight)) + "\,\,\,"
        print levelstring + " \\\\"

# Takes a set of weights and simpleroots and generates the weight diagram as a graph object.
def generateWeightGraph(weights,simpleroots,networkx=False):
    hw = findHighestWeight(weights,simpleroots)
    G = None
    if networkx: G = nx.Graph() # N.B. never implemented the networkx properly, I don't think
    else: G = Graph()

    for weight in weights: 
        if networkx: G.add_node(str(weight))
        else: G.add_vertex(str(weight))
        
    for J in weights:
        for root in simpleroots:
            Jproot, Jmroot = [J[k] + root[k] for k in range(len(root))] , [J[k] - root[k] for k in range(len(root))]
            # the next two lines were causing a double counting.
            # since we start with highest root, only need to subtract
            #if Jproot in weights:
            #    G.add_edge(str(J),str(Jproot),"nolabel")
            if Jmroot in weights:
                G.add_edge(str(J),str(Jmroot),"nolabel")

    return G

# Generates the positive roots from the given set of roots and simpleroots. 
def getPositiveRoots(roots,simpleroots):
    ### some initialization
    highestweight = findHighestWeight(roots,simpleroots)
    if highestweight == None: return 
    weightlevels = [[highestweight]] # contains level information and mults, will output this
    weightdict = [highestweight] # 
    ignoreweights = [] # keep track of weights we know to ignore (e.g. mult 0)
    level, keepGoing = 0, True
    zerovec = [0 for i in range(len(roots[0]))]
    posroots = [highestweight]
    
   
    # loop through
    while keepGoing:
        clevel = weightlevels[level] # will subtract pos roots sroots from current level
        nlevel = [] # building weights at next level
        for w in clevel: # go through weights at level
            for sr in simpleroots: # select an sr to subtract
                weight = [w[i] - sr[i] for i in range(len(w))] # form new weight
                
                # make sure not computed already
                if weight == zerovec: # done
                    keepGoing = False
                elif weight not in ignoreweights and str(weight) not in weightdict: 
                     if weight in roots and weight not in posroots: 
                         weightdict.append(weight)
                         posroots.append(weight)
                         nlevel.append(weight)
                     else: ignoreweights.append(weight)
        if clevel == []: keepGoing=False # prevents from choking if simpleroots aren't simple roots
        if keepGoing:
            weightlevels.append(nlevel)
            level = level + 1
    
    return posroots

                
# Given the roots and the rank of the algebra, find the possible sets of simple roots.
# This is crude and could definitely be optimized.
def getSimpleRootSets(roots,r, numpos):
    sets = []
    box = []
    generateBox([],r, len(roots)-1, box, onlyPos=True)
    cand = []
    for b in box:
        broots = [roots[i] for i in b]
        isCand = True
        for i in b:
            if multV(-1,roots[i]) in broots or b.count(i) > 1:
                isCand = False
                break

        if isCand: cand.append(b)

    cut1 = cand
    cand = []
    for i in range(len(cut1)):
        # if i % 100 == 0: print i, len(cut1)
        hasSame = False
        for j in range(i+1,len(cut1)):
            if set(cut1[i]) == set(cut1[j]): 
                hasSame = True
                break
        if hasSame == False: cand.append(cut1[i])
            

    #print "There are ", len(cand), " candidate sets of simple roots"
    for can in cand:
        #print cand.index(can), [roots[i] for i in can]
        temppos = getPositiveRoots(roots,[roots[i] for i in can])
        if temppos != None and len(temppos) == numpos:
            #print len(temppos)
            sets.append([roots[i] for i in can])
    return sets


# Given the highest weight of a representation, the positive roots, and the simple roots,
# this method uses Freudenthal's recursion formula to generate the level diagram of the
# representation, including multiplicities.
def freudenthalRecursionFormula(highestweight,proots,sroots):
    ### some initialization
    weightlevels = [[[highestweight,1]]] # contains level information and mults, will output this
    weightdict = {str(highestweight):1} # keep mults in easy to access dict
    ignoreweights = [] # keep track of weights we know to ignore (e.g. mult 0)
    weylvec = [0 for i in range(len(highestweight))] # compute weyl vector
    for pr in proots:
        for i in range(len(weylvec)): weylvec[i] = weylvec[i]+1/2*pr[i]
    level, keepGoing = 0, True
   
    # go through level by level
    while keepGoing:
        clevel = weightlevels[level] # will subtract pos roots sroots from current level
        nlevel = [] # building weights at next level
        for w in clevel: # go through weights at level
            for sr in sroots: # select an sr to subtract
                weight = [w[0][i] - sr[i] for i in range(len(w[0]))] # form new weight
                
                # make sure not computed already
                if weight not in ignoreweights and str(weight) not in weightdict.keys(): 
                     ### BEGIN COMPUTATION OF WEIGHT MULTIPLICITY
                     dim = 0
                     for pr in proots:
                         j, nonzerodim = 1, True
                         while nonzerodim:
                             nextup = [weight[i] + j * pr[i] for i in range(len(weight))]
                             if str(nextup) in weightdict.keys(): # if here, this has nextup has non-zero dim
                                 dim = dim + 2*intersect(nextup,pr)*weightdict[str(nextup)]
                                 j = j + 1
                             else: nonzerodim=False
                     v1 = [highestweight[i] + weylvec[i] for i in range(len(weight))]
                     v2 = [weight[i] + weylvec[i] for i in range(len(weight))]
                     den =  (selfIntersect(v1) - selfIntersect(v2))
                     if den != 0: 
                         dim = dim/den
                         if dim != 0:
                             weightdict[str(weight)] = dim
                             nlevel.append([weight,dim])
                         else: ignoreweights.append(weight)
                     else: ignoreweights.append(weight)
                #     ### END COMPUTATION OF WEIGHT MULTIPLICITY
                    
        if len(nlevel) == 0: keepGoing = False # no new weights at the new level, done        
        else: # found new weights. append and keep going
            weightlevels.append(nlevel)
            level = level + 1

    return weightlevels

# Prints the output of the method freudenthalRecursionFormula() in a nice manner.
def printFreudenthalData(freud,selfint=True):
    for lev in freud:
        print "Level ", freud.index(lev)
        for r in lev:
            if selfint: print "Junction: ", r[0], "       (J,J): ",selfIntersect(r[0]), "       Multiplicity: ", r[1]
            else: print "Junction: ", r[0], "       Multiplicity: ", r[1]

# Prints out the output of the Freudenthal recursion formula data in a nice manner for LaTeX.
def printFreudenthalDataLatex(freud,selfint=True):
    for lev in freud:
        for r in lev:
            if selfint: print "$",freud.index(lev), "$ & $", r[0], "$ & $",selfIntersect(r[0]), "$ & $", r[1], "$ \\\\"
            else: "$", printfreud.index(lev), "$ & $", r[0], "$ & $", r[1], "$ \\\\"

# Counts the dimension of the representation given the output of freudenthalRecursionFormula().
def countFreudenthal(freud):
    count = 0
    for lev in freud:
        for r in lev: count = count + r[1]
    return count

# returns the matrix that maps the junction basis to the Dynkin basis
# in the standard representation of ADE groups.
def getDynkinMap(sroots):
    return -1*matrix(sroots)*matrix(iMat)

        
def mapToDynkin(J,sroots):
    m1 = getDynkinMap(sroots)*matrix(QQ,J).transpose()
    return list(list(m1.transpose())[0])

def findHighestWeights(boxsize):
    junctionbox = []
    generateJunctionBox([],len(branes),boxsize,junctionbox)
    myDict = {}
    for J in junctionbox:
        temp = list(mapToDynkin(J))[0]
        allNonNeg = True
        for i in temp:
            if i < 0: allNonNeg = False
        if allNonNeg:
            if str(temp) in myDict.keys():
                myDict[str(temp)].append(J)
            else: myDict[str(temp)] = [J]

    return myDict

# takes output of freudenthal and counts states
def countRepDim(f):
    count = 0
    for lev in f:
        for r in lev:
            count = count + r[1]
    return count

# given a set of one-cycles in branes, e.g. [(1,0),(0,1),(1,0)], computes
# the Picard-Lefschetz action on the one-cycle gam
def picardLefschetz(gam,branes):
    g0, g1 = gam[0],gam[1]
    for brane in branes:
        g0old, g1old = g0,g1
        prod = oneCycleProduct((g0,g1),brane)
        g0 = g0old + prod*brane[0]
        g1 = g1old + prod*brane[1]
    return (g0,g1)

# The next four are just simple algebraic methods for lists
def addV(v1,v2):
    return [v1[i]+v2[i] for i in range(len(v1))]

def subV(v1,v2):
    return [v1[i]-v2[i] for i in range(len(v1))]

def minV(v):
    return [-1*v[i] for i in range(len(v))]

def multV(r,v):
    return [r*v[i] for i in range(len(v))]

def sumV(vecs):
    vec = [0 for i in range(len(vecs[0]))]
    for v in vecs:
        vec = addV(v,vec)
    return vec

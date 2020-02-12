#################
# post-processing routines for zobov
# Last modified Jan 28, 2008, by Mark Neyrinck

import pylab as M
import numpy as N

def logpofr_old(r):
    """
    Returns the nat. logarithm of a fit to the probability of finding
    a zone of density contrast r in a Poisson particle simulation
    """
    
    # 2D
    #p = -2.6*(r-1.)
    # 3D
    p = -5.12*(r-1) - 0.8*(r-1)**2.8
    
    return p

def logpofr(r,whichpofr = '3D'):
    """
    Returns the nat. logarithm of a fit to the probability of finding
    a zone of density contrast r in a Poisson particle simulation
    """

    if (whichpofr == '3Dfudged'):
        p = -5.12*(r-1) - 0.8*(r-1)**4.7
    elif (whichpofr == '2D'):
        p = -2.6*(r-1.)
    elif (whichpofr == '3D'):
        p = -5.12*(r-1) - 0.8*(r-1)**2.8
    else:
        print 'Unrecognized whichpofr in logpofr:',whichpofr
        p = 0.
    
    return p

def mostProbableVoidExtents(prefix, extraoutsuffix = '', densthresh = 0.,whichpofr='3D',plotlabel='',ls='',justdovoids=[], plotcurves=False):
    """
    Truncates large voids at their most probable level.  Outputs a
    plot showing the curve of probabilities at each zone-accretion
    event.

    If the program's thinking about joining zones 1 and 2 together,
    for instance, it compares P(1)P(2) to P(1+2), i.e. the probability
    that both zones 1 and 2 arose from Poisson noise to the
    probability that their union arose from Poisson noise.

    densthresh = an optional parameter keeping links between zones
    below a certain density, e.g. 0.2.  This can still be necessary to
    halt the growth of the largest void, since the ratio of max to min
    densities can be quite large.

    Outputs: .mpve.void returns an ASCII zone list for each output void.
             .mpve.txt returns a file in the same format as .txt (w/o header)
    """
    
    if (densthresh == 0.):
        densthresh = 1e30
    voidfile = prefix+'.void'
    txtfile = prefix+'.txt'
    mpvevoidfile = prefix+extraoutsuffix+'.mpve.void'
    mpvetxtfile = prefix+extraoutsuffix+'.mpve.txt'

    Fvoid = open(voidfile,'r')
    nvoids = int(Fvoid.readline())
    print nvoids,'voids (including possible crazy border voids).'
    Fvoid.close()

    voidsread = M.load(txtfile,skiprows=2)
    nvoids_txt = len(voidsread[:,0])
    vid = N.zeros(nvoids,dtype='int') # void id
    for v in range(nvoids_txt):
        vid[voidsread[v,1]] = v
    lodens_init = voidsread[vid,3]
    vol_init = voidsread[vid,4]
    np_init = voidsread[vid,5]
    r_init = voidsread[vid,9]

    # Read the .void file
    Fvoid = open(voidfile,'r')
    nvoids = int(Fvoid.readline())
    print nvoids,' voids.'

    colors = ['b','g','r','c','m','y','k']

    mostprobadd = N.zeros(nvoids,dtype='int')
    r_mpve = 1.*r_init
    numcurvesplotted = 0
    big_rlist = []
    big_denslist = []
    big_zonestoaddlist = []
    big_totalzonelist = []
    # "big" lists over all voids/zones
    naddsarray = N.zeros(nvoids,dtype=int)

    # read the file
    for v in range(nvoids):
        col = colors[v % 7]
        voidnums = (Fvoid.readline()).split()
        pos = 1
        numzonestoadd, r = int(voidnums[pos]), float(voidnums[pos+1])
        denslist = [lodens_init[v]]
        rlist = [r]
        zonestoaddlist = [] # we're not including [v] in the zones to add
        totalzonelist = [v]

        while (numzonestoadd > 0):
            zonestoadd = map(int, voidnums[pos+2:pos+2+numzonestoadd])
            rnext = float(voidnums[pos+numzonestoadd+3])
            dens = lodens_init[v]*r

            if (dens < densthresh):
                rlist.append(rnext)
                denslist.append(dens)
                zonestoaddlist.append(zonestoadd) # not extend, since we want a list of lists
                totalzonelist.extend(zonestoadd) # for this one, we just lump all the zones together

            pos += numzonestoadd+2
            numzonestoadd, r = int(voidnums[pos]), float(voidnums[pos+1])
            #when we write it, r will be the next r
                    
        big_rlist.append(rlist)
        big_denslist.append(denslist)
        big_zonestoaddlist.append(zonestoaddlist)
        big_totalzonelist.append(totalzonelist)
        naddsarray[v] = len(denslist)-1
    Fvoid.close()

    if len(justdovoids) == 0:
        justdovoids = range(nvoids)

    for v in justdovoids:
        denslist = big_denslist[v]
        if ((naddsarray[v] > 0)*(denslist[0] < densthresh)):
            col = colors[v % 7]
            rlist = big_rlist[v]
            
            zonestoaddlist = 1*big_zonestoaddlist[v]
           
            score = logpofr(rlist[0],whichpofr=whichpofr)
            scorelist = [score]
            
            for adds in range(naddsarray[v]):
                dens = denslist[adds+1] # This is the linking density for add "adds."
                
                if (dens < densthresh):
                    zonestoadd = zonestoaddlist[adds]
                    
                    # see if a subset of zonestoadd exists as a whole void v1;
                    # if so, we don't want to include v1's subvoids' probabilities
                    
                    for v1 in range(nvoids):
                        # There might be a more elegant way to do this.
                        # We see if listtocompare is entirely in zonestoadd;
                        # if so, we take out all but v1 (the first element in zonestoadd.
                        listtocompare = big_totalzonelist[v1]
                        lenlisttocompare = len(listtocompare)
                        if (lenlisttocompare > 1): #otherwise, we might as well leave it in
                            zonesincommon = sum(int(z in zonestoadd) for z in listtocompare)
                            if (zonesincommon == lenlisttocompare):
                                #remove listtocompare from zonestoadd, except for first elt.
                                for i in range(1,len(listtocompare)):
                                    zonestoadd.remove(listtocompare[i])
                                
                    score += logpofr(rlist[adds+1],whichpofr=whichpofr) - logpofr(rlist[adds],whichpofr=whichpofr)
                    print score, rlist[adds+1],rlist[adds],zonestoadd
                    for z in zonestoadd:
                        score -= logpofr(r_init[z],whichpofr=whichpofr)
                    print score
                
                    scorelist.append(score)
                    
            scorearray = M.array(scorelist)
            densarray = M.array(denslist)
            #rnextarray = M.array(rnextlist)
            mostprobadd[v] = N.where(scorearray == min(scorearray))[0]
            r_mpve[v] = rlist[mostprobadd[v]]
            print v,': most probadd:',mostprobadd[v]
            print ls
            if (plotcurves):
                M.plot(scorearray,ls,label=plotlabel,linewidth=2)
                numcurvesplotted += 1

    if (numcurvesplotted > 0):
        M.show()

    # Do it again, but save the most-probable void extents; output

    Fvoid = open(voidfile,'r')
    nvoids = int(Fvoid.readline())
    Fmpvevoid = open(mpvevoidfile,'w')
    vol_mpve = 1.*vol_init
    np_mpve = 1*np_init
    nz_mpve = 0*np_init

    Fmpvevoid.write(str(nvoids)+'\n')
    for v in range(nvoids):
        zonestoaddlist = 1*big_zonestoaddlist[v]
        zonelist = [v]
        for adds in range(mostprobadd[v]):
            zonestoadd = zonestoaddlist[adds]
            zonelist.extend(zonestoadd)
            vol_mpve[v] += M.sum(vol_init[M.array(zonestoadd)])
            np_mpve[v] += M.sum(np_init[M.array(zonestoadd)])
                    
        nz_mpve[v] = len(zonelist)
        for z in zonelist:
            Fmpvevoid.write(str(z)+' ')
        Fmpvevoid.write('\n')

    Fvoid.close()
    Fmpvevoid.close()

    #replace the old array for rewriting
    voidsread[vid,6] = 1*nz_mpve
    voidsread[vid,7] = 1.*vol_mpve
    voidsread[vid,8] = 1*np_mpve
    voidsread[vid,9] = 1.*r_mpve
    voidsread[vid,10] = M.exp(logpofr(r_mpve,whichpofr=whichpofr))

    index = (-voidsread[:,9]).argsort()

    #read x,y,z
    

    nvout = 1
    Fmpvetxt = open(mpvetxtfile,'w')
    for i in index: # index is the output of the sort by prob. threshold
        Fmpvetxt.writelines(str(nvout)+' '+str(int(voidsread[i,1]))+' '+\
                            str(int(voidsread[i,2]))+' '+str(voidsread[i,3])+' '+\
                            str(voidsread[i,4])+' '+str(int(voidsread[i,5]))+' '+\
                            str(int(voidsread[i,6]))+' '+str(voidsread[i,7])+' '+\
                            str(int(voidsread[i,8]))+' '+str(voidsread[i,9])+' '+\
                            str(voidsread[i,10])+'\n')
        nvout += 1
    Fmpvetxt.close()
                                                                     
    

def useRThreshold(prefix, extraoutsuffix = '', densthresh = 0., rthresh = 1.,whichpofr='3D'):
    """
    In this method, you set a density-contrast threshold rthresh.
    Voids not exceeding that threshold are not mentioned in the output
    file.  Voids stop accreting zones if they encounter another zone
    with r>rthresh.

    densthresh is an optional parameter below which link densities are constrained.

    Outputs: .rthresh.%f.void returns an ASCII zone list for each output void, w/ threshold %f.
             .mpve.%f.txt returns a file in the same format as .txt (w/o header)
    """
    
    if (densthresh == 0.) :
        densthresh = 1e30
    voidfile = prefix+'.void'
    txtfile = prefix+'.txt'
    rthreshvoidfile = prefix+extraoutsuffix+'.rthresh.'+str(rthresh)+'.void'
    rthreshtxtfile = prefix+extraoutsuffix+'.rthresh.'+str(rthresh)+'.txt'

    Fvoid = open(voidfile,'r')
    nvoids = int(Fvoid.readline())
    print nvoids,'voids (including possible crazy border voids).'
    Fvoid.close()

    voidsread = M.load(txtfile,skiprows=2)
    nvoids_txt = len(voidsread[:,0])
    vid = N.zeros(nvoids,dtype='int') # void id
    for v in range(nvoids_txt):
        vid[voidsread[v,1]] = v
    lodens_init = voidsread[vid,3]
    vol_init = voidsread[vid,4]
    np_init = voidsread[vid,5]
    r_init = voidsread[vid,9]

    # Read the .void file

    r_rthresh = 0.*r_init
    vol_rthresh = 1.*vol_init
    np_rthresh = 1*np_init
    nz_rthresh = 0*np_init
    
    Fvoid = open(voidfile,'r')
    nvoids = int(Fvoid.readline())
    Frthreshvoid = open(rthreshvoidfile,'w')

    numpassingrthresh = len(N.where(r_init > rthresh)[0])

    print nvoids,' voids,', numpassingrthresh, ' exceed rthresh.'
    Frthreshvoid.write(str(numpassingrthresh)+'\n')

    for v in range(nvoids):
        voidnums = (Fvoid.readline()).split()
        pos = 1
        numzonestoadd, r = int(voidnums[pos]), float(voidnums[pos+1])
        nadds = 0
        if (r_init[v] > rthresh):
            zonelist = [v]
            r_rthresh[v] = r_init[v] #set it to its original value
        else:
            zonelist = []
        
        while (numzonestoadd > 0):
            zonestoadd = M.array(voidnums[pos+2:pos+2+numzonestoadd]).astype(int)
            dens = lodens_init[v]*r
            rnext = float(voidnums[pos+numzonestoadd+3])
            if (r_rthresh[v] == r_init[v]):
                # If it's eligible for growing, butwe haven't already quit adding zones
                if ((max(r_init[zonestoadd]) > rthresh) or (dens > densthresh)):
                    r_rthresh[v] = r
                    print 'r of ',max(r_init[zonestoadd]),' found; stopping accretion'
                else:
                    zonelist.extend(zonestoadd)
                    vol_rthresh[v] += M.sum(vol_init[M.array(zonestoadd)])
                    np_rthresh[v] += M.sum(np_init[M.array(zonestoadd)])

            pos += numzonestoadd+2
            numzonestoadd, r = int(voidnums[pos]), float(voidnums[pos+1])

        if (r_init[v] >= rthresh):
            for z in zonelist:
                Frthreshvoid.write(str(z)+' ')
            Frthreshvoid.write('\n')

            nz_rthresh[v] = len(zonelist)
            print zonelist

    Fvoid.close()
    Frthreshvoid.close()

    #replace the old array for rewriting
    voidsread[vid,6] = 1*nz_rthresh
    voidsread[vid,7] = 1.*vol_rthresh
    voidsread[vid,8] = 1*np_rthresh
    voidsread[vid,9] = 1.*r_rthresh
    voidsread[vid,10] = M.exp(logpofr(r_rthresh,whichpofr=whichpofr))
    
    index = (-voidsread[:,9]).argsort()

    nvout = 1
    Frthreshtxt = open(rthreshtxtfile,'w')
    for i in index: # index is the output of the sort by prob. threshold
        if (r_rthresh[int(voidsread[i,1])] > 0.):
            Frthreshtxt.writelines(str(nvout)+' '+str(int(voidsread[i,1]))+' '+\
                            str(int(voidsread[i,2]))+' '+str(voidsread[i,3])+' '+\
                            str(voidsread[i,4])+' '+str(int(voidsread[i,5]))+' '+\
                            str(int(voidsread[i,6]))+' '+str(voidsread[i,7])+' '+\
                            str(int(voidsread[i,8]))+' '+str(voidsread[i,9])+' '+\
                            str(voidsread[i,10])+'\n')
            nvout += 1
    Frthreshtxt.close()

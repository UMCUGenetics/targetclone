import random
import math

def randvec(vecsum, N, maxval, minval):
    if N*minval > vecsum or N*maxval < vecsum:
        raise ValueError ('Cannot create desired vector')

    indices = list(range(N))
    vec = [random.randint(minval,maxval) for i in indices]
    diff = sum(vec) - vecsum # we were off by this amount.

    #Iterate through, incrementing/decrementing a random index 
    #by 1 for each value we were off.
    while diff != 0:  
        addthis = 1 if diff > 0 else -1 # +/- 1 depending on if we were above or below target.
        diff -= addthis

        ### IMPLEMENTATION 1 ###
        idx = random.choice(indices) # Pick a random index to modify, check if it's OK to modify
        while not (minval < (vec[idx] - addthis) < maxval):  #operator chaining.  If you don't know it, look it up.  It's pretty cool.
            idx = random.choice(indices) #Not OK to modify.  Pick another.

        vec[idx] -= addthis #Update that index.

        ### IMPLEMENTATION 2 ###
        # random.shuffle(indices)
        # for idx in indices:
        #    if minval < (vec[idx] - addthis) < maxval:
        #        vec[idx]-=addthis
        #        break
        #
        # in situations where (based on choices of N, minval, maxval and vecsum)
        # many of the values in vec MUST BE minval or maxval, Implementation 2
        # may be superior.

    return vec

minorCloneFractions = randvec(10, 3,10,1)
print minorCloneFractions

minorCloneFractions2 = random.sample(minorCloneFractions, 2)
print minorCloneFractions2

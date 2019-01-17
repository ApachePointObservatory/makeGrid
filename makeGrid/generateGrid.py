from __future__ import division, absolute_import

import numpy
import random

from .anneal import Annealer

numpy.seterr(invalid='raise')
# numpy.random.seed(0)

CoverHemisphere = False # if True then generate points across the whole sphere,
# but only use the ones in altitude limits. This leads to a more even distribution
# but fewer points at the altitude limit. It also results in the final # of points being an approximation
# if False then points near the altitude limits tend to be AT the altitude limits,
# which may be an advantage to make sure the the telescope goes below lowest fiducial.

class GenerateAzAltGrid(object):
    def __init__(self, nPts, minAlt=10, maxAlt=70, nIterDistribute=5000, nIterOrder=5000, azMix=60, fracMix=0.5):
        """Run a MCMC routine (specifically Metropolis-Hastings) to find
        A good distribution of points in alt, az space

        @param[in] nPts  number of points to find within the requested altitude bounds. Exact number not guarenteed in solution.
        @param[in] minAlt  degrees, minimum altitude to generate points for
        @param[in] maxAlt  degrees, max altitude to generate points for
        @param[in] azMix  delta-az to mix at -90 boundary
        @param[in] fracMix  fraction of points to mix

        The results are available as two attributes:
        - azAltList: a list of (az, alt) points
        - unorderedAzAltList: a list of (az, alt) points before ordering is applied
        """

        if CoverHemisphere:
            # determine density based on nPts and altitude bounds.
            # solve the double integral- see area by integration: http://en.wikipedia.org/wiki/Sphere
            # Integrate(0, 2pi) Integrate(minPhi, maxPhi) sin theta dTheta dPhi
            # note: phi = [0, pi] radians, but max/min alt are input in degrees with +z = 90!
            # be careful! note r=1 for all calculations
            phi = numpy.pi/2.0 + numpy.radians([minAlt, maxAlt])
            # below is analytic solution to integral
            boundedArea = 2*numpy.pi*(numpy.cos(phi[0]) - numpy.cos(phi[1]))
            areaHemisphere = 2*numpy.pi
            # determine additional points required to maintain desired density on the hemisphere
            self.nPts = nPts * (areaHemisphere/boundedArea)
        else:
            self.nPts = nPts
        self.maxAlt = float(maxAlt)
        self.minAlt = float(minAlt)
        self.azMix = float(azMix)
        self.fracMix = float(fracMix)
        self._orderPoints = OrderPoints()

        self.posVecList = self._getInitialSample(10)
        self.azAltList = self.cart2sph(self.posVecList)
        self.didWrap = False

    def _getInitialSample(self, nSamples=10):
        """Get a sufficiently good starting point for distributePoints

        @param[in] nSamples, int. Number of distributions to try, return the best one
        @return a list of position vectors points
        """
        sample = self.getRandomSample(self.nPts)
        cost, foo = self.computeDistributionCost(sample)
        for ii in xrange(nSamples):
            newSample = self.getRandomSample(self.nPts)
            newCost, foo = self.computeDistributionCost(newSample)
            if newCost < cost:
                sample = newSample[:]
                cost = newCost
        return sample

    def getRandomSample(self, nPts):
        """Randomly sample the surface of a hemisphere in cartesian coordinates

        @param[in] nPts  integer, number of points to generate
        @return 2D numpy array nPts x position vector

        http://mathworld.wolfram.com/SpherePointPicking.html shows how to do it.
        """
        posVecList = []
        while len(posVecList) < nPts:
            # using Marsaglia's method (1972), http://mathworld.wolfram.com/SpherePointPicking.html
            x1 = 2*numpy.random.sample()-1.
            x2 = 2*numpy.random.sample()-1.
            if x1**2+x2**2>=1:
                continue
            x = 2*x1*numpy.sqrt(1-x1**2-x2**2)
            y = 2*x2*numpy.sqrt(1-x1**2-x2**2)
            z = 1-2*(x1**2+x2**2)
            if z < 0:
                # reflect all points on lower hemisphere to upper one
                # so we don't have to solve the whole sphere.
                z = -z
            if not CoverHemisphere:
                alt = 90 - numpy.degrees(numpy.arccos(z))
                if alt < self.minAlt or alt > self.maxAlt:
                    # reject point as out of bounds
                    continue
            posVecList.append([x,y,z])

        return numpy.asarray(posVecList)

    def cart2sph(self, cartArray):
        """Convert cartesian points to az, alt points, throw away points outside altitude bounds

        @param[in] cartArray  numpy 2D array nPts x [x, y, z]
        @return a numpy 2D array of nPts x [az, alt] points

        Use spherical convention with +z = 90 degrees altitude (rather than zero)
        """
        out = []
        for x,y,z in cartArray:
            # atan2 defined on -pi to pi remap such that in degrees south = 0, east = 90
            az = numpy.degrees(numpy.arctan2(y,x) + numpy.pi/2)
            # set +z to 90 degrees this is what the plotting routine expects
            # altitude values range from 0 to pi/2
            # subtracting 90 will put in range -90 to 0
            # multiply by -1 will put in range 0 to 90
            alt = 90 - numpy.degrees(numpy.arccos(z))
            if CoverHemisphere and alt < self.minAlt or alt > self.maxAlt:
                # reject point as out of bounds
                continue
            out.append([az, alt])
        return numpy.asarray(out)

    def distributePoints(self, nIter):
        """Generate a list of points in semi-random order

        Use MCMC (Metropolis-Hastings) iteration to generate a set of points based on minimizing a cost,
        where the cost is proportional to the inverse square of the cartesian distance between one point
        and all others, total cost is the sum of costs for all points.

        @param[in] posVecList  numpy 2D array shape nPts x [x, y, z] initial starting point for algorithm

        The cost of any given set of points (or state) is defined:
            -The cost function between two points is defined as 1/r^2 where r is the distance between points
            -The total cost will be the sum of costs between every combination of points on the hemisphere
            -The costliest point is determined as the point with the highest single cost

        The algorithm works in the following way:

        1. Start with a random, uniform distribution of points on a hemisphere,
            compute it's cost, determine the costliest point (point with higest cost).
            Throughout the routine, keep track of the best seen configuration.
        -- iterate n times --
        2. Replace the costliest point with a random point on the sphere
        3. Measure cost of new configuration
        4. If this cost is lower than the previous, move to this new state
        5. If the cost is higher, move to this new state with a probability chosen to be:
                if random sample [0, 1] > (newCost - currCost)/newCost
                so if the cost is only slightly larger than the previous, it is more likely to move there
        -- stop iteration ---
        6. return the best configuration visited throughout the iteration.
        """
        currCoords = self.posVecList
        currCost, worstIndex = self.computeDistributionCost(currCoords)
        bestCoords = currCoords[:] # will hold the best set of points seen
        bestCost = currCost
        numTookWorst = 0
        for ii in xrange(nIter):
            testCoords = currCoords[:]
            # first get a new random sample
            newCoord = self.getRandomSample(1)
            # replace the previously worst coord with the new coord
            testCoords[worstIndex,:] = newCoord
            testCost, testWorstIndex = self.computeDistributionCost(testCoords)
            # is the cost lower?
            takeNewState = False
            if testCost <= currCost:
                # move to this new state
                takeNewState = True
            else:
                # cost was higher, but still allow some
                # probability of accepting the proposed state
                # build a stochastic model for acceptance based on
                # how much worse this new state is
                costRatio = (testCost - currCost)/testCost
                # lower the cost ratio, the more likely to accept state
                if numpy.random.sample() > costRatio:
                    # accept the proposed state
                    numTookWorst += 1
                    takeNewState = True
            if takeNewState:
                currCost = testCost
                currCoords = testCoords[:]
                worstIndex = testWorstIndex
                # check to see if this is the best state yet
                if currCost < bestCost:
                    # print 'new best cost:', currCost, 'at iteration: ', ii
                    bestCoords = currCoords[:]
                    bestCost = currCost
        # print "took worst cast %0.2f of the time" % (numTookWorst / nIter,)
        self.posVecList = bestCoords
        self.azAltList = self.cart2sph(self.posVecList)
        self.didWrap = False

    def computeDistributionCost(self, posVecList):
        """Determine how good this configuration of coords is based on a cost.

        @param[in] posVecList  2D numpy array nPts x [x, y, z]
        @return totalCost, costlyPointIndex

        the cost metric:
        cost = 1/r^2 where r is determined to be (note, linear not spherical!) distance between any two points

        totalCost = sum of all costs between all points
        costlyPointIndex = the index in posVecList corresponding to the most costly point
        """
        n = len(posVecList)
        di = numpy.diag_indices(n) # indices of diagonal elements
        # big dot product, generates an array nPts x nPts containing cartesian dot product between every point
        bigDot = numpy.dot(posVecList, posVecList.T)
        costMatrix = self.arccosCost(bigDot) # you may choose to use the arccosCost method instead
        # enter zero on the diagonal (no cost from being close to yourself)
        costMatrix[di] = 0.0
        totalCost = numpy.sum(costMatrix)
        # now determine the index in posVecList with the highest cost
        flattenedCost = numpy.sum(costMatrix, 0)
        costlyPointIndex = numpy.argmax(flattenedCost)
        return totalCost, costlyPointIndex

    def cosApproxCost(self, bigDot):
        """Cost function is approximately 0.5/angular separation^2

        posVecA dot posVecB = cos(angular separation)

        cos(theta) ~= 1 - 0.5*theta^2 for small theta

        So angular separation^2 ~= 2 * (1 - cos(angular separation))
        (and we ignore the factor of 2, since it has no effect on the final answer)

        This avoids taking an expensive arc-cosine to measure distance, and it seems to do just fine.
        """
        # bigDot may be slightly larger than 1 due to numerical error;
        # use a value slightly larger than 1 to avoid negative numbers
        halfDistSqMat = (1.0001 - bigDot) # approximately
        costMatrix = 1.0 / halfDistSqMat
        return costMatrix

    def arccosCost(self, bigDot):
        """Solve using arccos, slower, more exact, computes actual cartesian distance betweeen points
        """
        # Contition the matrix. Some of the dot products are numerically > 1 or < -1 so shrink a bit to avoid
        distance = numpy.arccos(bigDot * 0.999)
        costMatrix = distance**-4
        return costMatrix

    def orderPoints(self, nIter):
        """Improve the order of self.azAltList

        The idea is to get an order that is reasonably efficient, but not so nicely ordered
        that effects such as hysteresis are hidden

        The wrap limits are handled as follows:
        - The points are put in order -90, 270
        - Some points in the range azMix about -90 are put into the other wrap,
          so that points are now in the range -90 - azMix/2 to 270 + azMix/2
        - Those points are ordered, treating wrap as a large distance

        @param[in] nIter  the number of iterations
        """
        azAltList = numpy.copy(self.azAltList)

        # if not already done, wrap az into range -90, 270 and mix points around -90
        if not self.didWrap:
            wrappedAzList = azAltList[:,0]
            wrappedAzList = numpy.where(wrappedAzList < -90, wrappedAzList + 360, wrappedAzList)
            wrappedAzList = numpy.where(wrappedAzList > 270, wrappedAzList - 360, wrappedAzList)

            # mix points about -90
            incrAz = -90 + (self.azMix / 2.0)
            indToIncr = numpy.nonzero(wrappedAzList < incrAz)[0]
            numToIncr = int(len(indToIncr) * self.fracMix)
            if numToIncr > 0:
                indToIncr = numpy.random.choice(indToIncr, numToIncr, replace=False)
                wrappedAzList[indToIncr] += 360

            decrAz = 270 - (self.azMix / 2.0)
            indToDecr = numpy.nonzero(wrappedAzList > decrAz)[0]
            numToDecr = int(len(indToDecr) * self.fracMix)
            if numToDecr > 0:
                indToDecr = numpy.random.choice(indToDecr, numToDecr, replace=False)
                wrappedAzList[indToDecr] -= 360

            azAltList[:,0] = wrappedAzList
            self.didWrap = True

        self.azAltList = self._orderPoints(azAltList, nIter=nIter)


class OrderPoints(object):
    """Reorder a collection of points so they can be traversed quickly.

    To use:
    - Construct an OrderPoints
    - Call the object as a function with your list of points

    Credits:
    Thanks to Richard J. Wagner for the simulated annealing code

    Largely copied from generateCMMData
    """
    def __init__(self, minTemp=1.0e-9, nToPrint=0):
        """Construct a OrderPoints

        Inputs
        - minTemp: minimum temperature of the simulated annealer; the results seem fairly insensitive
            to this value as long as it is small enough
        - nToPrint: number of intermediate results from the simulated annealer to print to stdout
        """
        self.minTemp = float(minTemp)
        self.nToPrint = int(nToPrint)

    def computeEnergy(self, azAltList):
        """Compute the energy of a set of points.

        Distance between two adjacent points = max(abs(delta-az), abs(delta-alt))

        @param[in] azAltList  2-d array: (az, alt) x num points
        """
        absDiff = numpy.abs(numpy.diff(azAltList, axis=0))
        maxAbsDiff = numpy.max(absDiff, 1)
        totalCost = numpy.sum(maxAbsDiff)
        return totalCost

    def changeState(self, azAltList):
        """Change the state of azAltList in place.

        Randomly swap two points
        """
        nPts = len(azAltList)
        ind0 = random.randint(1, nPts-1)
        ind1 = random.randint(1, nPts-1)
        while ind1 == ind0:
            ind1 = random.randint(1, nPts-1)
        # make copy of the sources to make sure the swap works correctly
        azAltList[ind0], azAltList[ind1] = tuple(azAltList[ind1]), tuple(azAltList[ind0])

    def __call__(self, azAltList, nIter):
        """Reorder the points in azAltList to make a short path from the first to the last.

        Start at the point closest to az=90 (center of wrap), so the telescope
        moves equal amounts in both direction in azimuth.

        Inputs:
        - azAltList: a set of (x, y, ...?) data where x, y are positions on a plane
        - nIter: number of iterations
        """
        # numPoints = len(azAltList)
        begAzAltList = numpy.array(azAltList, dtype=float, copy=True)

        ctrAz = 90
        minInd = numpy.argmin(numpy.square(begAzAltList[:, 0] - ctrAz))
        if minInd != 0:
            begAzAltList[0], begAzAltList[minInd] = tuple(begAzAltList[minInd]), tuple(begAzAltList[0])
        if len(begAzAltList) > 2:
            initialEnergy = self.computeEnergy(begAzAltList)
            annealer = Annealer(self.computeEnergy, self.changeState)
            azAltList = annealer.anneal(begAzAltList, initialEnergy, self.minTemp, nIter, self.nToPrint)[0]
        else:
            azAltList = begAzAltList
        return azAltList


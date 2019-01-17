#!/usr/bin/env python2
from __future__ import division, absolute_import

import itertools
import Tkinter
import matplotlib
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy
import random

from anneal import Annealer

numpy.seterr(invalid='raise')
# numpy.random.seed(0)


class GenerateAzAltGrid(object):
    def __init__(self, nPts, minAlt=10, maxAlt=70, nIterDistribute=5000, nIterOrder=5000, azMix=60, fracMix=0.5):
        """Run a MCMC routine (specifically Metropolis-Hastings) to find
        A good distribution of points in alt, az space

        @param[in] nPts  number of points to find within the requested altitude bounds. Exact number not guarenteed in solution.
        @param[in] minAlt  degrees, minimum altitude to generate points for
        @param[in] maxAlt  degrees, max altitude to generate points for
        @param[in] nIterDistribute  number of iterations to run the MCMC algorithm for point distribution
        @param[in] nIterOrder  number of iterations to run the simulated annealer for order determination.
        @param[in] azMix  delta-az to mix at -90 boundary
        @param[in] fracMix  fraction of points to mix

        The results are available as two attributes:
        - azAltList: a list of (az, alt) points
        - unorderedAzAltList: a list of (az, alt) points before ordering is applied
        """
        # determine density based on nPts and altitude bounds.
        # solve the double integral- see area by integration: http://en.wikipedia.org/wiki/Sphere
        # Integrate(0, 2pi) Integrate(minPhi, maxPhi) sin theta dTheta dPhi
        # note: phi = [0, pi] radians, but max/min alt are input in degrees with +z = 90!
        # be careful! note r=1 for all calculations
        phi = numpy.pi/2. + numpy.radians([minAlt, maxAlt])
        # below is analytic solution to integral
        boundedArea = 2*numpy.pi*(numpy.cos(phi[0]) - numpy.cos(phi[1]))
        areaHemisphere = 2*numpy.pi
        # determine additional points required to maintain desired density on the hemisphere

        self.nPts = nPts * (areaHemisphere/boundedArea) * 1.08 # fudge found to be helpful via experimentation.
        self.maxAlt = float(maxAlt)
        self.minAlt = float(minAlt)
        self.nIterDistribute = int(nIterDistribute)
        self.nIterOrder = int(nIterOrder)
        self.azMix = float(azMix)
        self.fracMix = float(fracMix)
        self._orderPoints = OrderPoints(self.nIterOrder)

        # generate position vector point distribution
        xyzPoints = self.distributePoints(self.getInitialSample(10))
        # put in az, alt, throw out points outsize of altitude limits
        self.unorderedAzAltList = self.cart2sph(xyzPoints)
        # solve the order in az, alt coords
        if self.nIterOrder > 0:
            self.azAltList = self.orderPoints(self.unorderedAzAltList)
        else:
            self.azAltList = self.unorderedAzAltList

    def getInitialSample(self, nSamples=10):
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
                z = z*-1.
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
            alt = -1*(numpy.degrees(numpy.arccos(z))-90)
            if alt < self.maxAlt and alt > self.minAlt:
                # apply altitude bounds
                out.append([az, alt])
        return numpy.asarray(out)

    def distributePoints(self, posVecList):
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
        # start with a set of uniformly distributed points, over a sphere
        # put them in cartesian coordinate system
        currCoords = posVecList[:]
        currCost, worstIndex = self.computeDistributionCost(currCoords)
        bestCoords = currCoords[:] # will hold the best set of points seen
        bestCost = currCost
        for ii in xrange(self.nIterDistribute):
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
                    # print 'took worse state!!!!'
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
        return bestCoords

    def computeDistributionCost(self, posVecList):
        """Determine how good this configuration of coords is based on
        a cost.

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
        costMatrix = self.cosApproxCost(bigDot) # you may choose to use the arccosCost method instead
        # enter zero on the diagonal (no cost from being close to yourself)
        costMatrix[di] = 0.
        totalCost = numpy.sum(costMatrix)
        # now determine the index in posVecList with the highest cost
        flattenedCost = numpy.sum(costMatrix, 0)
        costlyPointIndex = numpy.argmax(flattenedCost)
        return totalCost, costlyPointIndex

    def cosApproxCost(self, bigDot):
        """
        Solve using metric based on cosine**2. 2x faster, less exact, but pretty good
        """
        distance = (1.000001 - bigDot**2)
        # cost function is (angular separation)^-2
        costMatrix = distance**-2
        return costMatrix

    def arccosCost(self, bigDot):
        """Solve using arccos, slower, more exact, computes actual cartesian distance betweeen points
        """
        # Contition the matrix. Some of the dot products are numerically > 1 or < -1. Enforce that that
        # cannot happen (it will if you don't). The arccos is only valid for values in [-1, 1]
        bigDot = bigDot / (numpy.max(numpy.abs(bigDot))+0.0000001)
        distance = numpy.arccos(bigDot)
        # cost function is (angular separation)^-2
        costMatrix = distance**-2
        return costMatrix

    def orderPoints(self, azAltList):
        """Break azAltList into a subsets of 3 sections, evenly spaced in azimuth.

        The idea is to get an order that is reasonably efficient, but not so nicely ordered
        that effects such as hysteresis are hidden

        The wrap limits are handled as follows:
        - The points are put in order -90, 270
        - Some points in the range azMix about -90 are put into the other wrap,
          so that points are now in the range -90 - azMix/2 to 270 + azMix/2
        - Those points are ordered, treating wrap as a large distance

        @param[in] azAltList  2D numpy array nPts x [az, alt]
        """
        azAltList = numpy.copy(azAltList)

        # wrap az into range -90, 270
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

        orderedCoords = self._orderPoints(azAltList)
        return orderedCoords


class OrderPoints(object):
    """Reorder a collection of points so they can be traversed quickly.

    To use:
    - Construct an OrderPoints
    - Call the object as a function with your list of points

    Credits:
    Thanks to Richard J. Wagner for the simulated annealing code

    Largely copied from generateCMMData
    """
    def __init__(self, nIter=500000, minTemp=1.0e-9, nToPrint=0):
        """Construct a OrderPoints

        Inputs
        - nIter: number of iterations of the simulated annealer
        - minTemp: minimum temperature of the simulated annealer; the results seem fairly insensitive
            to this value as long as it is small enough
        - nToPrint: number of intermediate results from the simulated annealer to print to stdout
        """
        self.nIter = int(nIter)
        self.minTemp = float(minTemp)
        self.nToPrint = int(nToPrint)

    def computeEnergy(self, azAltList):
        """Compute the energy of a set of points.

        Distance between two adjacent points = max(abs(delta-az), abs(delta-alt)),
        where delta-az is wrapped into the range (-180, 180]

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

    def __call__(self, azAltList):
        """Reorder the points in azAltList to make a short path from the first to the last.

        Start at the point closest to az=90 (center of wrap), so the telescope has
        unambiguous wrap for the first point.

        Inputs:
        - azAltList: a set of (x, y, ...?) data where x, y are positions on a plane
        - minToMax: if True, start from the smallest azimuth, else start at the largest
        """
        # numPoints = len(azAltList)
        begAzAltList = numpy.array(azAltList, dtype=float, copy=True)

        # find point closest to center of azimuth range and start with that
        ctrAz = 90
        minInd = numpy.argmin(numpy.square(begAzAltList[:, 0] - ctrAz))
        if minInd != 0:
            begAzAltList[0], begAzAltList[minInd] = tuple(begAzAltList[minInd]), tuple(begAzAltList[0])
        initialEnergy = self.computeEnergy(begAzAltList)
        annealer = Annealer(self.computeEnergy, self.changeState)
        azAltList = annealer.anneal(begAzAltList, initialEnergy, self.minTemp, self.nIter, self.nToPrint)[0]
        return azAltList


class AzOrderGraph(Tkinter.Frame):
    """Display Azmiuth against point order, to see what wrapping is looking like
    """
    def __init__(self, master):
        Tkinter.Frame.__init__(self, master)
        plotFig = matplotlib.figure.Figure(figsize=(6, 6), frameon=False)
        self.figCanvas = FigureCanvasTkAgg(plotFig, self)
        self.figCanvas.get_tk_widget().grid(row=0, column=0, sticky="news")
        self.grid_rowconfigure(0, weight=1)
        self.grid_columnconfigure(0, weight=1)
        self.axis = plotFig.add_subplot(1, 1, 1, autoscale_on=False)
        # self._setLimits()

    def _setLimits(self):
        """Update plot limits; must be done for every call to plotAzAltList
        """
        self.axis.set_xticklabels(['90', '135', '180', '-135', '-90', '-45', '0', '45'])
        self.axis.set_ylim(0, 90)
        self.axis.set_yticklabels([]) # ['80', '', '60', '', '40', '', '20', '', '0'])

    def plotAz(self, azAltList):
        """Plot az/alt points
        """
        self.axis.clear()

        # convert az, alt to r, theta, where r is 0 in the middle and theta is 0 right, 90 up
        az, alt = zip(*azAltList)
        # r = numpy.subtract(90, alt)
        # theta = numpy.deg2rad(numpy.subtract(az, 90))
        # az = numpy.subtract(az, 90)
        self.axis.plot(range(len(azAltList)), az, linestyle="solid", color="r", marker="o", markersize=4)
        self.axis.set_xlabel("Point In Order")
        self.axis.set_ylabel("Azimuth")
        # self._setLimits()
        self.figCanvas.draw()


class AzAltGraph(Tkinter.Frame):
    """Display points in an Az/Alt grid

    az 0 deg is down, 90 deg is right
    alt 90 deg is in the center
    """
    def __init__(self, master):
        Tkinter.Frame.__init__(self, master)
        plotFig = matplotlib.figure.Figure(figsize=(6, 6), frameon=False)
        self.figCanvas = FigureCanvasTkAgg(plotFig, self)
        self.figCanvas.get_tk_widget().grid(row=0, column=0, sticky="news")
        self.grid_rowconfigure(0, weight=1)
        self.grid_columnconfigure(0, weight=1)
        self.axis = plotFig.add_subplot(1, 1, 1, polar=True, autoscale_on=False)
        self._setLimits()

    def _setLimits(self):
        """Update plot limits; must be done for every call to plotAzAltList
        """
        self.axis.set_xticklabels(['90', '135', '180', '-135', '-90', '-45', '0', '45'])
        self.axis.set_ylim(0, 90)
        self.axis.set_yticklabels([]) # ['80', '', '60', '', '40', '', '20', '', '0'])

    def plotAzAltList(self, azAltList):
        """Plot az/alt points
        """
        self.axis.clear()

        # convert az, alt to r, theta, where r is 0 in the middle and theta is 0 right, 90 up
        az, alt = zip(*azAltList)
        r = numpy.subtract(90, alt)
        theta = numpy.deg2rad(numpy.subtract(az, 90))
        # theta = numpy.radians(az)
        self.axis.plot(theta, r, linestyle="none", color="r", marker="o", markersize=4)
        self._setLimits()
        self.figCanvas.draw()

    def plotAzAltListOrdered(self, azAltList):
        """Plot az/alt points, connect them with a line
        """
        # convert az, alt to r, theta, where r is 0 in the middle and theta is 0 right, 90 up
        az, alt = zip(*azAltList)
        r = numpy.subtract(90, alt)
        theta = numpy.deg2rad(numpy.subtract(az, 90))
        for r1, r2, t1, t2 in itertools.izip(r[:-1], r[1:], theta[:-1], theta[1:]):
            self.axis.plot([t1, t2], [r1, r2], linestyle="solid", color="r", marker="o", markersize=4)
        self._setLimits()
        self.figCanvas.draw()

    def plotErrors(self, azAltList, errPoints=()):
        """Plot error data

        Inputs:
        - azAltList: az,alt points at which errors were measured
        - errPoints: az,alt error at each az, alt point
        """
        if len(errPoints) != len(azAltList):
            raise RuntimeError("len(azAltList) = %s != %s = len(errPoints)" % (len(azAltList), len(errPoints)))

        # convert az, alt to r, theta, where r is 0 in the middle and theta is 0 right, 90 up
        az, alt = zip(*azAltList)
        r = numpy.subtract(90, alt)
        theta = numpy.deg2rad(numpy.subtract(az, 90))

        # quiver on a polar plot takes these strange arguments: theta, r, dx, dy
        sinAz = numpy.sin(numpy.deg2rad(az))
        cosAz = numpy.cos(numpy.deg2rad(az))
        azErr, altErr = zip(*errPoints)
        xErr = (cosAz * azErr) - (sinAz * altErr)
        yErr = (sinAz * azErr) + (cosAz * altErr)
        self.axis.quiver(theta, r, xErr, yErr, width=0.005)
        self.figCanvas.draw()

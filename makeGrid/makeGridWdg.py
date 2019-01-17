from __future__ import division, absolute_import

import Tkinter
import tkFileDialog

import matplotlib
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy
import RO.Wdg
import RO.Constants

from .generateGrid import GenerateAzAltGrid

class MakeGridWdg(Tkinter.Frame):
    def __init__(self, master):
        Tkinter.Frame.__init__(self, master)
        self.azAltList = []

        self.azAltGraph = AzAltGraph(master=self)
        self.azAltGraph.grid(row=0, column=0, sticky="nwes")

        entryFrame = Tkinter.Frame(self)
        gridder = RO.Wdg.Gridder(master=entryFrame)

        self.desNumPtsWdg = RO.Wdg.IntEntry(
            master = entryFrame,
            minValue = 1,
            defValue = 100,
            helpText = "number of points desired (approximate)",
        )
        gridder.gridWdg("Num Points", self.desNumPtsWdg)

        self.minAltWdg = RO.Wdg.FloatEntry(
            master = entryFrame,
            minValue = 0,
            maxValue = 90,
            defValue = 13,
            helpText = "minimum altitude",
        )
        gridder.gridWdg("Min Alt", self.minAltWdg, "deg")

        self.maxAltWdg = RO.Wdg.FloatEntry(
            master = entryFrame,
            minValue = 0,
            maxValue = 90,
            defValue = 82,
            helpText = "maximum altitude",
        )
        gridder.gridWdg("Max Alt", self.maxAltWdg, "deg")

        gridder.startNewCol()

        self.numIterDistributeWdg = RO.Wdg.IntEntry(
            master = entryFrame,
            minValue = 0,
            defValue = 5000,
            helpText = "number of iterations to distribute points",
        )
        gridder.gridWdg("Num Iter Distrib", self.numIterDistributeWdg)

        self.numIterOrderWdg = RO.Wdg.IntEntry(
            master = entryFrame,
            minValue = 0,
            defValue = 2500,
            helpText = "number of iterations to order points",
        )
        gridder.gridWdg("Num Iter Order", self.numIterOrderWdg)

        entryFrame.grid(row=1, column=0, sticky="w")

        self.statusBar = RO.Wdg.StatusBar(self)
        self.statusBar.grid(row=2, column=0, sticky="ew")

        ctrlFrame = Tkinter.Frame(self)
        self.generateWdg = RO.Wdg.Button(
            text = "Generate",
            master = ctrlFrame,
            command = self.doGenerate,
            helpText = "generate a new grid (overwriting the current grid)",
        )
        self.generateWdg.pack(side="left")

        self.saveWdg = RO.Wdg.Button(
            text = "Save",
            master = ctrlFrame,
            command = self.doSave,
            helpText = "save the current grid",
        )
        self.saveWdg.setEnable(False)
        self.saveWdg.pack(side="left")

        self.loadWdg = RO.Wdg.Button(
            text = "Load",
            master = ctrlFrame,
            command = self.doLoad,
            helpText = "load a grid from a file",
        )
        self.loadWdg.pack(side="left")
        ctrlFrame.grid(row=3, column=0, sticky="w")

        self.grid_rowconfigure(0, weight=1)
        self.grid_columnconfigure(0, weight=1)


    def doGenerate(self):
        """Generate a grid of points
        """
        nPts = self.desNumPtsWdg.getNum()
        nIterDistribute = self.numIterDistributeWdg.getNum()
        nIterOrder = self.numIterOrderWdg.getNum()
        minAlt = self.minAltWdg.getNum()
        maxAlt = self.maxAltWdg.getNum()

        self.statusBar.setMsg("Distributing points", isTemp=True)
        self.update_idletasks()

        generator = GenerateAzAltGrid(
            nPts = nPts,
            minAlt = minAlt,
            maxAlt = maxAlt,
        )
        totIter = 0
        while totIter < nIterDistribute:
            nIter = min(500, nIterDistribute - totIter)
            generator.distributePoints(nIter = nIter)
            totIter += nIter
            self.azAltList = generator.azAltList
            self.statusBar.setMsg("Distributing points: %d of %d iterations" % (totIter, nIterDistribute), isTemp=True)
            self.showPoints(showOrder=False)

        if len(self.azAltList) < 3:
            # just call once to set starting point
            nIterOrder = 1

        totIter = 0
        while totIter < nIterOrder:
            nIter = min(500, nIterOrder - totIter)
            generator.orderPoints(nIter = nIter)
            totIter += nIter
            self.azAltList = generator.azAltList
            self.statusBar.setMsg("Ordering points: %d of %d iterations" % (totIter, nIterOrder), isTemp=True)
            self.showPoints()
        self.statusBar.setMsg("Generated %s points" % (len(self.azAltList),), isTemp=True)

    def doSave(self):
        """Save a grid
        """
        if len(self.azAltList) < 1:
            self.saveWdg.setEnable(False)
            self.statusBar.setMsg("No grid to save", isTemp=True, severity=RO.Constants.sevError)
            return
        outPath = tkFileDialog.asksaveasfilename(
            defaultextension=".dat",
        )
        if not outPath:
            return
        with open(outPath, "w") as outFile:
            outFile.write("""! azimuth/altitude grid; # points = %s
 
! azimuth   altitude
!  (deg)      (deg)
""" % (len(self.azAltList),)
            )
            for azAltPoint in self.azAltList:
                outFile.write("%9.2f %9.2f\n" % (azAltPoint[0], azAltPoint[1]))
        self.statusBar.setMsg("Saved to %r" % (outPath,), isTemp=True)

    def doLoad(self):
        inPath = tkFileDialog.askopenfilename(
            defaultextension = ".dat",
        )
        if not inPath:
            return
        azAltList = []
        with open(inPath, "rU") as inFile:
            for i, line in enumerate(inFile):
                line = line.strip()
                if not line or line[0] in ("!", "#"):
                    continue
                dataList = line.split()
                try:
                    az, alt = [float(val) for val in dataList]
                except Exception:
                    self.statusBar.setMsg("Bad line %d: %r" % (i+1, line), isTemp=True, severity=RO.Constants.sevWarning)
                    continue
                azAltList.append((az, alt))
        if len(azAltList) < 1:
            self.statusBar.setMsg("No data found in %r" % (inPath,), isTemp=True, severity=RO.Constants.sevError)
            return

        self.azAltList = numpy.array(azAltList)
        self.statusBar.setMsg("Loaded %s points from %r" % (len(self.azAltList), inPath), isTemp=True)
        self.showPoints()

    def showPoints(self, showOrder=True):
        """Show the points on the graph

        @param[in] showOrder  show lines connecting the points?
        """
        if len(self.azAltList) > 0:
            self.saveWdg.setEnable(True)
            self.azAltGraph.plotAzAltList(self.azAltList, showOrder=showOrder)


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
        self.axis.set_xticklabels(['-90', '-135', 'N 180', '135', 'E 90', '45', '0', '-45'])
        self.axis.set_ylim(0, 90)
        self.axis.set_yticks((0, 15, 30, 45, 60, 75, 90))
        self.axis.set_yticklabels([]) # ['75', '60', '45', '15', '0'])

    def plotAzAltList(self, azAltList, showOrder):
        """Plot az/alt points

        @param[in]: azAltList: a collection of (az, alt) points (deg)
        @param[in] showOrder  show a line connecting the points in order?
        """
        self.axis.clear()

        # convert az, alt to r, theta, where r is 0 in the middle and theta is 0 left, 90 up
        az, alt = zip(*azAltList)
        r = numpy.subtract(90, alt)
        theta = numpy.deg2rad(numpy.subtract(270, az))
        # theta = numpy.radians(az)
        if showOrder:
            linestyle = "-"
        else:
            linestyle = ""
        self.axis.plot(theta, r, linestyle=linestyle, linewidth=0.5, color="gray", marker="o", markerfacecolor="black", markersize=4)
        if showOrder:
            self.axis.plot(theta[[0]], [r[0]], color="green", linestyle="", marker="^", markersize=8, label="start")
            self.axis.plot(theta[[-1]], [r[-1]], color="red", linestyle="", marker="v", markersize=8, label="end")
            self.axis.legend(loc="lower right", numpoints=1, bbox_to_anchor=(1.07, -0.12), frameon=False)
        self._setLimits()
        self.axis.text(-0.1, -0.08, "%s points" % (len(azAltList),), transform=self.axis.transAxes)
        self.figCanvas.draw()

#!/usr/bin/env python2
from __future__ import division #, absolute_import

import itertools
import Tkinter
import tkFileDialog
import matplotlib
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy
import RO.Wdg
import RO.Constants

from generateGrid import GenerateAzAltGrid

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

class MakeGridWdg(Tkinter.Frame):
    def __init__(self, master):
        Tkinter.Frame.__init__(self, master)
        self.azAltGraph = AzAltGraph(master=self)
        self.azAltGraph.grid(row=0, column=0, sticky="nwes")
        self.azAltList = []

        entryFrame = Tkinter.Frame(self)
        gridder = RO.Wdg.Gridder(master=entryFrame)

        self.numPtsWdg = RO.Wdg.IntEntry(
            master = entryFrame,
            minValue = 1,
            defValue = 200,
            helpText = "number of points desired (approximate)",
        )
        gridder.gridWdg("Num Points", self.numPtsWdg)

        self.minAltWdg = RO.Wdg.FloatEntry(
            master = entryFrame,
            minValue = 0,
            maxValue = 90,
            defValue = 8,
            helpText = "minimum altitude",
        )
        gridder.gridWdg("Min Alt", self.minAltWdg, "deg")

        self.maxAltWdg = RO.Wdg.FloatEntry(
            master = entryFrame,
            minValue = 0,
            maxValue = 90,
            defValue = 80,
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
        nPts = self.numPtsWdg.getNum()
        nIterDistribute = self.numIterDistributeWdg.getNum()
        nIterOrder = self.numIterOrderWdg.getNum()
        minAlt = self.minAltWdg.getNum()
        maxAlt = self.maxAltWdg.getNum()

        self.statusBar.setMsg("PLEASE WAIT: generating grid", isTemp=True)
        self.update_idletasks()

        generator = GenerateAzAltGrid(
            nPts = nPts,
            minAlt = minAlt,
            maxAlt = maxAlt,
            nIterDistribute = nIterDistribute,
            nIterOrder = nIterOrder,
        )
        self.azAltList = generator.azAltList
        if self.azAltList is None:
            raise RuntimeError("OOPS")
        self.showPoints()

    def doSave(self):
        """Save a grid
        """
        if len(self.azAltList) > 0:
            self.saveWdg.setEnable(False)
            self.statusBar.setMsg("No grid to save", isTemp=True, severity=RO.Constants.sevError)
            return
        outPath = tkFileDialog.asksaveasfile(
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
        self.statusBar.setMsg("Loaded %r points from" % (len(self.azAltList), inPath,), isTemp=True)
        self.showPoints()

    def showPoints(self):
        """Show the points on the graph
        """
        if len(self.azAltList) > 0:
            self.saveWdg.setEnable(True)
            self.azAltGraph.plotAzAltList(self.azAltList)


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
        self.axis.plot(theta, r, linestyle="solid", color="r", marker="o", markersize=4)
        self._setLimits()
        self.figCanvas.draw()

    def plotAzAltListOrdered(self, azAltList):
        """Plot az/alt points, connect them with a line
        """
        self.axis.clear()

        # convert az, alt to r, theta, where r is 0 in the middle and theta is 0 right, 90 up
        az, alt = zip(*azAltList)
        r = numpy.subtract(90, alt)
        theta = numpy.deg2rad(numpy.subtract(az, 90))
        for r1, r2, t1, t2 in itertools.izip(r[:-1], r[1:], theta[:-1], theta[1:]):
            self.axis.plot([t1, t2], [r1, r2], linestyle="", color="r", marker="o", markersize=4)
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


if __name__ == "__main__":
    root = Tkinter.Tk()
    wdg = MakeGridWdg(root)
    wdg.pack(expand=True, fill="both")
    root.mainloop()

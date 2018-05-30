# AstroCyles
version = 0.1
# Written by Tomass Wilson in python 3.6

# Import packages
import AstroLib as al
from astropy.stats import LombScargle
import sys
import linecache
import numpy as np
import wx
from pathlib import Path
from scipy import optimize as optimise
from scipy.signal import find_peaks_cwt
import math
# Matplotlib
import matplotlib
matplotlib.use('WXAgg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg
from matplotlib.figure import Figure


""" GLOBAL VARIABLES, CHANGE THESE TO FIT YOUR OBJECT """

# Coordinates of your object, used for HJD correction
# hours, minutes, seconds
rightAscension = '05 57 59.30'

# degrees, arcminutes, arcseconds
declination = '+53 53 44.9'

# JD constraints, most likely not needing changing, unless this program is still used in like the year 3000
# Minimum JD
minJD = 2400000

# Maximum JD
maxJD = 2500000

# Var constraints, should be changed to the ROUGH magnitudes of your object.
# This is mostly to make finding the output value, var, easier. If your data is good then this shouldnt be needed. Also changes how the program handles calib stars
# Minimum var
minVar = 5

# Maximum Var
maxVar = 20

# estimated periodicity
periodEst = 6

""" """

# Global variables DONT CHANGE
sizerBorder = 10
windowHeight = 800
windowWidth = 1400
fileList = []
borderLarge = 20
borderSmall = 5
linebreak = 50

# Code for HJD Correction
ra_hr, ra_min, ra_sec = rightAscension.split()
dec_deg, dec_min, dec_sec = declination.split()

obj = [ra_hr, ra_min, ra_sec, dec_deg, dec_min, dec_sec]


# Mostly just for luls. this class is a box for data
class data:

    def __init__(self, parent):
        self.reset()

    def reset(self):
        self.all = []
        self.JD = []
        self.HJD = []
        self.var = []
        self.error = []
        self.airmass = []


# Here we store the info for each file and modulate it
class fileConstruct:

    # Initiate yourself damnit
    def __init__(self, fileDir, exists=True):

        self.exists = exists

        # Initiate variables
        self.fileDir = fileDir
        self.format = al.getExtension(fileDir)
        self.extensionLength = len(self.format)

        self.identified = False
        self.oneLine = False
        self.observer = ""
        self.JDOffset = 0
        self.JDColumn = 0
        self.calibColumn = None
        self.varColumn = 1
        self.errorColumn = None
        self.airmassColumn = None
        self.totalColumns = 4

        self.SMAPower = 200
        self.foldingPeriod = 158.395

        self.SMAEnable = False
        self.sinusoid = False
        self.normalise = True
        self.detrend = False
        self.removeMacro = False
        self.foldingEnable = False
        self.showPeaks = False

        # files for storing fold informations, because it can be calculation intensive
        self.foldX = []
        self.foldY = []

        self.rawData = data(self)

        # Analyze yourself
        if self.exists:
            if self.format == ".astred":
                self.readAstred()
                self.gatherData()
            elif self.format == (".txt" or ".TXT"):
                self.analyze()
            else:
                raise TypeError("This is not a .txt or .astred file")

    # Function to read custom made ASTRED files
    def readAstred(self):
        myLines = self.getLines()

        self.observer = myLines[0][12:-1]
        self.JDOffset = float(myLines[1][17:])
        self.SMAPower = float(myLines[2][13:])
        self.SMAEnable = True if myLines[3][14:] == "True" else False
        self.detrend = True if myLines[4][11:] == "True" else False
        self.removeMacro = True if myLines[5][16:] == "True" else False
        self.oneLine = True if myLines[6][11:] == "True" else False
        self.errorColumn = 2
        self.identified = True

    # Return the lines for this file
    def getLines(self):
        myself = open(self.fileDir)
        myLines = list(myself)
        myself.close()
        return myLines

    # This is all of the original (Heavily modified) code to do with analyzing
    def analyze(self):

        # get our lines
        fileLines = self.getLines()
        # Iterate through each line
        for lineIndex, line in enumerate(fileLines):

            # This is the fun part

            # To begin, make sure we arent working with a new header, if so, reset format. Handles multiple files in one
            # NOT WORKING

            # Otherwise, if format is unidentified and this is not a header, identify format
            if self.identified is False and line[0] != "#" and len(line) >= 3:

                # step one, check if this is a lone line monstrosity. We do this by eing if there are more than 20 numbers in a line. If so then JD will have to be our marker for each datapoint
                numbersInLine = al.findNumbers(line)
                if len(numbersInLine) > 20:
                    self.oneLine = True

                # Create an array with possible JD on this line
                JDBox = al.findNumbers(line, True, minJD, maxJD)

                # If nothing is found then JD is likely formatted, find other values
                if len(JDBox) == 0:
                    JDBox = al.findNumbers(line, True)
                    # in this case we need to guess that JD is the first value
                    self.JDColumn = 0
                else:
                    # Lets note the column otherwise
                    self.JDColumn = JDBox[0][2]

                # Find JD Damnit, normally the first number. If it isnt, god save us
                JD = JDBox[0][0]
                print("JD of", JD, "found")

                # Analyze if JD is formatted
                if JD < minJD:
                    self.JDOffset = al.findNumOffset(
                        JD, fileLines, "JD", minJD, maxJD)

                varBox = al.findNumbers(line, True, minVar, maxVar)

                # If we find a var, it is either #1: Raw, in which case we need to get a calib column, or its in format #2, comp + (var-comp).
                if len(varBox) != 0:
                    var = varBox[0][0]
                    print("var of", var, "found")

                    # Save the column too, for later use
                    self.varColumn = varBox[0][2]

                    # Find the calib column, hopefully the second var-like number (eg betweeen minvar and maxvar). If it doesnt exist then we assume its already var-comp
                    if len(varBox) > 1:
                        self.calibColumn = varBox[1][2]

                # if we dont find a var though, its probably in format #3, var-comp, and we now have to guess that it's the second number
                else:
                    varBox = al.findNumbers(line, True)
                    self.varColumn = 1
                    var = varBox[1][0]
                    print("var of", var, "found")

                # Rinse and repeat for error and airmass, and hope to god they arent offset
                errorBox = al.findNumbers(line, False, 0, 0.1)
                if len(errorBox) != 0:
                    self.errorColumn = errorBox[0][2]

                airmassBox = al.findNumbers(line, False, 1, 2)
                if len(airmassBox) != 0:
                    self.airmassColumn = airmassBox[0][2]

                # We've done it, we worked out what the hell is going on in this file
                self.identified = True
                self.gatherData()

            elif line[0] == "#":
                sys.stdout.write(line)

    # A function to gather data from the txt file
    def gatherData(self):

        # Reset rawData
        self.rawData.reset()
        # set things to none if they should be None
        if self.errorColumn is None:
            self.rawData.error = None
        if self.airmassColumn is None:
            self.rawData.airmass = None
        # get our lines
        fileLines = self.getLines()
        # Iterate through each line
        for lineIndex, line in enumerate(fileLines):

            # Gather Data if we have a format and this is not a header
            if self.identified is True and line[0] != "#" and len(line) >= 3:

                # Initiate variables, for safety
                JD = 0
                var = 0
                error = None
                airmass = None

                # Collect a box of the numbers
                numberBox = al.findNumbers(line, True)

                # Different handling depending on file format
                if self.oneLine:
                    for numberIndex, number in enumerate(numberBox):

                        # IDK but this should work lol. Only works in oneline files where there are 4 columns. super convoluted but it worked once at least so yeah sorry
                        if (numberIndex + 4) % self.totalColumns == 0:

                            JD = numberBox[numberIndex + self.JDColumn][0] + self.JDOffset

                            # no handling for calib stars in this format, sorry
                            var = numberBox[numberIndex + self.varColumn][0]

                            if self.errorColumn is not None:
                                error = numberBox[numberIndex + self.errorColumn][0]
                            else:
                                error = None

                            if self.airmassColumn is not None:
                                airmass = numberBox[numberIndex + self.airmassColumn][0]
                            else:
                                airmass = None

                            # get HJD
                            HJD = al.trans(JD, obj)

                            # Save file data
                            self.rawData.all.append([JD, var, error, airmass])
                            self.rawData.JD.append(JD)
                            self.rawData.var.append(var)
                            self.rawData.error.append(error)
                            self.rawData.airmass.append(airmass)

                            # lol this should just work. HJD
                            self.rawData.HJD.append(HJD)

                else:
                    # Find JD
                    JD = numberBox[self.JDColumn][0] + self.JDOffset

                    # Find the variable
                    var = numberBox[self.varColumn][0]

                    # Remove calibstar if necessary
                    if var > minVar and var < maxVar and self.calibColumn is not None:
                        comp = numberBox[self.calibColumn][0]
                        var = comp + (var - comp)

                    # Rinse and repeat for error and airmass
                    if self.errorColumn is not None:
                        error = numberBox[self.errorColumn][0]
                        self.rawData.error.append(error)

                    if self.airmassColumn is not None:
                        airmass = numberBox[self.airmassColumn][0]
                        self.rawData.airmass.append(airmass)

                    # get HJD
                    HJD = al.trans(JD, obj)

                    # Save file data
                    self.rawData.all.append([JD, var, error, airmass])
                    self.rawData.JD.append(JD)
                    self.rawData.var.append(var)

                    # lol this should just work. HJD
                    self.rawData.HJD.append(HJD)

    # A function to generate a folded dataset
    def fold(self, period, bins=750, cycles=2):

        # period is in cycles/day so we need to convert it
        period = 1 / period

        # LOL check if the period is longer than the actuall dataset # FAIL
        if period > (self.rawData.HJD[-1] - self.rawData.HJD[1]):
            return self.rawData.HJD, self.rawData.var

        # Gather data
        xBox = self.rawData.HJD
        yBox = self.rawData.var

        # initiate new xBox and yBox
        foldedX = []
        foldedY = []

        # create a big array for sorting
        sortBox = []

        # Find the starting JD
        JDStart = xBox[0]

        # subtract the starting JD from each value (fold starts where it should and modulus calculations are faster)
        for xIndex, x in enumerate(xBox):
            xBox[xIndex] = x - JDStart

        # Iterate through every point and modulus it by the period to get its new X value
        for xIndex, x in enumerate(xBox):
            # add them to the sort box with corresponding Y values
            sortBox.append([x % period, yBox[xIndex]])

        # sort the sortbox
        sortBox.sort(key=al.sort_key)

        # split the values again, for simplicity
        for coord in sortBox:
            foldedX.append(coord[0])
            foldedY.append(coord[1])

        if bins is not None:
            # find the space between x values, which is the length of the fold devided by the number of points we want
            interval = period / bins

            # create our new x values
            xBins = []
            for i in range(bins):
                xBins.append(interval * i)

            # generate their corresponding y values
            binnedY = al.getSMA(foldedX, foldedY, 100000, bins=xBins)

            foldedX = xBins
            foldedY = binnedY

        if cycles is not None:
            multiplier = 1 / foldedX[-1]
            for xIndex, x in enumerate(foldedX):
                foldedX[xIndex] = x * multiplier

            npX = np.array(foldedX)
            npY = np.array(foldedY)

            i = 0
            while i < (cycles - 1):
                npX = np.append(npX, npX + i + 1)
                npY = np.append(npY, npY)
                i += 1

            foldedX = npX.tolist()
            foldedY = npY.tolist()

        # Save folded values for future use
        self.foldX = foldedX
        self.foldY = foldedY

        # return fold
        return foldedX, foldedY

    # a function to find peaks (not very good atm)
    def findPeaks(self, min, max):

        # this should be half width at half max (HWAHM), which is half the width of a gaussian curve at half its maximum when distributed around a peak !!!NOT USED RN!!!
        # width = 20

        # use findpeaks to find peaks in 1d arrays of the var values, inverted bc mag is inverted
        self.peakIndexes = find_peaks_cwt(np.array([x * -1 for x in self.rawData.var]), np.arange(min, max))

        # set showpeaks to true
        self.showPeaks = True

    # a function to do observed minus computed cycles (needs work)
    def OminusC(self, zeroPeriod):
        # O - C
        indexes = self.peakIndexes
        x = self.rawData.HJD
        # convert the array to anumpy array
        x = np.array(x)

        # use the first peak as zero
        zeroPeak = x[indexes[0]]

        # create our own period (that will be bad) if we arent given one
        zeroPeriod = x[indexes[1]] - zeroPeak if zeroPeriod is None else zeroPeriod

        # instantiate variables, cycle through peaks
        cycles = []
        observedMinusCalculated = []
        for i in range(len(indexes) - 1):

            # Get the next and current peaks
            currentPeak = x[indexes[i]]
            # nextPeak = x[indexes[i + 1]]

            # get the observed period
            # observed = currentPeak

            # calcualte the calculated cycle
            cycle = math.floor((currentPeak - zeroPeak) / zeroPeriod)
            # calculated = (cycle * zeroPeriod) + zeroPeak

            # dayDisparity = observed - calculated
            ######### this is essentially just a mod function, so
            dayDisparity = (currentPeak - zeroPeak) % zeroPeriod
            cycleDisparity = dayDisparity / zeroPeriod

            # idk títs wierd. this exists because it jumps around and stuff, but its not great
            if cycleDisparity > 0.5:
                cycleDisparity -= 1

            cycles.append(cycle)
            observedMinusCalculated.append(cycleDisparity)

        plt.scatter(cycles, observedMinusCalculated, color=(0, 0, 0.1, 0.75), marker="x")
        plt.show()


# Define File Drop Target class, no idea how this works
class FileDropTarget(wx.FileDropTarget):
    """ This object implements Drop Target functionality for Files """

    def __init__(self, obj):
        """ Initialize the Drop Target, passing in the Object Reference to
          indicate what should receive the dropped files """
        # Initialize the wxFileDropTarget Object
        wx.FileDropTarget.__init__(self)
        # Store the Object Reference for dropped files
        self.obj = obj

    def OnDropFiles(self, x, y, fileNames):
        """ Implement File Drop """
        # Here comes experimenta stuff
        """ Implement File Drop """
        # For Demo purposes, this function appends a list of the files dropped at the end of the widget's text
        # Move Insertion Point to the end of the widget's text
        # self.obj.SetInsertionPointEnd()
        # append a list of the file names dropped
        # self.obj.WriteText("%d file(s) dropped at %d, %d:\n" % (len(fileNames), x, y))
        # for file in fileNames:
        #     self.obj.WriteText(file + '\n')
        # self.obj.WriteText('\n')

        addFiles(fileNames)
        return len(fileNames)


# Detected values side Panel for the mainframe. Allows the user to fine tune the data collection from wierd TXT files
class DVPanel(wx.Panel):
    """This Panel is for changing detection mechanics"""

    def __init__(self, parent):
        """Create the panel"""
        wx.Panel.__init__(self, parent=parent, id=wx.ID_ANY)

        # create its sizer
        self.panelSizer = wx.BoxSizer(wx.VERTICAL)
        self.itemSizer = wx.BoxSizer(wx.VERTICAL)

        # Detected values header
        self.detectedValues = wx.StaticText(
            self, 0, "Detected values", style=wx.ALIGN_CENTRE_HORIZONTAL)
        self.detectedValues.SetFont(headerFont)
        # Add to right column
        self.itemSizer.Add(self.detectedValues, wx.SizerFlags(
            0).Border(wx.BOTTOM | wx.TOP, borderLarge))

        # JD gridsizer
        self.JDSizer = wx.GridSizer(
            rows=2, cols=2, hgap=borderSmall, vgap=borderSmall)

        # JD Column
        self.JDColumnText = wx.StaticText(self, 0, "JD column:")
        self.JDColumnText.SetFont(mainFont)
        self.JDColumn = wx.TextCtrl(self, 0, "")
        self.JDColumn.Disable()
        self.JDSizer.Add(self.JDColumnText, 0)
        self.JDSizer.Add(self.JDColumn, 0)

        # JD Offset
        self.JDOffsetText = wx.StaticText(self, 0, "JD offset:")
        self.JDOffsetText.SetFont(mainFont)
        self.JDOffset = wx.TextCtrl(self, 0, "")
        self.JDOffset.Disable()
        self.JDSizer.Add(self.JDOffsetText, 0)
        self.JDSizer.Add(self.JDOffset, 0)
        # Add to right column
        self.itemSizer.Add(self.JDSizer, wx.SizerFlags(
            0).Expand().Border(wx.BOTTOM, borderLarge))

        # var gridsizer
        self.varSizer = wx.GridSizer(
            rows=2, cols=2, hgap=borderSmall, vgap=borderSmall)

        # Var column
        self.varColumnText = wx.StaticText(self, 0, "Var column:")
        self.varColumnText.SetFont(mainFont)
        self.varColumn = wx.TextCtrl(self, 0, "")
        self.varColumn.Disable()
        self.varSizer.Add(self.varColumnText, 0)
        self.varSizer.Add(self.varColumn, 0)

        # calib star column
        self.calibColumnText = wx.StaticText(self, 0, "Calib column:")
        self.calibColumnText.SetFont(mainFont)
        self.calibColumn = wx.TextCtrl(self, 0, "")
        self.calibColumn.Disable()
        self.varSizer.Add(self.calibColumnText, 0)
        self.varSizer.Add(self.calibColumn, 0)
        # Add to right column
        self.itemSizer.Add(self.varSizer, wx.SizerFlags(
            0).Expand().Border(wx.BOTTOM, borderLarge))

        # error gridsizer
        self.errorSizer = wx.GridSizer(
            rows=1, cols=2, hgap=borderSmall, vgap=borderSmall)

        # error column
        self.errorColumnText = wx.StaticText(self, 0, "Error column:")
        self.errorColumnText.SetFont(mainFont)
        self.errorColumn = wx.TextCtrl(self, 0, "")
        self.errorColumn.Disable()
        self.errorSizer.Add(self.errorColumnText, 0)
        self.errorSizer.Add(self.errorColumn, 0)
        # Add to right column
        self.itemSizer.Add(self.errorSizer, wx.SizerFlags(
            0).Expand().Border(wx.BOTTOM, borderLarge))

        # airmass gridsizer
        self.airmassSizer = wx.GridSizer(
            rows=1, cols=2, hgap=borderSmall, vgap=borderSmall)

        # airmass column
        self.airmassColumnText = wx.StaticText(self, 0, "Airmass column:")
        self.airmassColumnText.SetFont(mainFont)
        self.airmassColumn = wx.TextCtrl(self, 0, "")
        self.airmassColumn.Disable()
        self.airmassSizer.Add(self.airmassColumnText, 0)
        self.airmassSizer.Add(self.airmassColumn, 0)
        # Add to right column
        self.itemSizer.Add(self.airmassSizer, wx.SizerFlags(
            0).Expand().Border(wx.BOTTOM, borderLarge))

        # Observer gridsizer
        self.observerSizer = wx.GridSizer(
            rows=1, cols=2, hgap=borderSmall, vgap=borderSmall)

        # Observer
        self.observerText = wx.StaticText(self, 0, "Observer:")
        self.observerText.SetFont(mainFont)
        self.observer = wx.TextCtrl(self, 0, "")
        self.observer.Disable()
        self.observerSizer.Add(self.observerText, 0)
        self.observerSizer.Add(self.observer, 0)
        # Add to right column
        self.itemSizer.Add(self.observerSizer, wx.SizerFlags(
            0).Expand().Border(wx.BOTTOM, linebreak))

        self.panelSizer.Add(self.itemSizer, wx.SizerFlags(
            0).Border(wx.LEFT, borderLarge))

        self.SetSizer(self.panelSizer)

        """ """

        # Create the detectionFields array
        self.detectionFields = [self.JDColumn, self.JDOffset, self.varColumn,
                                self.calibColumn, self.errorColumn, self.airmassColumn, self.observer]

        # Bindings for all the Detection fields
        for detectionIndex, detector in enumerate(self.detectionFields):
            # second time i ever use lambda, WTF
            detector.Bind(wx.EVT_TEXT, lambda evt,
                          temp=detectionIndex: self.onChangesToDetectionField(evt, temp))

    # detection field handler
    def onChangesToDetectionField(self, evt, detectionIndex):
        # A mental reminder for the fields:
        # 0 = self.JDColumn
        # 1 = self.JDOffset
        # 2 = self.varColumn
        # 3 = self.calibColumn
        # 4 = self.errorColumn
        # 5 = self.airmassColumn
        # 6 = self.observer

        # Modify the properties of this file
        if mainFrame.file is not None:
            shouldUpdate = False

            if detectionIndex == 0:
                source = al.getField(self.JDColumn, "int")
                if source is not None:
                    mainFrame.file.JDColumn = source
                    shouldUpdate = True
                else:
                    None

            elif detectionIndex == 1:
                source = al.getField(self.JDOffset, "float")
                if source is not None:
                    mainFrame.file.JDOffset = source
                    shouldUpdate = True
                else:
                    None

            elif detectionIndex == 2:
                source = al.getField(self.varColumn, "int")
                if source is not None:
                    mainFrame.file.varColumn = source
                    shouldUpdate = True
                else:
                    None

            elif detectionIndex == 3:
                source = al.getField(self.calibColumn, "int")
                if source is not None:
                    mainFrame.file.calibColumn = source
                    shouldUpdate = True
                else:
                    None

            elif detectionIndex == 4:
                source = al.getField(self.errorColumn, "int")
                if source is not None:
                    mainFrame.file.errorColumn = source
                    shouldUpdate = True
                else:
                    None

            elif detectionIndex == 5:
                source = al.getField(self.airmassColumn, "int")
                if source is not None:
                    mainFrame.file.airmassColumn = source
                    shouldUpdate = True
                else:
                    None

            elif detectionIndex == 6:
                mainFrame.file.observer = self.observer.GetValue()

            # Make sure the source has something valid in it and should update the plots n stuff
            if shouldUpdate:
                # update the data and display
                mainFrame.file.gatherData()
                mainFrame.updateFileInfo()

                # only update plot if we are plotting
                if mainFrame.plotOnMain:
                    mainFrame.updatePlot()


# Graph Modulation side Panel. Adds functions for a bunch of cool stuff like SMA, normalisation, folding, finding peaks and O - C
class GMPanel(wx.Panel):
    """This Panel is for changing how the graph is displayed/modulated"""

    def __init__(self, parent):
        """Create the panel"""
        wx.Panel.__init__(self, parent=parent, id=wx.ID_ANY)

        # create its sizer
        self.panelSizer = wx.BoxSizer(wx.VERTICAL)
        self.itemSizer = wx.BoxSizer(wx.VERTICAL)

        # Detected values header
        self.graphSettings = wx.StaticText(
            self, 0, "Graph Settings", style=wx.ALIGN_CENTRE_HORIZONTAL)
        self.graphSettings.SetFont(headerFont)
        # Add to right column
        self.itemSizer.Add(self.graphSettings, wx.SizerFlags(
            0).Border(wx.BOTTOM | wx.TOP, borderLarge))

        # SMA gridsizer
        self.SMASizer = wx.GridSizer(
            rows=2, cols=2, hgap=borderSmall, vgap=borderSmall)

        # SMA Enable
        self.SMAEnableText = wx.StaticText(self, 0, "Show SMA:")
        self.SMAEnableText.SetFont(mainFont)
        self.SMAEnable = wx.CheckBox(self, 0, "")
        self.SMAEnable.Disable()
        self.SMASizer.Add(self.SMAEnableText, 0)
        self.SMASizer.Add(self.SMAEnable, 0)

        # SMA Power
        self.SMAPowerText = wx.StaticText(self, 0, "SMA power:")
        self.SMAPowerText.SetFont(mainFont)
        self.SMAPower = wx.TextCtrl(self, 0, "")
        self.SMAPower.Disable()
        self.SMASizer.Add(self.SMAPowerText, 0)
        self.SMASizer.Add(self.SMAPower, 0)

        # Add to right column
        self.itemSizer.Add(self.SMASizer, wx.SizerFlags(
            0).Expand().Border(wx.BOTTOM, borderLarge))

        # Fit Sinusoid gridsizer
        self.sinusoidSizer = wx.GridSizer(
            rows=1, cols=2, hgap=borderSmall, vgap=borderSmall)

        # Fit Sinusoid
        self.sinusoidText = wx.StaticText(self, 0, "Fit Sinusoid:")
        self.sinusoidText.SetFont(mainFont)
        self.sinusoid = wx.CheckBox(self, 0, "")
        self.sinusoid.Disable()
        self.sinusoidSizer.Add(self.sinusoidText, 0)
        self.sinusoidSizer.Add(self.sinusoid, 0)

        # Add to right column
        self.itemSizer.Add(self.sinusoidSizer, wx.SizerFlags(
            0).Expand().Border(wx.BOTTOM, borderLarge))

        # Normalise gridsizer
        self.normaliseSizer = wx.GridSizer(
            rows=1, cols=2, hgap=borderSmall, vgap=borderSmall)

        # Normalise
        self.normaliseText = wx.StaticText(self, 0, "Normalise:")
        self.normaliseText.SetFont(mainFont)
        self.normalise = wx.CheckBox(self, 0, "")
        self.normalise.Disable()
        self.normaliseSizer.Add(self.normaliseText, 0)
        self.normaliseSizer.Add(self.normalise, 0)

        # Add to right column
        self.itemSizer.Add(self.normaliseSizer, wx.SizerFlags(
            0).Expand().Border(wx.BOTTOM, borderLarge))

        # Detrend gridsizer
        self.detrendSizer = wx.GridSizer(
            rows=1, cols=2, hgap=borderSmall, vgap=borderSmall)

        # Detrend
        self.detrendText = wx.StaticText(self, 0, "Detrend:")
        self.detrendText.SetFont(mainFont)
        self.detrend = wx.CheckBox(self, 0, "")
        self.detrend.Disable()
        self.detrendSizer.Add(self.detrendText, 0)
        self.detrendSizer.Add(self.detrend, 0)

        # Add to right column
        self.itemSizer.Add(self.detrendSizer, wx.SizerFlags(
            0).Expand().Border(wx.BOTTOM, borderLarge))

        # removeMacro gridsizer
        self.removeMacroSizer = wx.GridSizer(
            rows=1, cols=2, hgap=borderSmall, vgap=borderSmall)

        # Remove Macro
        self.removeMacroText = wx.StaticText(self, 0, "Remove Macro:")
        self.removeMacroText.SetFont(mainFont)
        self.removeMacro = wx.CheckBox(self, 0, "")
        self.removeMacro.Disable()
        self.removeMacroSizer.Add(self.removeMacroText, 0)
        self.removeMacroSizer.Add(self.removeMacro, 0)

        # Add to right column
        self.itemSizer.Add(self.removeMacroSizer, wx.SizerFlags(
            0).Expand().Border(wx.BOTTOM, borderLarge))

        # folding gridsizer
        self.foldingSizer = wx.GridSizer(
            rows=2, cols=2, hgap=borderSmall, vgap=borderSmall)

        # folding Enable
        self.foldingEnableText = wx.StaticText(self, 0, "Show fold:")
        self.foldingEnableText.SetFont(mainFont)
        self.foldingEnable = wx.CheckBox(self, 0, "")
        self.foldingEnable.Disable()
        self.foldingSizer.Add(self.foldingEnableText, 0)
        self.foldingSizer.Add(self.foldingEnable, 0)

        # folding Period
        self.foldingPeriodText = wx.StaticText(self, 0, "Folding period:")
        self.foldingPeriodText.SetFont(mainFont)
        self.foldingPeriod = wx.TextCtrl(self, 0, "")
        self.foldingPeriod.Disable()
        self.foldingSizer.Add(self.foldingPeriodText, 0)
        self.foldingSizer.Add(self.foldingPeriod, 0)

        # Add to right column
        self.itemSizer.Add(self.foldingSizer, wx.SizerFlags(
            0).Expand().Border(wx.BOTTOM, borderLarge))

        # functionButtons Gridsizer
        self.functionButtons = wx.GridSizer(
            rows=1, cols=2, hgap=borderSmall, vgap=borderSmall)

        # findPeaks button
        self.findPeaks = wx.Button(self, 0, "Find Peaks")
        self.findPeaks.Disable()
        self.functionButtons.Add(self.findPeaks, wx.SizerFlags(0).Align(wx.CENTRE))

        # OminusC button
        self.OminusC = wx.Button(self, 0, "O - C")
        self.OminusC.Disable()
        self.functionButtons.Add(self.OminusC, wx.SizerFlags(0).Align(wx.CENTRE))

        # Add to right column
        self.itemSizer.Add(self.functionButtons, wx.SizerFlags(
            0).Expand().Border(wx.BOTTOM, borderLarge))

        # finish up
        self.panelSizer.Add(self.itemSizer, wx.SizerFlags(
            0).Border(wx.LEFT, borderLarge))

        self.SetSizer(self.panelSizer)

        """ """

        # Create the detectionFields array
        self.detectionFields = [self.SMAPower, self.foldingPeriod]

        # Bindings for all the Detection fields
        for detectionIndex, detector in enumerate(self.detectionFields):
            # third time i ever use lambda, WTF
            detector.Bind(wx.EVT_TEXT, lambda evt, temp=detectionIndex: self.onChangesToDetectionField(evt, temp))

        # Create the Checkbox Array
        self.checkboxes = [self.SMAEnable, self.sinusoid, self.normalise, self.detrend, self.removeMacro, self.foldingEnable]

        # Bindings for checkboxes
        for checkboxIndex, checkbox in enumerate(self.checkboxes):
            # fourth time i ever use lambda, WTF
            checkbox.Bind(wx.EVT_CHECKBOX, lambda evt,
                          temp=checkboxIndex: self.onChecked(evt, temp))

        # create the buttons array
        self.buttons = [self.findPeaks, self.OminusC]

        # Bindings for buttons
        for checkboxIndex, checkbox in enumerate(self.buttons):
            # eighth time i ever use lambda, WTF
            checkbox.Bind(wx.EVT_BUTTON, lambda evt,
                          temp=checkboxIndex: self.onButton(evt, temp))

    # detection field handler
    def onChangesToDetectionField(self, evt, detectionIndex):
        # A mental reminder for the fields:
        # 0 = self.SMAPower
        # 1 = self.foldingPeriod

        # Modify the properties of this file
        if mainFrame.file is not None:
            shouldUpdate = False
            reFold = False

            if detectionIndex == 0:
                source = al.getField(self.SMAPower, "float")
                if source is not None:
                    mainFrame.file.SMAPower = source
                    shouldUpdate = True
                else:
                    None

            elif detectionIndex == 1:
                source = al.getField(self.foldingPeriod, "float")
                if source is not None:
                    mainFrame.file.foldingPeriod = source
                    shouldUpdate = True
                    reFold = True
                else:
                    None

            # Make sure the source has something valid in it and should update the plots n stuff
            if shouldUpdate and mainFrame.plotOnMain:
                mainFrame.updatePlot(reFold)

    # Checkbox Handler
    def onChecked(self, evt, checkboxIndex):
        # A mental reminder for the checkboxes:
        # 0 = self.SMAEnable
        # 1 = self.sinusoid
        # 2 = self.normalise
        # 3 = self.detrend
        # 4 = self.removeMacro
        # 5 = self.foldingEnable

        if mainFrame.file is not None:
            reFold = False

            if checkboxIndex == 0:
                if self.SMAEnable.IsChecked():
                    mainFrame.file.SMAEnable = True
                    self.SMAPower.Enable()
                else:
                    mainFrame.file.SMAEnable = False
                    self.SMAPower.Disable()

            elif checkboxIndex == 1:
                if self.sinusoid.IsChecked():
                    mainFrame.file.sinusoid = True
                else:
                    mainFrame.file.sinusoid = False

            elif checkboxIndex == 2:
                if self.normalise.IsChecked():
                    mainFrame.file.normalise = True
                else:
                    mainFrame.file.normalise = False

            elif checkboxIndex == 3:
                if self.detrend.IsChecked():
                    mainFrame.file.detrend = True
                else:
                    mainFrame.file.detrend = False

            elif checkboxIndex == 4:
                if self.removeMacro.IsChecked():
                    mainFrame.file.removeMacro = True
                else:
                    mainFrame.file.removeMacro = False

            elif checkboxIndex == 5:
                if self.foldingEnable.IsChecked():
                    mainFrame.file.foldingEnable = True
                    mainFrame.file.sinusoid = True
                    self.sinusoid.SetValue(True)
                    reFold = True
                else:
                    mainFrame.file.foldingEnable = False

            # Make sure we should update the plots n stuff
            if mainFrame.plotOnMain:
                mainFrame.updatePlot(reFold)

    # button handler
    def onButton(self, evt, detectionIndex):
        # A mental reminder for the fields:
        # 0 = self.findPeaks
        # 1 = self.OminusC

        # you know the drill
        if mainFrame.file is not None:
            shouldUpdate = False
            reFold = False

            if detectionIndex == 0:

                # either hide peaks or generate and show new peaks
                if mainFrame.file.showPeaks is True:
                    mainFrame.file.showPeaks = False
                    self.findPeaks.SetLabel("Find Peaks")
                    self.OminusC.Disable()
                    shouldUpdate = True
                else:
                    # create a dialog to get the settings for the peaks function
                    peaksDialog = functionDialog(parent=self, title="Find Peaks", fields=[{"title": "min", "type": "int"}, {"title": "max", "type": "int"}])
                    peaksDialog.ShowModal()
                    if mainFrame.dlgReturn is not None:
                        # get the data and reset the return
                        result = mainFrame.dlgReturn
                        mainFrame.dlgReturn = None
                        # something here should get min and max from the dialog
                        min = result[0]["result"]
                        max = result[1]["result"]

                        mainFrame.file.findPeaks(min, max)
                        self.findPeaks.SetLabel("Hide Peaks")
                        self.OminusC.Enable()
                        shouldUpdate = True

            elif detectionIndex == 1 and mainFrame.file.showPeaks is True:
                # create a dialog to get the settings for the O - c function
                OminusCDialog = functionDialog(parent=self, title="Observed Minus Calculated Period Graph", fields=[{"title": "Observed Period", "type": "float"}])
                OminusCDialog.ShowModal()
                if mainFrame.dlgReturn is not None:
                    # Get the data and reset the return
                    result = mainFrame.dlgReturn
                    mainFrame.dlgReturn = None
                    # Get observed from the dialog
                    observed = result[0]["result"]
                    # Run the graph
                    mainFrame.file.OminusC(observed)

            # Make sure the source has something valid in it and should update the plots n stuff
            if shouldUpdate and mainFrame.plotOnMain:
                mainFrame.updatePlot(reFold)


# A side panel for period graph settings. Not much here right now
class periodPanel(wx.Panel):
    """This side Panel is for changing how the periodogram is displayed"""

    def __init__(self, parent):
        """Create the panel"""
        wx.Panel.__init__(self, parent=parent, id=wx.ID_ANY)

        # create its sizer
        self.panelSizer = wx.BoxSizer(wx.VERTICAL)
        self.itemSizer = wx.BoxSizer(wx.VERTICAL)

        # Detected values header
        self.graphSettings = wx.StaticText(
            self, 0, "Graph Settings", style=wx.ALIGN_CENTRE_HORIZONTAL)
        self.graphSettings.SetFont(headerFont)
        # Add to right column
        self.itemSizer.Add(self.graphSettings, wx.SizerFlags(
            0).Border(wx.BOTTOM | wx.TOP, borderLarge))

        # minMax gridsizer, for the min and max x values of the periodogram
        self.minMaxSizer = wx.GridSizer(
            rows=2, cols=2, hgap=borderSmall, vgap=borderSmall)

        # minX
        self.minXFieldText = wx.StaticText(self, 0, "Min x:")
        self.minXFieldText.SetFont(mainFont)
        self.minXField = wx.TextCtrl(self, 0, "")
        self.minXField.ChangeValue(str(parent.minX))
        self.minMaxSizer.Add(self.minXFieldText, 0)
        self.minMaxSizer.Add(self.minXField, 0)

        # minY
        self.maxXFieldText = wx.StaticText(self, 0, "Max x:")
        self.maxXFieldText.SetFont(mainFont)
        self.maxXField = wx.TextCtrl(self, 0, "")
        self.maxXField.ChangeValue(str(parent.maxX))
        self.minMaxSizer.Add(self.maxXFieldText, 0)
        self.minMaxSizer.Add(self.maxXField, 0)

        # Add to right column
        self.itemSizer.Add(self.minMaxSizer, wx.SizerFlags(
            0).Expand().Border(wx.BOTTOM, borderLarge))

        self.panelSizer.Add(self.itemSizer, wx.SizerFlags(
            0).Border(wx.LEFT, borderLarge))

        self.SetSizer(self.panelSizer)

        """ """

        # Create the detectionFields array
        self.detectionFields = [self.minXField, self.maxXField]

        # Bindings for all the Detection fields
        for detectionIndex, detector in enumerate(self.detectionFields):
            # third time i ever use lambda, WTF
            detector.Bind(wx.EVT_TEXT, lambda evt, temp=detectionIndex: self.onChangesToDetectionField(parent, evt, temp))

        # Create the Checkbox Array
        self.checkboxes = None

        # Bindings for checkboxes
        if self.checkboxes is not None:
            for checkboxIndex, checkbox in enumerate(self.checkboxes):
                # fourth time i ever use lambda, WTF
                checkbox.Bind(wx.EVT_CHECKBOX, lambda evt,
                              temp=checkboxIndex: self.onChecked(evt, temp))

    # detection field handler
    def onChangesToDetectionField(self, parent, evt, detectionIndex):
        # A mental reminder for the fields:
        # 0 = self.minXField
        # 1 = self.maxXField

        # Modify the properties of this file
        if mainFrame.file is not None:
            shouldUpdate = False

            if detectionIndex == 0:
                source = al.getField(self.minXField, "int")
                if source is not None:
                    parent.minX = source
                    shouldUpdate = True
                else:
                    None

            elif detectionIndex == 1:
                source = al.getField(self.maxXField, "int")
                if source is not None:
                    parent.maxX = source
                    shouldUpdate = True
                else:
                    None

            # Make sure the source has something valid in it and should update the plots n stuff
            if shouldUpdate and mainFrame.plotOnMain:
                parent.updatePlot()

    # Checkbox Handler
    def onChecked(self, evt, checkboxIndex):
        # A mental reminder for the checkboxes:
        # THERE ARE NONE

        if mainFrame.file is not None:

            if checkboxIndex == 0:
                None


# a class for dialog windows when the user wants to use a complex algorythm. My first constructor yay!
class functionDialog(wx.Dialog):

    def __init__(self, parent, title, fields, id=-1):

        # get rows and cols
        if len(fields) > 6:
            self.cols = 6
        else:
            self.cols = 4
        # for every field we need 2 columns
        self.rows = math.ceil((len(fields) * 2) / self.cols)

        # set height and width
        height = 60 + (self.rows * 50)
        width = 500 if self.cols == 4 else 750
        size = (width, height)

        wx.Dialog.__init__(self, parent, id, title, size)

        # Fields come in form {"title": "something", "type": "int" or "float" or "check" or "string"}

        self.SetTitle(title)
        self.SetSize(size)
        self.fields = fields

        # create the supersizer
        self.superSizer = wx.BoxSizer(wx.VERTICAL)

        # create the field sizer
        self.fieldSizer = wx.GridSizer(
            rows=self.rows, cols=self.cols, hgap=borderSmall, vgap=borderSmall)

        # create each field with label
        for fieldIndex, field in enumerate(self.fields):
            fieldTitle = field["title"]
            fieldType = field["type"]
            ctrl = None

            # Instantiate the label
            text = wx.StaticText(self, 0, fieldTitle)
            text.SetFont(mainFont)
            self.fieldSizer.Add(text)

            # instantiate the correct control type
            if fieldType == "check":
                ctrl = wx.CheckBox(self, 0, "")
            else:
                ctrl = wx.TextCtrl(self, 0, "")

            self.fieldSizer.Add(ctrl)
            self.fields[fieldIndex]["ctrl"] = ctrl

        # Create the ok and close buttons
        self.buttonBox = wx.BoxSizer(wx.HORIZONTAL)
        self.okButton = wx.Button(self, label='Ok')
        self.closeButton = wx.Button(self, label='Close')
        self.buttonBox.Add(self.okButton, wx.SizerFlags(0).Align(wx.RIGHT))
        self.buttonBox.Add(self.closeButton, wx.SizerFlags(0).Border(wx.LEFT, borderSmall).Align(wx.RIGHT))

        # add the fields to the supersizer
        self.superSizer.Add(self.fieldSizer, wx.SizerFlags(0).Border(wx.ALL, borderSmall))
        self.superSizer.Add(self.buttonBox, wx.SizerFlags(0).Border(wx.LEFT | wx.RIGHT | wx.BOTTOM, borderSmall).Align(wx.RIGHT))
        self.SetSizer(self.superSizer)
        self.Layout()

        self.okButton.Bind(wx.EVT_BUTTON, self.onOK)
        self.closeButton.Bind(wx.EVT_BUTTON, self.onClose)

    def onOK(self, evt):
        # return field data
        for fieldIndex, field in enumerate(self.fields):
            fieldType = field["type"]
            ctrl = field["ctrl"]

            if fieldType == "check":
                self.fields[fieldIndex]["result"] = ctrl.GetValue()
            else:
                self.fields[fieldIndex]["result"] = al.getField(ctrl, fieldType)

        # IK its messy but RN this is what works
        mainFrame.dlgReturn = self.fields

        # destroy self
        self.Destroy()

    def onClose(self, evt):
        mainFrame.dlgReturn = None
        self.Destroy()


# MAinframe Right Column Notebook class
class RCNotebook(wx.Notebook):

    # ----------------------------------------------------------------------
    def __init__(self, parent):
        wx.Notebook.__init__(self, parent, id=wx.ID_ANY, style=wx.BK_DEFAULT)

        # Create the first tab and add it to the notebook
        self.DVPanel = DVPanel(self)
        self.DVPanel.SetBackgroundColour("light gray")
        self.AddPage(self.DVPanel, "Detection")

        # Create and add the second tab
        self.GMPanel = GMPanel(self)
        self.GMPanel.SetBackgroundColour("light gray")
        self.AddPage(self.GMPanel, "Graph")

        # self.Bind(wx.EVT_NOTEBOOK_PAGE_CHANGED, self.OnPageChanged)
        # self.Bind(wx.EVT_NOTEBOOK_PAGE_CHANGING, self.OnPageChanging)

    # def OnPageChanged(self, event):
    #     old = event.GetOldSelection()
    #     new = event.GetSelection()
    #     sel = self.GetSelection()
    #     event.Skip()

    # def OnPageChanging(self, event):
    #     old = event.GetOldSelection()
    #     new = event.GetSelection()
    #     sel = self.GetSelection()
    #     event.Skip()


# Gui for showing periodogram, is its own frame and has its own sidepanel plus gizmos like gaussian fits and monte carlo simulations
class periodFrame(wx.Frame):

    # Initiate
    def __init__(self):

        # initiate starting variables
        self.showGuass = False
        self.fidelity = 100000
        self.minX = 0
        self.maxX = 1000

        # some more initiation
        self.title = "Lomb Scargle Periodogram"
        wx.Frame.__init__(self, wx.GetApp().TopWindow, title=self.title, size=(1200, 800))
        wx.Frame.SetMinSize(self, (600, 400))
        self.CreateStatusBar()  # A Statusbar in the bottom of the window
        self.SetBackgroundColour("light gray")

        # grab the file from the mainwindow
        self.file = mainFrame.file

        ''' layout data '''

        # Create the outerSizer
        self.outerSizer = wx.BoxSizer(wx.VERTICAL)

        # Create the SuperSizer
        self.superSizer = wx.BoxSizer(wx.HORIZONTAL)

        # Create the left column
        self.leftColumn = wx.BoxSizer(wx.VERTICAL)

        #### GRAPH CANVAS

        # Create the matplotlib canvas
        self.figure = Figure()
        self.axes = self.figure.add_subplot(111)
        self.canvas = FigureCanvas(self, -1, self.figure)
        self.canvas.mpl_connect('motion_notify_event', self.on_move)

        # Add it to the left column
        self.leftColumn.Add(self.canvas, 1, wx.EXPAND)

        # Create the toolbar
        self.toolbar = NavigationToolbar2WxAgg(self.canvas)
        self.toolbar.Realize()

        # Add it to the bottom of the frame
        self.leftColumn.Add(self.toolbar, 0, wx.LEFT | wx.EXPAND)

        #### CONTROLS

        # create the control sizer
        self.controlSizer = wx.BoxSizer(wx.HORIZONTAL)

        # create the coordinate text
        self.xCoord = wx.StaticText(self, -1, "X:                   ")
        self.yCoord = wx.StaticText(self, -1, "Y:       ")

        # add them to the sizer
        self.controlSizer.Add(self.xCoord, wx.SizerFlags(0).Border(
            wx.LEFT, borderLarge).Align(wx.ALIGN_CENTRE_VERTICAL | wx.ALIGN_RIGHT))
        self.controlSizer.Add(self.yCoord, wx.SizerFlags(0).Border(
            wx.LEFT, borderLarge).Align(wx.ALIGN_CENTRE_VERTICAL | wx.ALIGN_RIGHT))

        # Add the controlsizer to the left column
        self.leftColumn.Add(self.controlSizer, wx.SizerFlags(0).Border(
            wx.TOP, borderSmall).Align(wx.ALIGN_CENTRE_HORIZONTAL))

        # Add the leftcolumn to the supersizer
        self.superSizer.Add(self.leftColumn, wx.SizerFlags(1).Expand())

        # Create the right Column
        self.rightColumn = wx.BoxSizer(wx.VERTICAL)

        #### RIGHT PANEL

        # Create the periodPanel
        self.periodPanel = periodPanel(self)

        # Add it to the right column
        self.rightColumn.Add(self.periodPanel, wx.SizerFlags(1).Expand())

        #### INFO

        # Info Header
        self.info = wx.StaticText(
            self, -1, "Info", style=wx.ALIGN_CENTRE_HORIZONTAL)
        self.info.SetFont(headerFont)
        # Add to right column
        self.rightColumn.Add(self.info, wx.SizerFlags(
            0).Expand().Border(wx.TOP, borderLarge))

        # info gridsizer
        self.infoGridSizer = wx.GridSizer(
            rows=3, cols=2, hgap=borderSmall, vgap=borderSmall)

        self.peakFreqLabelText = wx.StaticText(self, 0, "Peak Frequency:")
        self.peakFreqLabelText.SetFont(mainFont)
        self.peakFreqLabel = wx.StaticText(self, 0, "")
        self.peakFreqLabel.SetFont(mainFont)
        self.infoGridSizer.Add(self.peakFreqLabelText, 0)
        self.infoGridSizer.Add(self.peakFreqLabel, 0)

        self.sigmaLabelText = wx.StaticText(self, 0, "Standard Deviation:")
        self.sigmaLabelText.SetFont(mainFont)
        self.sigmaLabel = wx.StaticText(self, 0, "")
        self.sigmaLabel.SetFont(mainFont)
        self.infoGridSizer.Add(self.sigmaLabelText, 0)
        self.infoGridSizer.Add(self.sigmaLabel, 0)

        self.missingPeaksLabelText = wx.StaticText(self, 0, "Missed Peaks:")
        self.missingPeaksLabelText.SetFont(mainFont)
        self.missingPeaksLabel = wx.StaticText(self, 0, "")
        self.missingPeaksLabel.SetFont(mainFont)
        self.infoGridSizer.Add(self.missingPeaksLabelText, 0)
        self.infoGridSizer.Add(self.missingPeaksLabel, 0)

        # Add to right column
        self.rightColumn.Add(self.infoGridSizer, wx.SizerFlags(
            0).Expand().Border(wx.LEFT | wx.RIGHT, borderLarge))

        # Add the right column to the supersizer
        self.superSizer.Add(self.rightColumn, wx.SizerFlags(0))

        # Add the supersizer to the outerSizer
        self.outerSizer.Add(self.superSizer, wx.SizerFlags(
            1).Expand().Border(wx.ALL, borderLarge))

        #### TOOLS

        # Creater the toolsSizer
        self.toolSizer = wx.GridSizer(
            rows=1, cols=2, hgap=borderSmall, vgap=borderSmall)

        # Create the analyseSizer
        self.analyseSizer = wx.BoxSizer(wx.HORIZONTAL)

        # fitGuassButton
        self.fitGuassButton = wx.Button(self, -1, "Fit Guassian")
        self.analyseSizer.Add(self.fitGuassButton, wx.SizerFlags(0))

        # monteCarloButton
        self.monteCarloButton = wx.Button(self, -1, "Run Monte Carlo")
        self.analyseSizer.Add(self.monteCarloButton, wx.SizerFlags(0))

        # Add the analyseSizer to the toolbar
        self.toolSizer.Add(self.analyseSizer, wx.SizerFlags(0))

        # create the exportSizer
        self.exportSizer = wx.BoxSizer(wx.HORIZONTAL)

        # savePNGAllButton
        # self.savePNGAllButton = wx.Button(self, -1, "Save all plots")
        # self.savePNGAllButton.Disable()
        # self.exportSizer.Add(self.savePNGAllButton, wx.SizerFlags(0).Border(wx.RIGHT, borderSmall))

        # exportButton
        self.exportButton = wx.Button(self, -1, "Export")
        self.exportSizer.Add(self.exportButton, wx.SizerFlags(0))

        # Add the exportSizer to the toolbar
        self.toolSizer.Add(self.exportSizer, wx.SizerFlags(0).Align(wx.ALIGN_RIGHT))

        # Add the toolSizer to the outerSizer
        self.outerSizer.Add(self.toolSizer, wx.SizerFlags(0).Expand().Border(wx.LEFT | wx.BOTTOM | wx.RIGHT, borderLarge))

        ''' IT IS DONE '''
        ''' No more worries, the layout is behind us '''

        # Tool buttons. I realise this is a terrible way of manualy adding each tool to an array, but it works
        self.toolButtons = [self.fitGuassButton, self.monteCarloButton, self.exportButton]
        # Bind them toi the button handler
        for buttonIndex, button in enumerate(self.toolButtons):
            # sixth time i ever use lambda, WTF
            button.Bind(wx.EVT_BUTTON, lambda evt, temp=buttonIndex: self.onButton(evt, temp))

        # update the plot
        self.updatePlot()

        # Layout sizers
        self.SetSizer(self.outerSizer)
        self.SetAutoLayout(1)
        # self.superSizer.Fit(self)

    # Show coordinates for mouse cursor
    def on_move(self, event):
        # get the x and y pixel coords
        x, y = event.xdata, event.ydata

        if event.inaxes:
            self.xCoord.SetLabel("X: " + str(x))
            self.yCoord.SetLabel("Y: " + str(y))

    # Update Periodogram
    def updatePlot(self):
        # lock down the screen to stop awkward visuals
        self.Freeze()

        # clear plot
        self.axes.clear()

        time = self.file.rawData.HJD
        mag = self.file.rawData.var

        # Create a period graph with lomb scargle (using astropy)
        ls = LombScargle(time, mag)
        self.frequency, self.power = ls.autopower()

        # offset has been disabled in current instance. No periodogram stacking sorry :(

        self.figure.tight_layout()

        # Plot Time!!!
        # Plot a line and Xs for each datapoint
        self.axes.clear()
        self.axes.plot(self.frequency, self.power)

        # set the x limits
        limMin = 2
        limMax = 400
        self.axes.set_xlim(limMin, limMax)
        self.axes.set_ylim(0, 0.1)

        # show Guass!!
        if self.showGuass:
            self.axes.plot(self.gaussX, self.gaussY)

        # Error bars (not implemented properly)
        # self.axes.errorbar(x, yNorm, yerr=0.01, xerr=0.0001, fmt="ro")

        # Cosmetic things
        self.axes.set_xlabel('Cycles Per Day', fontsize=12)
        self.axes.set_ylabel('Significance', fontsize=12)

        # some significance lines, why not!
        probabilities = [0.1, 0.001, 0.00001]
        sigLines = ls.false_alarm_level(probabilities)
        self.axes.axhline(sigLines[0], label=str(int((1 - probabilities[0]) * 100)) + "%", color="r")
        self.axes.axhline(sigLines[1], label=str((1 - probabilities[1]) * 100) + "%", color="y")
        self.axes.axhline(sigLines[2], label=str(100 - (probabilities[2] * 100)) + "%", color="g")

        # lets make a period in minutes ticker along the top bc why not (BC ITS REALLY FUCKING HARD NVM THEN)
        # self.axes2 = self.axes.twiny()
        # # convert cycles/day into minutes/cycle
        # conversionFactor = (1 / np.array(self.frequency)) * 24 * 60
        # self.axes2.set_xlim(limMin*conversionFactor, limMax*conversionFactor)
        # self.axes2.set_xlabel("Minutes Per Cycle", fontsize=12)

        # Add a nice legend
        self.axes.legend()

        # update the toolbar
        self.toolbar.update()

        # finish up!
        self.Layout()
        self.Thaw()

    # Update the info stuffs
    def updateFileInfo(self, headerUpdate=False):

        self.peakFreqLabel.SetLabel(str(self.MCPeak))
        self.sigmaLabel.SetLabel(str(self.MCSD))
        self.missingPeaksLabel.SetLabel(str(self.missingPeaks))
        self.Layout()

    # button handler
    def onButton(self, evt, detectionIndex):
        # A mental reminder for the fields:
        # 0 = self.fitGuassButton
        # 1 = self.monteCarloButton
        # 2 = self.exportButton  ###### IDK WHAT TO DO HERE

        # you know the drill. run code for each button
        if mainFrame.file is not None:

            if detectionIndex == 0:
                self.fitGuassDialog()

            elif detectionIndex == 1:
                self.monteCarloDialog()

            elif detectionIndex == 2:
                None

    # function to create the gauss fit dialog and fit a gaussian using fitGuass
    def fitGuassDialog(self):

        # either hide the curve or generate and show a new curve
        if self.showGuass is True:
            self.showGuass = False
            self.fitGuassButton.SetLabel("Fit Guassian")
        else:
            # create a dialog to get the settings for the gauss function
            gaussDialog = functionDialog(parent=self, title="Fit Guassian", fields=[{"title": "Start", "type": "float"}, {"title": "Stop", "type": "float"}, {"title": "Center", "type": "float"}])
            gaussDialog.ShowModal()
            if mainFrame.dlgReturn is not None:
                # get the data and reset the return
                result = mainFrame.dlgReturn
                mainFrame.dlgReturn = None
                # get settings from the dialog
                start = result[0]["result"]
                stop = result[1]["result"]
                center = result[2]["result"]

                # get x and y
                bestValues, self.gaussX, self.gaussY = fitGauss(self.frequency, self.power, start, stop, center)
                self.showGuass = True
                self.fitGuassButton.SetLabel("Hide Guassian")

        # then update the plot
        self.updatePlot()

    # function to create the monte Carlo dialog and run lomb scargles on newly redistributed datapoints to nail down averages.
    def monteCarloDialog(self):
        # first, we need some input data, which we will then pass to the proper montecarlo function
        # create a dialog to get the range we will search for periods (start and stop) and the number of simulations to run
        MCDialog = functionDialog(parent=self, title="Run Periodicity Monte Carlo", fields=[{"title": "Minimum Frequency", "type": "float"}, {"title": "Maximum Frequency", "type": "float"}, {"title": "Runs", "type": "int"}])
        MCDialog.ShowModal()

        # if dialog returns
        if mainFrame.dlgReturn is not None:
            # get the data and reset the return
            result = mainFrame.dlgReturn
            mainFrame.dlgReturn = None
            # get settings from the dialog
            minFreq = result[0]["result"]
            maxFreq = result[1]["result"]
            runs = result[2]["result"]

            # get time and mag
            time = self.file.rawData.HJD
            mag = self.file.rawData.var
            error = self.file.rawData.error

            # run the simulation
            self.MCPeak, self.MCSD, self.missingPeaks = periodMC(time, mag, error, minFreq, maxFreq, runs)
            # update the button
            self.monteCarloButton.SetLabel("Rerun Monte Carlo")
            # update the info tab to display the new values
            self.updateFileInfo()

        # then update the plot
        self.updatePlot()


# GUI Data
class mainWindow(wx.Frame):
    """ We simply derive a new class of Frame. """

    def __init__(self, parent, title):
        self.file = None
        self.plotOnMain = False
        wx.Frame.__init__(self, parent, title=title,
                          size=(windowWidth, windowHeight))
        wx.Frame.SetMinSize(self, (1200, 800))
        self.CreateStatusBar()  # A Statusbar in the bottom of the window
        self.SetBackgroundColour("light gray")

        # Setting up the menu.
        filemenu = wx.Menu()

        # wx.ID_ABOUT and wx.ID_EXIT are standard IDs provided by wxWidgets.
        menuOpen = filemenu.Append(
            wx.ID_OPEN, "&Open", " Open a file to analyze")
        filemenu.AppendSeparator()
        menuItem = filemenu.Append(
            wx.ID_ABOUT, "&About", " Information about astroRedux")
        filemenu.AppendSeparator()
        menuExit = filemenu.Append(
            wx.ID_EXIT, "E&xit", " Terminate the program")

        # Creating the menubar.
        menuBar = wx.MenuBar()
        # Adding the "filemenu" to the MenuBar
        menuBar.Append(filemenu, "&File")
        self.SetMenuBar(menuBar)  # Adding the MenuBar to the Frame content.

        # Events
        # Menu Items
        self.Bind(wx.EVT_MENU, self.OnAbout, menuItem)
        self.Bind(wx.EVT_MENU, self.OnExit, menuExit)
        self.Bind(wx.EVT_MENU, self.OnOpen, menuOpen)
        # X TO CLOSE ALL
        self.Bind(wx.EVT_CLOSE, self.OnExit)

        ''' HERE THERE BE DRAGONS '''
        ''' beyond this point is the layout data, tread with care '''

        # Create the outerSizer
        self.outerSizer = wx.BoxSizer(wx.VERTICAL)

        # Create the SuperSizer
        self.superSizer = wx.BoxSizer(wx.HORIZONTAL)

        # Create the left column
        self.leftColumn = wx.BoxSizer(wx.VERTICAL)

        #### HEADER TEXT BOX

        # Define a Text Control to receive Dropped Files
        # Create a read-only Text Control
        self.mainBox = wx.TextCtrl(self, -1, "Drag a .txt file here to begin processing",
                                   style=wx.TE_MULTILINE | wx.HSCROLL | wx.TE_READONLY)
        # Make this control a File Drop Target
        # Create a File Drop Target object
        dt3 = FileDropTarget(self.mainBox)
        # Link the Drop Target Object to the Text Control
        self.mainBox.SetDropTarget(dt3)
        # Add it to the left column
        self.leftColumn.Add(self.mainBox, 2, wx.EXPAND)

        #### GRAPH CANVAS

        # Create the matplotlib canvas
        self.figure = Figure()
        self.axes = self.figure.add_subplot(111)
        self.canvas = FigureCanvas(self, -1, self.figure)
        self.canvas.mpl_connect('motion_notify_event', self.on_move)

        # Add it to the left column
        self.leftColumn.Add(self.canvas, 1, wx.EXPAND)

        # Hide the canvas for now
        self.canvas.Hide()

        # Create the toolbar
        self.toolbar = NavigationToolbar2WxAgg(self.canvas)
        self.toolbar.Realize()

        # Add it to the bottom of the frame
        self.leftColumn.Add(self.toolbar, 0, wx.LEFT | wx.EXPAND)

        # Hide the toolbar for now
        self.toolbar.Hide()

        #### CONTROLS

        # create the control sizer
        self.controlSizer = wx.BoxSizer(wx.HORIZONTAL)

        # create directional buttons sizer
        self.directionButtonSizer = wx.BoxSizer(wx.HORIZONTAL)

        # Buttons, to change file/header
        self.directionButtons = []
        self.directionButtons.append(wx.Button(self, -1, "<<<"))
        self.directionButtonSizer.Add(self.directionButtons[0], 1, wx.EXPAND)
        self.directionButtons.append(wx.Button(self, -1, "<"))
        self.directionButtonSizer.Add(self.directionButtons[1], 1, wx.EXPAND)
        self.fileNumber = wx.TextCtrl(self, -1, "", style=wx.TE_CENTRE)
        self.directionButtonSizer.Add(self.fileNumber, 1, wx.EXPAND)
        self.directionButtons.append(wx.Button(self, -1, ">"))
        self.directionButtonSizer.Add(self.directionButtons[2], 1, wx.EXPAND)
        self.directionButtons.append(wx.Button(self, -1, ">>>"))
        self.directionButtonSizer.Add(self.directionButtons[3], 1, wx.EXPAND)

        # add buttons to the controlsizer
        self.controlSizer.Add(self.directionButtonSizer, wx.SizerFlags(
            0).Align(wx.ALIGN_CENTRE_HORIZONTAL))

        # create the coordinate text
        self.xCoord = wx.StaticText(self, -1, "X:                   ")
        self.yCoord = wx.StaticText(self, -1, "Y:       ")

        # add them to the sizer
        self.controlSizer.Add(self.xCoord, wx.SizerFlags(0).Border(
            wx.LEFT, borderLarge).Align(wx.ALIGN_CENTRE_VERTICAL | wx.ALIGN_RIGHT))
        self.controlSizer.Add(self.yCoord, wx.SizerFlags(0).Border(
            wx.LEFT, borderLarge).Align(wx.ALIGN_CENTRE_VERTICAL | wx.ALIGN_RIGHT))

        # Add the controlsizer to the left column
        self.leftColumn.Add(self.controlSizer, wx.SizerFlags(0).Border(
            wx.TOP, borderSmall).Align(wx.ALIGN_CENTRE_HORIZONTAL))

        # Add the leftcolumn to the supersizer
        self.superSizer.Add(self.leftColumn, wx.SizerFlags(1).Expand())

        # Create the right Column
        self.rightColumn = wx.BoxSizer(wx.VERTICAL)

        #### NOTEBOOK

        # Create the notebook
        self.notebook = RCNotebook(self)

        # Add it to the right column
        self.rightColumn.Add(self.notebook, wx.SizerFlags(1).Expand())

        #### INFO

        # Info Header
        self.info = wx.StaticText(
            self, -1, "Info", style=wx.ALIGN_CENTRE_HORIZONTAL)
        self.info.SetFont(headerFont)
        # Add to right column
        self.rightColumn.Add(self.info, wx.SizerFlags(
            0).Expand().Border(wx.TOP, borderLarge))

        # info gridsizer
        self.infoGridSizer = wx.GridSizer(
            rows=3, cols=2, hgap=borderSmall, vgap=borderSmall)

        self.linesOfDataText = wx.StaticText(self, 0, "Lines of data:")
        self.linesOfDataText.SetFont(mainFont)
        self.linesOfData = wx.StaticText(self, 0, "")
        self.linesOfData.SetFont(mainFont)
        self.infoGridSizer.Add(self.linesOfDataText, 0)
        self.infoGridSizer.Add(self.linesOfData, 0)

        self.startDateText = wx.StaticText(self, 0, "Start date:")
        self.startDateText.SetFont(mainFont)
        self.startDate = wx.StaticText(self, 0, "")
        self.startDate.SetFont(mainFont)
        self.infoGridSizer.Add(self.startDateText, 0)
        self.infoGridSizer.Add(self.startDate, 0)

        self.endDateText = wx.StaticText(self, 0, "End date:")
        self.endDateText.SetFont(mainFont)
        self.endDate = wx.StaticText(self, 0, "")
        self.endDate.SetFont(mainFont)
        self.infoGridSizer.Add(self.endDateText, 0)
        self.infoGridSizer.Add(self.endDate, 0)

        # Add to right column
        self.rightColumn.Add(self.infoGridSizer, wx.SizerFlags(
            0).Expand().Border(wx.LEFT | wx.RIGHT, borderLarge))

        # Add the right column to the supersizer
        self.superSizer.Add(self.rightColumn, wx.SizerFlags(0))

        # Add the supersizer to the outerSizer
        self.outerSizer.Add(self.superSizer, wx.SizerFlags(
            1).Expand().Border(wx.ALL, borderLarge))

        #### TOOLS

        # Creater the toolsSizer
        self.toolSizer = wx.GridSizer(
            rows=1, cols=2, hgap=borderSmall, vgap=borderSmall)

        # Create the analyseSizer
        self.analyseSizer = wx.BoxSizer(wx.HORIZONTAL)

        # plotButton
        self.plotButton = wx.Button(self, -1, "Plot Graph")
        self.plotButton.Disable()
        self.analyseSizer.Add(self.plotButton, wx.SizerFlags(0).Border(wx.RIGHT, borderSmall))

        # periodButton
        self.periodButton = wx.Button(self, -1, "Analyse Period")
        self.periodButton.Disable()
        self.analyseSizer.Add(self.periodButton, wx.SizerFlags(0).Border(wx.RIGHT, borderSmall))

        # combineButton
        self.combineButton = wx.Button(self, -1, "Combine Data")
        self.combineButton.Disable()
        self.analyseSizer.Add(self.combineButton, wx.SizerFlags(0))

        # Add the analyseSizer to the toolbar
        self.toolSizer.Add(self.analyseSizer, wx.SizerFlags(0))

        # create the exportSizer
        self.exportSizer = wx.BoxSizer(wx.HORIZONTAL)

        # savePNGAllButton
        # self.savePNGAllButton = wx.Button(self, -1, "Save all plots")
        # self.savePNGAllButton.Disable()
        # self.exportSizer.Add(self.savePNGAllButton, wx.SizerFlags(0).Border(wx.RIGHT, borderSmall))

        # savePNGButton
        self.savePNGButton = wx.Button(self, -1, "Save Plot")
        self.savePNGButton.Disable()
        self.exportSizer.Add(self.savePNGButton, wx.SizerFlags(0).Border(wx.RIGHT, borderSmall))

        # exportAllButton
        self.exportAllButton = wx.Button(self, -1, "Export All")
        self.exportAllButton.Disable()
        self.exportSizer.Add(self.exportAllButton, wx.SizerFlags(0).Border(wx.RIGHT, borderSmall))

        # exportButton
        self.exportButton = wx.Button(self, -1, "Export")
        self.exportButton.Disable()
        self.exportSizer.Add(self.exportButton, wx.SizerFlags(0))

        # Add the exportSizer to the toolbar
        self.toolSizer.Add(self.exportSizer, wx.SizerFlags(0).Align(wx.ALIGN_RIGHT))

        # Add the toolSizer to the outerSizer
        self.outerSizer.Add(self.toolSizer, wx.SizerFlags(0).Expand().Border(wx.LEFT | wx.BOTTOM | wx.RIGHT, borderLarge))

        ''' IT IS DONE '''
        ''' No more worries, the layout is behind us '''

        # Bindings for buttons and stuff
        for buttonIndex, button in enumerate(self.directionButtons):
            # First time i ever use lambda, WTF
            button.Bind(wx.EVT_BUTTON, lambda evt, temp=buttonIndex: self.onDirectionButton(evt, temp))

        # Tool buttons
        self.toolButtons = [self.plotButton, self.periodButton, self.combineButton, self.savePNGButton, self.exportAllButton, self.exportButton]
        # Bind them
        for buttonIndex, button in enumerate(self.toolButtons):
            # Fourth time i ever use lambda, WTF
            button.Bind(wx.EVT_BUTTON, lambda evt, temp=buttonIndex: self.onToolbarButton(evt, temp))

        # fileNumber changes detector
        self.fileNumber.Bind(wx.EVT_TEXT, self.onChangeFile)

        # Layout sizers
        self.SetSizer(self.outerSizer)
        self.SetAutoLayout(1)
        # self.superSizer.Fit(self)

        self.Show(True)

    # Time to brag
    def OnAbout(self, event):
        # A message dialog box with an OK button. wx.OK is a standard ID in wxWidgets.
        dlg = wx.MessageDialog(self, "AstroCycles version " + str(version) + "\nA tool for making my life easier.\nDeveloped by Tomass Wilson in Python 3.6.\n(C) 2017, credit to everyone behind astropy, wxPython, matplotlib, numpy and of course the python team.\nSpecial thanks to Helena Uthas", "About AstroCyles", style=wx.OK | wx.CENTRE)
        dlg.ShowModal()  # Show it
        dlg.Destroy()  # finally destroy it when finished.

    # Oof, someone used the ancient "open" command, what a nerd
    def OnOpen(self, event):
        """ Open a file"""
        dirname = ''
        dlg = wx.FileDialog(self, "Choose a file",
                            self.dirname, "", "*.*", wx.FD_OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            dirname = dlg.GetDirectory()
            addFiles(dirname)
        dlg.Destroy()

    # Direction button handler
    def onDirectionButton(self, event, buttonIndex):
        # A mental reminder for the buttons:
        # 0 = <<<
        # 1 = <
        # 2 = >
        # 3 = >>>

        currentFile = al.getField(self.fileNumber, "int")

        # What the buttons do, including a check to see if we have any files
        # Also make sure we are not out of bounds
        if len(fileList) != 0:

            if buttonIndex == 0:
                self.fileNumber.SetValue("1")
            elif buttonIndex == 1 and currentFile >= 2:
                self.fileNumber.SetValue(
                    str(al.getField(self.fileNumber, "int") - 1))
            elif buttonIndex == 2 and currentFile < len(fileList):
                self.fileNumber.SetValue(
                    str(al.getField(self.fileNumber, "int") + 1))
            elif buttonIndex == 3:
                self.fileNumber.SetValue(str(len(fileList)))

    # Toolbar button handler
    def onToolbarButton(self, event, buttonIndex):
        # A mental reminder for the buttons:
        # 0 = plotButton
        # 1 = periodButton
        # 2 = combineButton
        #### 3 = savePNGAllButton
        # 3 = savePNGButton
        # 4 = exportAllButton
        # 5 = exportButton

        # What the buttons do
        if buttonIndex == 0:
            self.onPlot()
        elif buttonIndex == 1:
            self.periodAnalyse()
        elif buttonIndex == 2:
            self.combineCurves()
        # elif buttonIndex == 3:
        #     self.onSave(True)
        elif buttonIndex == 3:
            self.onSave(False)
        elif buttonIndex == 4:
            self.onExport(True)
        elif buttonIndex == 5:
            self.onExport(False)

    # DATA! YAY
    def onPlot(self):

        # which button function do we want
        if not self.plotOnMain and self.file is not None:

            self.plotOnMain = True

            # update button
            self.plotButton.SetLabel("Show Header")

            # swap frames
            self.mainBox.Hide()
            self.canvas.Show()
            self.toolbar.Show()

            self.updatePlot(reFold=True)

        # Return to the header
        elif self.plotOnMain:

            self.plotOnMain = False

            self.plotButton.SetLabel("Plot Graph")
            self.canvas.Hide()
            self.toolbar.Hide()
            self.mainBox.Show()
            self.updateHeader()
            self.Layout()

    # Update the graph to show the latest data
    def updatePlot(self, reFold=False):
        # lock down the screen to stop awkward visuals
        self.Freeze()

        # clear plot
        self.axes.clear()

        # get data, depending on if we are folding or not
        if self.file.foldingEnable:
            # save on processing if we dont need to refold the entire thing
            if reFold:
                x, y = self.file.fold(self.file.foldingPeriod, bins=None)
            else:
                x = self.file.foldX
                y = self.file.foldY
        else:
            x = self.file.rawData.HJD
            y = self.file.rawData.var

        # modify data, dont bother normalising if we are already detrending
        if self.file.detrend:
            y = al.removeTrend(x, y)
        elif self.file.normalise:
            y = al.normalise(x, y)

        # Remove macro if necessary
        if self.file.removeMacro:
            y = al.removeMacro(x, y)

        # Plot Time!!!
        # Plot a line and Xs for each datapoint
        self.axes.clear()
        self.axes.plot(x, y, color=(0, 0, 0.1, 0.35), lineStyle=":")
        self.axes.scatter(x, y, color=(0, 0, 0.1, 0.35), marker="x", label="Data Points")

        # Plot SMA
        if self.file.SMAEnable:
            SMAY = al.getSMA(x, y, self.file.SMAPower)
            self.axes.plot(x, SMAY, "r-", label="SMA")

        # print peaks that we have detected
        if self.file.showPeaks:
            self.axes.plot(np.array(x)[self.file.peakIndexes], np.array(y)[self.file.peakIndexes], "*", ms=20, color="green", label="Peaks")

        # plot sinusoid
        if self.file.sinusoid:
            y = np.array(al.getSMA(x, y, 50))
            x = np.array(x)
            try:
                parameters, parametersCovariance, = optimise.curve_fit(sinusoid, x, y, p0=[0, periodEst, 0, 0, periodEst * 2, 0])
                self.axes.plot(x, sinusoid(x, parameters[0], parameters[1], parameters[2], parameters[3], parameters[4], parameters[5]), "b-", label="Fitted Sinusoid")
            except RuntimeError:
                wx.MessageBox('Sinusoid could not be fit!', 'Warning', wx.OK | wx.ICON_WARNING)
                self.file.sinusoid = False
                self.updateFileInfo()

        # Error bars
        if not self.file.foldingEnable:
            self.axes.errorbar(x, y, yerr=self.file.rawData.error, color=(0, 0, 0.1, 0.05))

        # Cosmetic things
        if self.file.foldingEnable:
            self.axes.set_xlabel('Cycles', fontsize=18)
        else:
            self.axes.set_xlabel('HJD', fontsize=18)

        if self.file.normalise or self.file.detrend or self.file.removeMacro:
            self.axes.set_ylabel('Normalised Magnitude', fontsize=18)
        else:
            self.axes.set_ylabel('Magnitude', fontsize=18)

        # Invert the x value because magnitude is brighter for smaller numbers
        self.axes.invert_yaxis()

        # add the legend
        self.axes.legend()

        # redraw the canvas
        self.canvas.draw()

        # update the toolbar
        self.toolbar.update()

        # finish up!
        self.Thaw()
        self.Layout()

    # Period graph function
    def periodAnalyse(self):
        periodFrame().Show()

    # Show coordinates for mouse cursor
    def on_move(self, event):
        # get the x and y pixel coords
        x, y = event.xdata, event.ydata

        if event.inaxes:
            self.xCoord.SetLabel("X: " + str(x))
            self.yCoord.SetLabel("Y: " + str(y))

    # Update the displayed header and whatnot, when the selected file is changed
    def onChangeFile(self, event):
        fileToGoTo = al.getField(self.fileNumber, "int")

        # Check if we are out of bounds
        if len(fileList) == 0:
            self.fileNumber.ChangeValue("")
        elif fileToGoTo > len(fileList):
            self.fileNumber.ChangeValue(str(len(fileList)))
        elif fileToGoTo < 0:
            self.fileNumber.ChangeValue("1")
        else:
            # This is the start of a huge process
            # 1. Get the file
            self.file = fileList[fileToGoTo - 1]

            # 2. depends on what is onscreen
            if self.plotOnMain:
                self.updatePlot(reFold=True)
                self.updateFileInfo(False)
            else:
                self.updateFileInfo(True)

    # Update the info stuffs
    def updateFileInfo(self, headerUpdate=False):

        # Only do this if we want to update the header, i.e the headerbox is on the main page
        if headerUpdate:
            self.updateHeader()

        # update the JD offset
        self.notebook.DVPanel.JDOffset.ChangeValue(str(self.file.JDOffset))

        # update the var offset
        self.notebook.DVPanel.calibColumn.ChangeValue(
            str(self.file.calibColumn))

        # update other things
        self.notebook.DVPanel.JDColumn.ChangeValue(str(self.file.JDColumn))
        self.notebook.DVPanel.varColumn.ChangeValue(
            str(self.file.varColumn))
        self.notebook.DVPanel.errorColumn.ChangeValue(
            str(self.file.errorColumn))
        self.notebook.DVPanel.airmassColumn.ChangeValue(
            str(self.file.airmassColumn))
        self.notebook.GMPanel.SMAPower.ChangeValue(str(self.file.SMAPower))
        self.notebook.GMPanel.SMAEnable.SetValue(self.file.SMAEnable)
        self.notebook.GMPanel.sinusoid.SetValue(self.file.sinusoid)
        self.notebook.GMPanel.normalise.SetValue(self.file.normalise)
        self.notebook.GMPanel.detrend.SetValue(self.file.detrend)
        self.notebook.GMPanel.removeMacro.SetValue(self.file.removeMacro)
        self.notebook.GMPanel.foldingEnable.SetValue(self.file.foldingEnable)
        self.notebook.GMPanel.foldingPeriod.ChangeValue(str(self.file.foldingPeriod))

        self.linesOfData.SetLabel(str(len(self.file.rawData.all)))
        self.startDate.SetLabel(str(al.JDToGC(self.file.rawData.all[0][0]))[:-7])
        self.endDate.SetLabel(
            str(al.JDToGC(self.file.rawData.all[(len(self.file.rawData.all) - 1)][0]))[:-7])
        self.Layout()

    # function to update the header
    def updateHeader(self):
        fileLines = self.file.getLines()
        # Go find the header for this file and print it onto screen
        header = ""

        if not self.file.oneLine:
            for lineIndex, line in enumerate(fileLines):
                if line[0] == "#" or len(line) < 3:
                    header = header + line
                else:
                    header = header + line
                    header = header + fileLines[lineIndex + 1]
                    header = header + fileLines[lineIndex + 2]

                    header = header + "..............\n"

                    # all this while stuff is to avoid annoying commented end areas
                    i = 1
                    lastLine = fileLines[len(fileLines) - i]
                    while lastLine[0] == "#":
                        i += 1
                        lastLine = fileLines[len(fileLines) - i]

                    header = header + fileLines[len(fileLines) - (i + 2)]
                    header = header + fileLines[len(fileLines) - (i + 1)]
                    header = header + lastLine
                    break

        # What to do if its a oneline
        else:
            for line in fileLines:
                if line[0] == "#" or len(line) < 3:
                    header = header + line
                else:
                    numbers = al.findNumbers(line)
                    for numberIndex in range(0, self.file.totalColumns):
                        number = numbers[numberIndex]
                        header = header + "  " + str(number[0])
                    header = header + "\n..............\n"
                    for numberIndex in range(0, self.file.totalColumns):
                        # Wierd necessary inversion
                        numberTemp = (
                            self.file.totalColumns) - numberIndex
                        number = numbers[-numberTemp]
                        header = header + "  " + str(number[0])
        self.mainBox.SetValue(header)

    # What to do on export button press
    def onExport(self, all):
        if not all and self.file is not None:
            export(self.file, self.getUserDirectory())

        elif all and self.file is not None:
            for file in fileList:
                export(file)

    # What to do on save plot button press
    def onSave(self, all):
        if not all and self.file is not None:

            # Open new file
            imageDir = self.file.fileDir[:-self.file.extensionLength] + ".png"
            newFile = Path(imageDir)

            # Make sure it has a unique name
            i = 1
            while True:
                if newFile.is_file():
                    imageDir = imageDir[:-self.file.extensionLength] + "(" + str(i) + ")" + ".png"
                    newFile = Path(imageDir)
                    i += 1
                else:
                    break

            # save it
            self.figure.savefig(imageDir, bbox_inches='tight')

        # elif all and self.file is not None:
        #     for file in fileList:

        #         imageDir = file.fileDir[:-4] + "_raw.astred"
        #         self.figure.savefig(imageDir, bbox_inches='tight')

    # You're leaving!?! What!?
    def OnExit(self, event):
        self.Destroy()  # Close the frame.
        quit()

    # The list of fields and buttons to be anabled at the first file entry
    def enableStuff(self):
        # Enable buttons and things
        if mainFrame.fileNumber.GetValue() == "":
            for directionButton in mainFrame.directionButtons:
                directionButton.Enable()
            for toolButton in mainFrame.toolButtons:
                toolButton.Enable()
            for detector in mainFrame.notebook.DVPanel.detectionFields:
                detector.Enable()
            for detector in mainFrame.notebook.GMPanel.detectionFields:
                detector.Enable()
            for checkbox in mainFrame.notebook.GMPanel.checkboxes:
                checkbox.Enable()
            for button in mainFrame.notebook.GMPanel.buttons:
                if button != mainFrame.notebook.GMPanel.OminusC:
                    button.Enable()

    # Code to combine all data
    def combineCurves(self):
        # Let the user specify the new Dir
        fileDir = self.getUserDirectory()

        if fileDir is None:
            return

        # instantiate our file as a file class
        self.superFile = fileConstruct(fileDir, False)
        self.superFile.foldingEnable = False

        # create a box to store our new data
        superFileData = []

        # Combine all the data
        for file in fileList:
            # normalise the data
            var = al.normalise(file.rawData.JD, file.rawData.var)

            # Collate it all, with HJD
            for HJDIndex, HJD in enumerate(file.rawData.HJD):
                superFileData.append([HJD, var[HJDIndex], file.rawData.error[HJDIndex], file.rawData.airmass[HJDIndex]])

        # sort it
        superFileData.sort(key=al.sort_key)

        # Give our file its data
        for dataPoint in superFileData:
            self.superFile.rawData.HJD.append(dataPoint[0])
            self.superFile.rawData.var.append(dataPoint[1])
            self.superFile.rawData.error.append(dataPoint[2])
            self.superFile.rawData.airmass.append(dataPoint[3])

        # save out the new superfile in astred format
        export(self.superFile, fileDir)

    # A function to allow the user to specify a save file
    def getUserDirectory(self, type=".astred"):

        with wx.FileDialog(self, "Save " + type + " file", wildcard="" + type + " files (*" + type + ")|*" + type + "",
                           style=wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT) as fileDialog:

            if fileDialog.ShowModal() == wx.ID_CANCEL:
                return None    # the user changed their mind

            # save the current contents in the file
            pathName = fileDialog.GetPath()
            return pathName


# a "simple" function for sinusoids used for fitting
def sinusoid(x, a, b, c, d, e, f):
    return (a * np.sin(b * x + c)) + (d * np.sin(e * x + f))


# a not-so simple function for a gaussian curve
def gaussian(x, amplitude, center, width):
    return amplitude * np.exp(-(x - center) ** 2 / width)


# another Gauss function??? in truth i dont know the difference between these two but it just seems better to me. I should have done more maths
def gaussianSigma(x, amplitude, mean, sigma):
    return amplitude * np.exp(-np.power(x - mean, 2.) / (2 * np.power(sigma, 2.)))


# A function to fit a gaussian curve to data
def fitGauss(xBox, yBox, start=None, stop=None, center=None, width=None, output="XandY", fidelity=500):
    # NUMPY ARRAYS FFS
    yBox = np.array(yBox)
    xBox = np.array(xBox)

    # first we need to localize our data to the given window, specified by start and stop. but only if they are provided
    if start is not None and stop is not None:
        indexesToRemove = []
        for xIndex, x in enumerate(xBox):
            if x < start or x > stop:
                indexesToRemove.append(xIndex)

        np.delete(xBox, indexesToRemove)
        np.delete(yBox, indexesToRemove)

        width = (stop - start) / 2  # lets just guess that the user put reasonable start stop values

    else:
        # our width will be based on 10% of the data width (idk can u come up with something better????)
        width = 0.1 * (xBox[-1] - xBox[0])
        start = xBox[0] if start is None else start
        stop = xBox[-1] if stop is None else stop

    # so, what if center isnt given? we just take the x value of the highest y value #HAXXOR
    if center is None:
        center = xBox[np.argmax(yBox)]

    # depending on the type of guassian fuicntion we will use
    # give back data depending on what you want
    if output == "variables":
        initialValues = [1, center, width / 2]  # for [amplitude, mean, sigma]
        try:
            bestValues, covar = optimise.curve_fit(gaussianSigma, xBox, yBox, p0=initialValues)
        except RuntimeError:
            print("No fit found!")
            return None
        return bestValues

    elif output == "XandY":
        initialValues = [1, center, width / 2]  # for [amplitude, mean, sigma]
        try:
            bestValues, covar = optimise.curve_fit(gaussianSigma, xBox, yBox, p0=initialValues, maxfev=10000)
        except RuntimeError:
            print("No fit found!")
            return None, None, None

        # lets also generate a x and y box
        # now we generate the x and y values using our fit
        xBox = np.arange(start, stop, (stop - start) / fidelity)  # new Xbox for known fidelity
        yBox = []
        for x in xBox:
            yBox.append(gaussianSigma(x, bestValues[0], bestValues[1], bestValues[2]))
        return bestValues, xBox, yBox
    else:
        return


# a function to redistribute data and rerun period analysis, to then find the true signal of a dataset. Returns the peak and standard deviation (error)
def periodMC(xBox, yBox, errorBox, minFreq, maxFreq, runs, fidelity=5000, display=True):
    # get some numpy arrays (I admit they are nice)
    xBox = np.array(xBox)
    yBox = np.array(yBox)

    # first, lets get the error. either we use errors provided or generate our own
    if errorBox is None:
        errorBox = np.zeros(yBox.size)
        # We want the errors to be 10% of the average magnitude
        errorBox.fill(0.1 * np.mean(yBox))
    else:
        errorBox = np.array(errorBox)

    # create a box where we will store our peaks
    peakFreqBox = []

    # lets also make a tally for the number of times a peak was not found in our range
    missingPeaks = 0

    # set up gaussX and y stores
    gaussX = np.array([])
    gaussY = np.array([])

    # now we can start iterating through a bunch of simulations
    for run in range(runs):
        # lets create a special box for y values in this run and their offsets
        runY = np.array([])
        offsetMultipliers = np.array([])

        # lets also make an array of multipliers, between -1 and 1, to apply to each error. the first run is displayed so we set the offsets to zero
        if run != 0:
            offsetMultipliers = np.random.uniform(-1, 1, yBox.size)
        else:
            offsetMultipliers = np.zeros(yBox.size)

        # we start by applying an offset to each y value, +- a random number within its corresponding error. ##### THIS COULD PROBABLY ALL BE DONE WITH NUMPY ARRAYS
        # for yIndex, y in enumerate(yBox):
        #     # the offset is the y values corresponding error times its own corresponding multiplier
        #     offset = errorBox[yIndex] * offsetMultipliers[yIndex]
        #     # then we just append the y value plus its offset
        #     np.append(runY, [y + offset])
        # I DID IT ALL WITH NUMPY ARRAYS
        runY = np.add(yBox, np.multiply(errorBox, offsetMultipliers))

        # now we run a lomb scargle on our x values with their new corresponding y values. We speed up the process by limitng our periodogram to the specified window
        frequency, power = LombScargle(xBox, runY).autopower(minimum_frequency=minFreq, maximum_frequency=maxFreq)

        # now we get the frequency of the peak
        peakFreq = frequency[np.argmax(power)]

        # Because Helena told me to, we will fine tune our peak by fitting a gaussian curve to it (now that I think about it this seems like a waste of time to me but oh well)
        # different handling depending on weather this is the first 5 runs or not
        if run < 5:
            gaussVariables, gaussX, gaussY = fitGauss(frequency, power, center=peakFreq, fidelity=(fidelity / 10), output="XandY")
        else:
            gaussVariables = fitGauss(frequency, power, center=peakFreq, fidelity=fidelity, output="variables")

        if gaussVariables is not None:
            # so now we set the "real" peak
            peakFreq = gaussVariables[1]  # this is the mean

            # save the peakfreq to the list of peaks
            peakFreqBox.append(peakFreq)

            # we want to display the first five runs for funsies
            if run < 5:
                # create the plt in the first run
                if run == 0:
                    plt.figure(1)
                    plt.subplot(2, 3, 1)
                    plt.plot(frequency, power, color="b")
                    plt.plot(gaussX, gaussY, color="r")
                else:
                    # add some more plots
                    plt.subplot(2, 3, run + 1)
                    plt.plot(frequency, power, color="b")
                    plt.plot(gaussX, gaussY, color="r")

                # and on the last run we show it!
                if run == 4:
                    plt.show()
        else:
            missingPeaks += 1

    # lets quickly turn peakfreqbox into a numpy array. (it wasnt one earlier bc i dont like numpys append feature)
    peakFreqBox = np.array(peakFreqBox)

    # Phew, that was a lot of hard work. now lets build our "brand new" gaussian
    gaussBins = np.linspace(minFreq, maxFreq, fidelity)
    # OMG numpy.digitize is so cool!
    # indexes = np.digitize(peakFreqBox, gaussBins)  # gives us the "bin" where every frequency should go, eg this goes in bin 3, this bin 5, etc ###### THIS DIDNT WORK; NUMPY:HUSTOGRAM IS BETTER
    # also numpy.bincount is AWESOME!
    # binY = np.bincount(indexes)  # counts the number of times a frequency was assigned to a bin, so for example there were 125 peaks assigned to the center bin, hooray!
    binY, bins = np.histogram(peakFreqBox, gaussBins)  # HA lol this is so much easier. This is clearly the coolest function

    # now for technichal correctness we need our binX to be halfway between each bin, but because they are linearly spaced this is easy. pop the last element too
    binX = gaussBins[:-1] + ((gaussBins[1] - gaussBins[0]) / 2)

    # now we generate our gaussian
    gaussVariables, gaussX, gaussY = fitGauss(binX, binY, output="XandY", fidelity=fidelity)
    if gaussVariables is None:
        print("error fitting to histogram")
        plt.subplot(2, 3, 6)
        plt.hist(peakFreqBox, gaussBins, color="b")
        plt.show()
        return None, None, None

    # now to get the one, the only, peak frequency
    peakFreq = gaussVariables[1]  # this is equal to the "center" or "mean" variable from the gauss function

    # we also want sigma
    sigma = gaussVariables[2]

    # lets shown them the histogram of peaks
    plt.subplot(2, 3, 6)
    plt.hist(peakFreqBox, gaussBins, color="b")
    plt.plot(gaussX, gaussY, color="r")
    # lets limit the width to 10 sigma
    plt.xlim(peakFreq - (5 * sigma), peakFreq + (5 * sigma))
    plt.show()

    # finally, hand those values back
    return peakFreq, sigma, missingPeaks


# A funtion for exporting
def export(file, newFileDir=None):

    checkName = False

    if newFileDir is None:
        # Open new file
        newFileDir = file.fileDir[:-file.extensionLength] + ".astred"
        checkName = True

    newFile = Path(newFileDir)

    if checkName:
        # Make sure it has a unique name if its been generated
        i = 1
        while True:
            if newFile.is_file():
                newFileDir = newFileDir[:-file.extensionLength] + "(" + str(i) + ")" + ".astred"
                newFile = Path(newFileDir)
                i += 1
            else:
                break

    outputFile = open(newFileDir, "w")

    # Write header
    outputFile.write("# Observer: " + str(file.observer) + "\n")
    outputFile.write("# JD Truncation: " + str(file.JDOffset) + "\n")
    outputFile.write("# SMA power: " + str(file.SMAPower) + "\n")
    outputFile.write("# SMA Enable: " + str(file.SMAEnable) + "\n")
    outputFile.write("# detrend: " + str(file.detrend) + "\n")
    outputFile.write("# Remove macro: " + str(file.removeMacro) + "\n")
    outputFile.write("# Oneline: " + str(file.oneLine) + "\n")
    outputFile.write("# isFold: " + str(file.foldingEnable) + "\n")
    outputFile.write("#" + "\n")

    # Get data
    if file.foldingEnable:
        xBox = file.foldX
        yBox = file.foldY
    else:
        xBox = file.rawData.HJD
        yBox = file.rawData.var

    # get error and airmass
    errorBox = file.rawData.error
    airmassBox = file.rawData.airmass

    # Write data
    for xIndex, x in enumerate(xBox):
        # format data
        x = '{0:.20f}'.format(x)
        y = '{0:.20f}'.format(yBox[xIndex])
        error = errorBox[xIndex]
        airmass = airmassBox[xIndex]

        # error and airmass
        error = '{0:.20f}'.format(error)
        if airmass is not None and airmass != "None":
            airmass = '{0:.20f}'.format(airmass)

        # Save it out
        outputFile.write(x + " " + y + " " + error + " " + str(airmass) + "\n")
    outputFile.close()


# Function that adds files to the fileList
def addFiles(fileDirectories):
    changesMade = False
    startingLength = len(fileList) + 1
    firstRun = True if len(fileList) == 0 else False

    # We try because perhaps there is only 1 file being added
    try:
        # Iterate through all the new dirs
        for fileDir in fileDirectories:
            alreadyExisting = False
            # Check to make sure we dont already have this file in our list
            for existingFile in fileList:
                existingFileDir = existingFile.fileDir
                if existingFileDir == fileDir:
                    alreadyExisting = True

            # if they dont exist, add them
            if not alreadyExisting:
                fileList.append(fileConstruct(fileDir))
                changesMade = True

    # Only 1 file
    except TypeError:
        try:
            fileDir = fileDirectories
            for existingFile in fileList:
                existingFileDir = existingFile.fileDir
                if existingFileDir == fileDirectories:
                    alreadyExisting = True

            if not alreadyExisting:
                fileList.append(fileConstruct(fileDir))
                changesMade = True
        except TypeError:
            return

    if changesMade:
        if firstRun:
            mainFrame.enableStuff()
        # Sort the fileList in chronological order
        fileList.sort(key=al.file_sort_key)
        # Place the mainwindow on the latest file
        mainFrame.fileNumber.SetValue(str(startingLength))


# HALP function
def PrintException():
    exc_type, exc_obj, tb = sys.exc_info()
    f = tb.tb_frame
    lineno = tb.tb_lineno
    filename = f.f_code.co_filename
    linecache.checkcache(filename)
    line = linecache.getline(filename, lineno, f.f_globals)
    print('EXCEPTION IN ({}, LINE {} "{}"): {}'.format(
        filename, lineno, line.strip(), exc_obj))
    # A message dialog box with an OK button. wx.OK is a standard ID in wxWidgets.
    dlg = wx.MessageDialog(mainFrame, 'EXCEPTION IN ({}, LINE {} "{}"): {}'.format(filename, lineno, line.strip(), exc_obj), style=wx.OK | wx.CENTRE)
    dlg.ShowModal()  # Show it
    dlg.Destroy()  # finally destroy it when finished.


try:

    # Initiate GUI
    # Create a new app, don't redirect stdout/stderr to a window.
    mainApp = wx.App(False)
    headerFont = wx.Font(20, wx.MODERN, wx.NORMAL, wx.BOLD)
    mainFont = wx.Font(12, wx.MODERN, wx.NORMAL, wx.NORMAL)
    mainFrame = mainWindow(None, "AstroCyles")
    mainApp.MainLoop()

except Exception as err:
    PrintException()
    input()

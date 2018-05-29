# AstroLib
version = 1.0
# Written by Tomass Wilson in python 3.6
# This is a library of functions used by the program AstroRedux

from astropy.time import Time
import numpy as np
from math import floor, pi, sin, cos, tan
from numpy import arctan
import sys

# Stop creating that pesky pycache folder
sys.dont_write_bytecode = True


# JD conversion function
def JDToGC(JD):
    JDTime = Time(JD, format="jd", scale="utc")

    # Convert it to Gregorian
    return JDTime.iso


# Funcion that gets the extension for a file directory
def getExtension(fileDir):
    extension = ""
    extensionLength = None

    for i in range(len(fileDir)):
        # iterates backwards through the directory
        character = fileDir[-i:-(i - 1)]

        # now we know where the extension ends
        if character == ".":
            extensionLength = i
            break

    if extensionLength is not None:
        extension = fileDir[-extensionLength:]
    else:
        return None

    return extension


# a function to normalise a dataset
def normalise(xBox, yBox):
    avg = np.mean(yBox)
    normalisedY = []
    for y in yBox:
        normalisedY.append(y - avg)

    return normalisedY


# Function to remove trends
def removeTrend(x, y):

    # Fit line
    fittedLine = np.polyfit(x, y, 1)
    fittedFunction = np.poly1d(fittedLine)
    lineX = []
    lineY = []
    for xVal in x:
        lineX.append(xVal)
        lineY.append(fittedFunction(xVal))

    # Remove line from data
    yNorm = []
    for yIndex, yVal in enumerate(y):
        yNorm.append(yVal - lineY[yIndex])

    return yNorm


# Function to remove macro
def removeMacro(x, y):

    # Fit macro line
    SMAY = getSMA(x, y, 100, 0.1)

    # Remove line from data
    yNorm = []
    for yIndex, yVal in enumerate(y):
        yNorm.append(yVal - SMAY[yIndex])

    return yNorm


# a function to get a reasonable fidelity, target points is the number of total points any given SMA will consider, so 41 ~> 20 either side
def getFidelity(xBox, targetPoints=41):
    xMin = xBox[0]
    xMax = xBox[-1]
    xRange = xMax - xMin    # gets the range of the dataset
    numOfPoints = len(xBox)

    # Gets the percentage of the dataset that each point should consider
    sectionSize = targetPoints / numOfPoints
    # gets the fidelity
    fidelity = xRange * sectionSize

    # divide by 2 as its per side
    return fidelity / 2


# a function to generate a simple moving average
def getSMA(xBox, yBox, power, fidelity=None, bins=None):
    SMAY = []

    # get fidelity if none supplied
    fidelity = getFidelity(xBox) if fidelity is None else fidelity

    # this is for conditional handling, so the function either uses a supplied set of bins, or the existing set of x values
    if bins is None:
        bins = xBox

    # Iterate through each X value that we want to get an average for
    for xToAv in bins:

        # Find our minimum and maximum x
        minimum = xToAv - fidelity
        maximum = xToAv + fidelity

        # initiate the list for the important y values and their weights
        yBoxWithWeight = []

        # fill that list, by finding y values within the min/max, and calculate their weights
        for xIndex, x in enumerate(xBox):
            if x > maximum:
                break
            elif x > minimum:
                distance = abs(xToAv - x)
                yBoxWithWeight.append([yBox[xIndex], calcWeight(distance, power)])

        # Normalise all the weights so that they sum to 1
        # first get the sum of all weights
        sum = 0
        for box in yBoxWithWeight:
            weight = box[1]
            sum += weight

        # then devide each weight by the sum
        for boxIndex, box in enumerate(yBoxWithWeight):
            weight = box[1]
            yBoxWithWeight[boxIndex][1] = weight / sum

        # then finally get the weighted average
        WA = 0
        for yIndex, yToAv in enumerate(yBoxWithWeight):
            # Value times its weight
            WA += yToAv[0] * yToAv[1]

        # now that we have the wieghted average, append it to the list of Y values
        SMAY.append(WA)

    return SMAY


# a function to calculate weight
def calcWeight(distance, power):
    return ((1 / (distance + 1)) ** power)


# Function to help sort lists of lists by the first nested object
def sort_key(list):
    return list[0]


# A function to help sort files in cronological order
def file_sort_key(file):
    return file.rawData.HJD[0]


# Number checking function
def isNumber(s, type="float"):
    if s == " ":
        return False
    try:
        f = float(s)
        if type == "int":
            return f.is_integer()
        else:
            return True
    except (ValueError, TypeError):
        return False


# Function that finds numbers in text, including negative ones
def findNumbers(line, includeNegative=False, min=None, max=None, offset=None):
    buffer = ""
    numList = []

    # Iterate through each character in a line and add it to the buffer if it is a number, a minus at the start or a period, signifying a decimal
    for characterIndex, character in enumerate(line):
        if isNumber(character) or (includeNegative is True and character == "-" and len(buffer) == 0) or (len(buffer) != 0 and character == "."):
            buffer = buffer + character

        # If the character is not a number, but the buffer has something in it, add that to the numlist and include the index of the first digit in the number
        elif len(buffer) != 0:

            # Make sure its a number first
            if isNumber(buffer):
                numList.append([float(buffer), characterIndex - len(buffer)])
                buffer = ""
            else:
                buffer = ""

    # Catch stragglers, if there isnt whitespace after a line
    if len(buffer) != 0:
        if isNumber(buffer):
            numList.append([float(buffer), characterIndex - len(buffer)])
            buffer = ""
        else:
            buffer = ""

    # Add a third item to each number, signifying "column", or which number in a bare list of numbers is this type of number, as in within the limitations
    for numberIndex, number in enumerate(numList):
        numList[numberIndex].append(numberIndex)

    # Evaluate conditionals
    if (min is not None or max is not None) and len(numList) > 0:
        numbersForRemoval = []

        # First, find numbers to be removed
        for numberBox in numList:
            number = numberBox[0]

            # Add Offset
            if offset is not None:
                number = number + offset

            # Remove number under min
            if min is not None:
                if number < min:
                    numbersForRemoval.append(numberBox)

            # Remove number over max
            if max is not None:
                if number > max:
                    numbersForRemoval.append(numberBox)

        # THEN, remove them. this is done seperately because otherwise the indexes are messed up.
        for numberBox in numbersForRemoval:
            numList.remove(numberBox)

    # Return the list that contains all the numbers in this line
    return numList


# Number offset finder, mostly just a guessing game
def findNumOffset(offsetNum, fileLines, name, min=None, max=None):
    numOffset = None

    # Step 1, check for clues in header
    possibleNumOffsets = []
    for lineIndex, line in enumerate(fileLines):

        if line[0] == "#":
            numbersOnLine = findNumbers(line, False, min, max, offsetNum)
            for number in numbersOnLine:
                possibleNumOffsets.append(number)

    if len(possibleNumOffsets) == 1:
        numOffset = possibleNumOffsets[0][0]
    elif len(possibleNumOffsets) > 1:
        print("ERROR, multiple offsets detected, earliest picked")
        numOffset = possibleNumOffsets[0][0]

    # Step 2, try to add random numbers

    # Check if this is reasonable with user
    if isNumber(numOffset):
        return numOffset
    else:
        print(name, "doesn't look right, but no offset was found")
        return None


# Function that looks for strings of text in a line
def findKeyword(line, keyWord):
    # initiate variables, weather we are iterating through a word (started), and iterator for that word (i), and a list of the "hits" for a keyword, for their indexes
    started = False
    i = 1
    hits = []

    # Go through the characters in  a line
    for characterIndex, character in enumerate(line):
        # Check if we succesfully found a keyword
        if i == len(keyWord):
            i = 1
            started = False
            hits.append(characterIndex - len(keyWord))

        # See if this is the first letter
        if not started and character == keyWord[0]:
            started = True

        # Go through each letter
        elif started and character == keyWord[i]:
            i += 1

        # This particular string isnt a keyword
        elif started:
            started = False
            i = 1

    # Return hits
    if len(hits) == 0:
        return None
    else:
        return hits


# A funtion to find numbers enclosed by prefixes or suffixes or both
def findEnclosedNum(fileLines, prefix=None, min=None, max=None, suffix=None, onlyInHeader=True):
    hits = []

    for line in fileLines:

        if (onlyInHeader and line[0] == "#") or not onlyInHeader:

            # Step 1, find numbers on this line
            numbers = findNumbers(line, True, min, max)
            if len(numbers) != 0:

                # Step 2, find the prefix, remove number if it doesnt fit
                if prefix is not None:
                    numbersForRemoval = []

                    for numberBoxIndex, numberBox in enumerate(numbers):
                        numberIndex = numberBox[1]
                        prefixArea = line[numberIndex -
                                          len(prefix):numberIndex]
                        if prefixArea != prefix:
                            numbersForRemoval.append(numberBox)

                    # remove numbers that dont fit
                    for numberForRemoval in numbersForRemoval:
                        numbers.remove(numberForRemoval)

                # Step 3, find the suffix, remove number if it doesnt fit
                if suffix is not None:
                    numbersForRemoval = []

                    for numberBoxIndex, numberBox in enumerate(numbers):
                        numberIndex = numberBox[1]
                        number = numberBox[0]
                        # this is really complicated so lets explain it
                        # we want to find the suffix, which starts at the end of the number, for that we use the numberindex + the length of the number as a string
                        # then we need the end, which is the same but adding the length of the suffix
                        suffixArea = line[numberIndex + len(
                            str(number)):numberIndex + len(str(number)) + len(suffix)]
                        if suffixArea != suffix:
                            numbersForRemoval.append(numberBox)

                    # remove numbers that dont fit
                    for numberForRemoval in numbersForRemoval:
                        numbers.remove(numberForRemoval)

                # add values to hits, omit their indexes
                for hit in numbers:
                    hits.append(hit[0])

    if len(hits) != 0:
        return hits
    else:
        return None


# gets the contents of a field, sensative to type
def getField(field, type="float"):

    fieldContents = field.GetValue()

    # This should never be used
    if type == "string":
        return fieldContents

    isCompatible = isNumber(fieldContents, type)

    if isCompatible:
        if type == "float":
            return float(fieldContents)
        elif type == "int":
            return int(fieldContents)
    else:
        return None


####################################

# NOT MINE REALLY BAD DONT FIX TOMASS PLZ
# DEFINITION FOR HJD CORRECTION

# Add jd - hjd converter translated by Marcus Levine:

'''
This routine calculates the relative positions of the planets and the Sun at the given jd based on the orbital parameters of the solar system,
transforms the given ra/dec coordinates into a heliocentric reference frame,
and then calculates the difference between the date on the Sun and the date on Earth.
Finally it subtracts this correction value from the original date and returns a full seven digit hjd.
'''


def trans(jd, obj):
    # Compute Positions of All Planets

    rahr = float(obj[0])
    ramin = float(obj[1])
    rasec = float(obj[2])
    decdeg = float(obj[3])
    decmin = float(obj[4])
    decsec = float(obj[5])

    Rads = pi / 180
    nm = []
    el = []
    p = 3

    D = jd - 2451545

    # the Planets

    nm.append("Mercury")
    nm.append("Venus")
    nm.append("Sun")
    nm.append("Mars")
    nm.append("Jupiter")
    nm.append("Saturn")
    nm.append("Uranus")
    nm.append("Neptune")
    nm.append("Pluto")

    # Mercury

    el.append((7.00487 - 0.000000178797 * D) * Rads)
    el.append((48.33167 - 0.0000033942 * D) * Rads)
    el.append((77.45645 + 0.00000436208 * D) * Rads)
    el.append(0.38709893 + 1.80698E-11 * D)
    el.append(0.20563069 + 0.000000000691855 * D)
    el.append((Rads * (252.25084 + 4.092338796 * D)))

    # Venus

    el.append((3.39471 - 0.0000000217507 * D) * Rads)
    el.append((76.68069 - 0.0000075815 * D) * Rads)
    el.append((131.53298 - 0.000000827439 * D) * Rads)
    el.append(0.72333199 + 2.51882E-11 * D)
    el.append(0.00677323 - 0.00000000135195 * D)
    el.append((Rads * (181.97973 + 1.602130474 * D)))

    # Earth

    el.append((0.00005 - 0.000000356985 * D) * Rads)
    el.append((-11.26064 - 0.00013863 * D) * Rads)
    el.append((102.94719 + 0.00000911309 * D) * Rads)
    el.append(1.00000011 - 1.36893E-12 * D)
    el.append(0.01671022 - 0.00000000104148 * D)
    el.append((Rads * (100.46435 + 0.985609101 * D)))

    # Mars

    el.append((1.85061 - 0.000000193703 * D) * Rads)
    el.append((49.57854 - 0.0000077587 * D) * Rads)
    el.append((336.04084 + 0.00001187 * D) * Rads)
    el.append(1.52366231 - 0.000000001977 * D)
    el.append(0.09341233 - 0.00000000325859 * D)
    el.append((Rads * (355.45332 + 0.524033035 * D)))

    # Jupiter

    el.append((1.3053 - 0.0000000315613 * D) * Rads)
    el.append((100.55615 + 0.00000925675 * D) * Rads)
    el.append((14.75385 + 0.00000638779 * D) * Rads)
    el.append(5.20336301 + 0.0000000166289 * D)
    el.append(0.04839266 - 0.00000000352635 * D)
    el.append((Rads * (34.40438 + 0.083086762 * D)))

    # Saturn

    el.append((2.48446 + 0.0000000464674 * D) * Rads)
    el.append((113.71504 - 0.0000121 * D) * Rads)
    el.append((92.43194 - 0.0000148216 * D) * Rads)
    el.append(9.53707032 - 0.0000000825544 * D)
    el.append(0.0541506 - 0.0000000100649 * D)
    el.append((Rads * (49.94432 + 0.033470629 * D)))

    # Uranus

    el.append((0.76986 - 0.0000000158947 * D) * Rads)
    el.append((74.22988 + 0.0000127873 * D) * Rads)
    el.append((170.96424 + 0.0000099822 * D) * Rads)
    el.append(19.19126393 + 0.0000000416222 * D)
    el.append(0.04716771 - 0.00000000524298 * D)
    el.append((Rads * (313.23218 + 0.011731294 * D)))

    # Neptune

    el.append((1.76917 - 0.0000000276827 * D) * Rads)
    el.append((131.72169 - 0.0000011503 * D) * Rads)
    el.append((44.97135 - 0.00000642201 * D) * Rads)
    el.append(30.06896348 - 0.0000000342768 * D)
    el.append(0.00858587 + 0.000000000688296 * D)
    el.append((Rads * (304.88003 + 0.0059810572 * D)))

    # Pluto

    el.append((17.14175 + 0.0000000841889 * D) * Rads)
    el.append((110.30347 - 0.0000002839 * D) * Rads)
    el.append((224.06676 - 0.00000100578 * D) * Rads)
    el.append(39.48168677 - 0.0000000210574 * D)
    el.append(0.24880766 + 0.00000000177002 * D)
    el.append((Rads * (238.92881 + 0.003931834 * D)))

    q = 6 * (p - 1)
    ip = el[q]
    op = el[q + 1]
    pp = el[q + 2]
    ap = el[q + 3]
    ep = el[q + 4]
    lp = el[q + 5]
    ie = el[12]
    oe = el[13]
    pe = el[14]
    ae = el[15]
    ee = el[16]
    le = el[17]

    # Get Earth's position using Kepler's equation

    Me1 = (le - pe)
    B = Me1 / (2 * pi)
    Me1 = 2 * pi * (B - floor(abs(B)))

    if B < 0:
        Me1 = 2 * pi * (B + floor(abs(B)))

    if Me1 < 0:
        Me1 = 2 * pi + Me1

    e = Me1

    delta = 0.05

    while abs(delta) >= pow(10, -12):
        delta = e - ee * sin(e) - Me1
        e = e - delta / (1 - ee * cos(e))

    ve = 2 * arctan(pow(((1 + ee) / (1 - ee)), 0.5) * tan(0.5 * e))

    if ve < 0:
        ve = ve + 2 * pi

    re = ae * (1 - ee * ee) / (1 + ee * cos(ve))

    xe = re * cos(ve + pe)

    ye = re * sin(ve + pe)

    ze = 0

    # Get planet's position using Kepler's equation

    mp = (lp - pp)

    B = mp / (2 * pi)

    mp = 2 * pi * (B - floor(abs(B)))

    if B < 0:
        mp = 2 * pi * (B + floor(abs(B)))

    if mp < 0:
        mp = 2 * pi + mp

    e = mp

    delta = 0.05

    while abs(delta) >= pow(10, -12):
        delta = e - ep * sin(e) - mp
        e = e - delta / (1 - ep * cos(e))

    vp = 2 * arctan(pow(((1 + ep) / (1 - ep)), 0.5) * tan(0.5 * e))

    if vp < 0:
        vp = vp + 2 * pi

    rp = ap * (1 - ep * ep) / (1 + ep * cos(vp))

    xh = rp * (cos(op) * cos(vp + pp - op) - sin(op)
               * sin(vp + pp - op) * cos(ip))

    yh = rp * (sin(op) * cos(vp + pp - op) + cos(op)
               * sin(vp + pp - op) * cos(ip))

    zh = rp * (sin(vp + pp - op) * sin(ip))

    xg = xh - xe

    yg = yh - ye

    zg = zh

    # compute RA and DEC

    ecl = 23.439292 * Rads  # Updated 1-29-2009 429 instead of 439

    xeq = xg

    yeq = yg * cos(ecl) - zg * sin(ecl)

    zeq = yg * sin(ecl) + zg * cos(ecl)

    ra = arctan(yeq / xeq) * 12 / pi

    if xeq < 0:
        ra = ra + 12

    if yeq < 0:
        if xeq > 0:
            ra = ra + 24

    dec = 180 * arctan(zeq / pow((xeq * xeq + yeq * yeq), 0.5)) / pi

    # Sun Coodinates

    xeq = xe

    yeq = ye * cos(ecl) - ze * sin(ecl)

    zeq = ye * sin(ecl) + ze * cos(ecl)

    rae = 12 + arctan(yeq / xeq) * 12 / pi

    if xe < 0:
        rae = rae + 12

    if ye < 0:
        if xe > 0:
            rae = rae + 24

    dece = -180 * arctan(zeq / pow((xeq * xeq + yeq * yeq), 0.5)) / pi

    if p == 3:
        ra = rae
        dec = dece

    if ra < 12:
        raa = 12 - ra

    else:
        raa = 36 - ra

    xpix = 32 + floor(raa * 965 / 24)

    ypix = 201 + floor(dec * -3)

    # compute HJD #

    # Sun

    rah = floor(ra)
    ram = floor((ra - floor(ra)) * 60)
    decd = floor(abs(dec))

    if dec < 0:
        decd = -1 * decd

    decm = floor((abs(dec) - floor(abs(dec))) * 60)

    # Earth

    dec = -1 * dec
    ra = ra + 12

    if (ra > 24):
        ra = ra - 24

    rah = floor(ra)
    ram = floor((ra - floor(ra)) * 60)
    decd = floor(abs(dec))

    if dec < 0:
        decd = -1 * decd

    decm = floor((abs(dec) - floor(abs(dec))) * 60)

    # Object

    ORAH = float(rahr)
    ORAM = float(ramin)
    ORAS = float(rasec)
    ORA = (ORAH + ORAM / 60 + ORAS / 3600) * 15
    ODECD = float(decdeg)
    ODECM = float(decmin)
    ODECS = float(decsec)
    ODEC = abs(ODECD) + ODECM / 60 + ODECS / 3600

    if ODECD < 0:
        ODEC = -1 * ODEC

    # Earth XYZ

    cel = cos(dec * pi / 180)
    earthx = cos(ra * pi / 12) * cel
    earthy = sin(ra * pi / 12) * cel
    earthz = sin(dec * pi / 180)

    # Object XYZ

    cel = cos(ODEC * pi / 180)
    objectx = cos(ORA * pi / 180) * cel
    objecty = sin(ORA * pi / 180) * cel
    objectz = sin(ODEC * pi / 180)

    # Light Time (Minutes per AU)
    ausec = 8.3168775

    correction = ausec * (earthx * objectx + earthy *
                          objecty + earthz * objectz)  # to HJD

    D = D + 2451545  # JD

    D = D + correction / (24 * 60)  # HJD

    return D

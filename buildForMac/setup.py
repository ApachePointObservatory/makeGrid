"""The Python portion of the script that builds TUI

Usage:
% python setup.py [--quiet] py2app
"""
import os
from plistlib import Plist
import shutil
import subprocess
import sys
from setuptools import setup

# add makeGridRoot to sys.path before importing makeGrid
makeGridRoot = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path = [makeGridRoot] + sys.path
from makeGrid import __version__

appName = "MakeGrid"
mainProg = os.path.join(makeGridRoot, "runMakeGrid.py")
# iconFile = "%s.icns" % appName
appPath = os.path.join("dist", "%s.app" % (appName,))
contentsDir = os.path.join(appPath, "Contents")
fullVersStr = __version__
shortVersStr = fullVersStr

inclModules = (
#    "FileDialog",
)
# packages to include recursively
inclPackages = (
#    "RO",
)

plist = Plist(
    CFBundleName                = appName,
    CFBundleShortVersionString  = shortVersStr,
    CFBundleGetInfoString       = "%s %s" % (appName, fullVersStr),
    CFBundleExecutable          = appName,
    LSMinimumSystemVersion      = "10.6.0",
    LSArchitecturePriority      = ("i386",) # force 32-bit mode;
        # this is needed for Tcl/TK 8.5.11 to run on MacOS X 10.9;
        # I'm stuck with 8.5.11 due to a crashing bug in Tcl/Tk 8.5.12 - 8.5.15.1
)

setup(
    app = [mainProg],
    setup_requires = ["py2app"],
    options = dict(
        py2app = dict (
            plist = plist,
            # iconfile = iconFile,
            # includes = inclModules,
            # packages = inclPackages,
        )
    ),
)

# Delete Tcl/Tk documentation
tclFrameworkDir = os.path.join(contentsDir, "Frameworks", "Tcl.framework")
tclDocDir = os.path.join(tclFrameworkDir, "Resources", "English.lproj", "ActiveTcl-8.4")
if os.path.isdir(tclFrameworkDir):
    print "*** Tcl/Tk Framework is part of the application package ***"
    if os.path.isdir(tclDocDir):
        # Delete extraneous files
        print "*** Removing Tcl/Tk help from the application package ***"
        shutil.rmtree(tclDocDir)
else:
    print "*** WARNING: Tcl/Tk Framework is NOT part of the application package ***"

print "*** Creating disk image ***"
appName = "%s_%s" % (appName, shortVersStr)
destFile = os.path.join("dist", appName)
args=("hdiutil", "create", "-srcdir", appPath, destFile)
retCode = subprocess.call(args=args)

print "*** Done building %s ***" % (appName,)

#!/usr/bin/env python2
from __future__ import division, absolute_import

import Tkinter

# the following is needed by makeGridWdg to avoid warnings
import RO.Comm.Generic
RO.Comm.Generic.setFramework("tk")

from makeGrid import MakeGridWdg

root = Tkinter.Tk()
wdg = MakeGridWdg(root)
wdg.pack(expand=True, fill="both")
root.mainloop()

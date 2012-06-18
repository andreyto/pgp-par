### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


from MGT.Options import *
import os
pjoin = os.path.join
pabs = os.path.abspath


class MGTOptions(Options):
    def __init__(self,*l,**kw):
        Options.__init__(self,*l,**kw)
        self.debug = 1

options = MGTOptions()


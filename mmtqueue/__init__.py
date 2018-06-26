# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
astroplan is an open source (BSD licensed) observation planning package for
astronomers that can help you plan for everything but the clouds.

It is an in-development `Astropy <http://www.astropy.org>`__
`affiliated package <http://www.astropy.org/affiliated/index.html>`__ that
seeks to make your life as an observational astronomer a little less
infuriating.

* Code: https://github.com/astropy/astroplan
* Docs: https://astroplan.readthedocs.io/
"""
from __future__ import (absolute_import, division, print_function,
                                    unicode_literals)

# Affiliated packages may add whatever they like to this file, but
# should keep this content at the top.
# ----------------------------------------------------------------------------
from ._astropy_init import *
# ----------------------------------------------------------------------------

# For egg_info test builds to pass, put package imports here.
# if not _ASTROPY_SETUP_:
#     get_IERS_A_or_workaround()



# Application specific modules
from .utilities import *

from .utils import *
from .observer import *
from .target import *
from .exceptions import *
from .moon import *
from .constraints import *
from .scheduling import *
from .periodic import *

# Combining all of the names
__all__ = (constraints.__all__
        + scheduling.__all__
        + utilities.__all__
        + utils.__all__
        + observer.__all__
        + target.__all__
        + exceptions.__all__
        + moon.__all__
        + constraints.__all__
        + scheduling.__all__
        + periodic.__all__)

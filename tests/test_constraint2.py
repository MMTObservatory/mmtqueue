#!/usr/bin/env python3

from astroplan import Constraint
from astroplan.constraints import _get_altaz

class AzimuthConstraint(Constraint):
    """Constrain the azimuth of the target.

    longish explanation: returns the square of a: :math:`a^2`
    
    .. note::
        This can misbehave if you try to constrain negative altitudes, as
        the `~astropy.coordinates.AltAz` frame tends to mishandle negative

    :param min: `~astropy.units.Quantity` or `None`
        Minimum altitude of the target (inclusive). `None` indicates no limit.
    :param max: `~astropy.units.Quantity` or `None`
        Maximum altitude of the target (inclusive). `None` indicates no limit.
    :param boolean_constraint: bool
        If True, the constraint is treated as a boolean (True for within the
        limits and False for outside).  If False, the constraint returns a
        float on [0, 1], where 0 is the min altitude and 1 is the max.
    
    :returns: a*a
    
    """
    def __init__(self, min=None, max=None, boolean_constraint=True):
        if min is None:
            self.min = -90*u.deg
        else:
            self.min = min
        if max is None:
            self.max = 90*u.deg
        else:
            self.max = max

        self.boolean_constraint = boolean_constraint

        # Change default behavior of gridding of results to True
        self.grid_times_targets=True


    def compute_constraint(self, times, observer, targets):
        """
        :param min: `~astropy.units.Quantity` or `None`
            Minimum altitude of the target (inclusive). `None` indicates no limit.
        :param max: `~astropy.units.Quantity` or `None`
            Maximum altitude of the target (inclusive). `None` indicates no limit.
        :param boolean_constraint: bool
            If True, the constraint is treated as a boolean (True for within the
            limits and False for outside).  If False, the constraint returns a
            float on [0, 1], where 0 is the min altitude and 1 is the max.
    
        :returns: a*a
        """
        
        cached_altaz = _get_altaz(times, observer, targets)
        az = cached_altaz['altaz'].az
        if self.boolean_constraint:
            lowermask = self.min <= az
            uppermask = az <= self.max
            return lowermask & uppermask
        else:
            return _rescale_minmax(az, self.min, self.max)


# From:
# http://astroplan.readthedocs.io/en/stable/tutorials/constraints.html

from astroplan import Observer, FixedTarget

from astropy.time import Time


subaru = Observer.at_site("Subaru")
time_range = Time(["2017-08-01 06:00", "2017-08-01 12:00"])

# Read in the table of targets
from astropy.table import Table
target_table = Table.read('targets.txt', format='ascii')

# Create astroplan.FixedTarget objects for each one in the table
from astropy.coordinates import SkyCoord
import astropy.units as u

from mmtqueue.targets import MMTFixedTarget
targets = [FixedTarget(coord=SkyCoord(ra=ra*u.deg, dec=dec*u.deg), name=name)
           for name, ra, dec in target_table]
#targets = [MMTFixedTarget(coord=SkyCoord(ra=ra*u.deg, dec=dec*u.deg), name=name)
#           for name, ra, dec in target_table]

from astroplan import (AltitudeConstraint, 
		       AirmassConstraint,
                       AtNightConstraint)

for target in targets:
    print("Target: ", target.name)
    if False:
        print("Posang: ", target.posang)


# constraints = [AltitudeConstraint(10*u.deg, 80*u.deg),
#	       AzimuthConstraint(0*u.deg, 180*u.deg),
#               AirmassConstraint(5), AtNightConstraint.twilight_civil()]
constraints = [AltitudeConstraint(10*u.deg, 80*u.deg),
               AirmassConstraint(5), 
               AzimuthConstraint(0*u.deg, 10*u.deg),	       
               AtNightConstraint.twilight_civil()]

from astroplan import is_observable, is_always_observable, months_observable
# Are targets *ever* observable in the time range?
ever_observable = is_observable(constraints, subaru, targets, time_range=time_range)

# Are targets *always* observable in the time range?
always_observable = is_always_observable(constraints, subaru, targets, time_range=time_range)

# During what months are the targets ever observable?
# best_months = months_observable(constraints, subaru, targets)

import numpy as np
observability_table = Table()
observability_table['targets'] = [target.name for target in targets]
observability_table['ever_observable'] = ever_observable
observability_table['always_observable'] = always_observable
print(observability_table)





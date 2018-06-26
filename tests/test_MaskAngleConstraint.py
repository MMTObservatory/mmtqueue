#!/usr/bin/env python3

from astroplan import Constraint
from astroplan.constraints import _get_altaz

from astroplan import (AltitudeConstraint, 
                AirmassConstraint,
                AtNightConstraint)
	
from mmtqueue.constraints import MaskAngleConstraint

# From:
# http://astroplan.readthedocs.io/en/stable/tutorials/constraints.html

from astroplan import Observer, FixedTarget

from astropy.time import Time


subaru = Observer.at_site("Subaru")
time_range = Time(["2017-08-01 06:00", "2017-08-01 12:00"])

# Read in the table of targets
# Read in the table of targets
from astropy.table import Table
# target_table = Table.read('targets.txt', format='ascii')
target_table = Table.read('targets_posang.txt', format='ascii')

# Create astroplan.FixedTarget objects for each one in the table

# Including the SkyCoordInfo mixin
from astropy.coordinates import SkyCoord, SkyCoordInfo
import astropy.units as u

from mmtqueue.targets import MMTFixedTarget

print("Hi1")

info=[SkyCoordInfo(unit=parang_design) for name, ra, dec, parang_design in target_table]
targets = [MMTFixedTarget(coord=SkyCoord(ra=ra*u.deg, dec=dec*u.deg), info=info, name=name) for name, ra, dec, parang_design in target_table]

print("Hi2")


for target in targets:
    print("Target: ", target.name)
    if False:
        print("Posang: ", target.posang)


print("Hi3")


# constraints = [AltitudeConstraint(10*u.deg, 80*u.deg),
#	       AzimuthConstraint(0*u.deg, 180*u.deg),
#               AirmassConstraint(5), AtNightConstraint.twilight_civil()]
constraints = [AltitudeConstraint(10*u.deg, 80*u.deg),
                MaskAngleConstraint(parang_limit=30*u.deg),
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





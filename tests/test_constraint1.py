#!/usr/bin/env python3

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

targets = [FixedTarget(coord=SkyCoord(ra=ra*u.deg, dec=dec*u.deg), name=name)
           for name, ra, dec in target_table]

from astroplan import (AltitudeConstraint,
                       AirmassConstraint,
                       AtNightConstraint)

for target in targets:
  print("Target: ", target.name)
  #if False:
      #print("Posang: ", target.posang)


constraints = [AltitudeConstraint(10*u.deg, 80*u.deg),
               AirmassConstraint(5), 
               AtNightConstraint.twilight_civil()]

#from astroplan import is_observable, is_always_observable, months_observable
## Are targets *ever* observable in the time range?
#ever_observable = is_observable(constraints, subaru, targets, time_range=time_range)

## Are targets *always* observable in the time range?
#always_observable = is_always_observable(constraints, subaru, targets, time_range=time_range)

## During what months are the targets ever observable?
## best_months = months_observable(constraints, subaru, targets)

#import numpy as np
#observability_table = Table()
#observability_table['targets'] = [target.name for target in targets]
#observability_table['ever_observable'] = ever_observable
#observability_table['always_observable'] = always_observable
#print(observability_table)





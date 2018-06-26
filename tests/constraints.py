"""
import pkg_resources
pkg_resources.require("astropy>=2.0")
pkg_resources.require("astroplan>=0.4"

from astroplan import Observer, FixedTarget
from astropy.time import Time, TimeDelta
from astroplan import Constraint, AtNightConstraint, AirmassConstraint
from astroplan import SequentialScheduler, ObservingBlock
from astroplan.constraints import _get_altaz, is_event_observable
from astropy.coordinates import Angle
from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np
import datetime
import requests
import json
import sqlite3
from sqlite3 import Error
"""
from astroplan import Constraint
import astropy.units as u
import sqlite3
from sqlite3 import Error


class Logger():
    
    def __init__(self, 
        database, 
        table,
        debug=False):

        self.database = database
        self.table = table
        self.debug = debug

        self.conn = sqlite3.connect(self.database)
        self.conn.row_factory = self.dict_factory


    def get_score(key):
        """

        """
        c = self.conn.cursor()
        c.execute("SELECT * FROM scores WHERE key='{}'".format(key)) 
        row = c.fetchone()
        return row


    def set_score(key, value):
        """

        """
        global conn
        c = conn.cursor()
        
        try:
            sql = "REPLACE INTO scores VALUES ('{}', {})".format(key, value)
            if True:
                print("sql: ", sql)
            
            # Can do this with new versions of sqlite.  It does commits automatically.
            with conn:
                conn.execute(sql)

            # Older version.
            # c.execute(sql)
            # Save (commit) the changes
            # conn.commit()
            # conn.close()
        except sqlite3.IntegrityError:
            print('ERROR: ID already exists in PRIMARY KEY column {}'.format("key"))

    def get_constraint_name(constraint):
        return  type(constraint).__name__

    def get_key(constraint, target, time):
        """
        
        """
        time.format = 'isot'
        constraint_name = get_constraint_name(constraint)
        key = "{}.{}.{}".format(target, 
                                str(time), 
                                constraint_name)
        if True:
            print("key:", key)
        return key

    def dict_factory(cursor, row):
        d = {}
        for idx,col in enumerate(cursor.description):
            d[col[0]] = row[idx]
        return d


# Define an astroplan constraint for the distance of the target from the meridian.
# The returned value is either a boolean [0,1] if the target is outside of an allowed time
# from meridian transit or a float from [0.0:1.0], where the value is 1.0
# when the target is o/typen the meridian to 0.0 when it is at the anti-meridian (12 hours from the meridian)
class MeridianConstraint(Constraint):
    """Constrains the time for targets from meridian transit.

    Principal investigators (PI's) are required to assigned an integer priority from 1 (highest) to 3 (lowest) to each of their targets.
    The targets should be equally divided into the three priorities (i.e., 1, 2, and 3) so that 1/3 of the requested time correspondes to each
    priority.
    This equal division into the three priorities by even time requested is needed to keep scheduling fair for all projects.
    Every effort will be made to observe all targets, but PI's should anticipate that at least part of their priority 3 targets will not be observed because of poor weather or other causes.


    """
    def __init__(self, mode="sunset",
                 min_alt_degrees=20 * u.deg,
                 max_solar_altitude=-12 * u.deg,
                 grid_times_targets=False, 
                 debug=False):
        """
        Parameters
        ----------
        max : `~astropy.units.Quantity` or `None` (optional)
            Maximum acceptable separation (in decimal hours) between meridian and target (inclusive).
            `None` indicates no constraint of how far the target can be from the meridian.
        boolean_constraint : bool
            If True, the constraint is treated as a boolean (True for within the
            limits and False for outside).  If False, the constraint returns a
            float on [0, 1], where 0 is when the target is on the anti-meridian and
            1 is when the target is on the meridian.
        """
        self.mode = mode 
        self.min_alt_degrees = min_alt_degrees
        self.max_solar_altitude = max_solar_altitude
        self.grid_times_targets = grid_times_targets
        self.debug = debug
        

    def compute_constraint(self, times, observer, targets):
        """
        The MeridianConstraint is calculated by: 1) determining the number of hours the target is from the meridian, and 2) calculating a constraint using Math.abs((12. - hours_from_meridan)/12.0) for the target's position at the beginning, middle, and end of the observing block.  
        The calculated scores for these three times will be different.  
        This causes the constraint to equal 1.0 on the meridian and 0.0 on the anti-meridian (12 hours away).  
        Since the absolute value is used, it doesn't matter which direction the target is from the meridan, i.e., positive hours or negative hours.  Values will always vary from 1.0 to 0.0
        It is possible that the target passes through the meridian during the observing block, i.e., it "transits". 
        Caution should be used in cases where the target transits in that azimuth velocities can be very large if the target is close to zenith.
        The maximum AltitudeConstraint should help prevent extremely large azimuth velocities.  

        It should be remembered that constraint scores are calculated at the beginning, middle, and end of each observing block as part of score for the block.
        This causes the constraint to be multiplied by itself three times and the constraint to vary as 1/X^^3 rather than 1/X.

        """
        
        # 12 hours: the maximum possible time for a target to be from the meridian 
        seconds_in_12hrs = 43200     # 12 hours ==> 12 * 60 * 60 = 43200 seconds
        
        # Set up to TimeDelta constatnt values for future use.
        dt_0hrs = TimeDelta(0,format='sec')
        dt_1hrs = TimeDelta(3600,format='sec')
        dt_1_5hrs = TimeDelta(3600*1.5,format='sec')
        dt_2hrs = TimeDelta(3600*2.0,format='sec')
        dt_2_5hrs = TimeDelta(3600*2.5,format='sec')
        dt_3hrs = TimeDelta(3600*3.0,format='sec')
        dt_3_5hrs = TimeDelta(3600*3.5,format='sec')
        dt_4hrs = TimeDelta(3600*4.0,format='sec')
                
        
        
        # This list will eventually be of length = len(targets) * len(times).
        # We'll reshape it when done.
        # To do:  See if there is a way to user more numpy functions for these
        # calculations to speed things up.
        mask=[]
        
        for target in targets:
            print("target: ", target)
            """
            if False:
                # This is a short-circuit to speed up computation.
                #
                # If the altitude is below the minimum allowed altitude,
                # assign zero's to all times for this target.
                # 
                # In reality, the scores for one or more times within the 
                # observing block could be greater than zero, but
                # the constraint will fail overall for the target since 
                # there is at least one zero score.
            
                cached_altaz = _get_altaz(times, observer, target)
                alt = cached_altaz['altaz'].alt
                
                # Step 1: Check if altitudes are below the minimum allowed
                #       altitude.  This allows short-circuiting of the constraint
                #       calculation to improve efficiency.
                #       If the target is below the minimum altitude,
                #       the overall score for the target will be zero anyway.
                alt_check = alt < self.min_alt_degrees

                # Step 2: This is the special case where the target is below 
                # the lower allowed altitude at any time of the time slot.  
                # The scores for all times for this target will be set to zero.
                if alt_check.any() == True:
                    if self.debug:
                        print("alt is less than minimum: %s, %s" % (alt.degree, self.min_alt_degrees))
                    
                    # Adding a list of zeros of length = len(times)
                    mask += [0.0] * len(times)
            
            # Otherwise, we need to look in detail for each time.
            else:
            """
            if True:
                for time in times:
                    key = get_key(self, target.name, time)
                    obj = get_score(key)
                    if obj is not None:
                        score = obj['value']
                        if True:
                            print("Using saved score.", score)  
                    
                    else:
                        # We do a series of tests to see if the observing block is 
                        # is setting early in the evening.  We want to give high
                        # priority to PIPriority == 1 observing blocks that are
                        # setting near sunset and that can still be observed.
                        #
                        # We take into account the duration of observing blocks
                        # when doing this special "sunset" mode.
                        #
                        # Step 3:  Get the time of the previous sunset.
                        #       This will be used to see if the time is close to sunset.
                        prev_sun_set_time = observer.sun_set_time(time, 
                                            which='previous', 
                                            horizon=self.max_solar_altitude)

                        # Step 4: Calculate the time from the previous sunset.  
                        #        This is an indication of how close we are to sunset.
                        #        It will be a small number if we are trying to observe
                        #        just after sunset.
                        td1 = time - prev_sun_set_time 

                        # Step 5: Get the time of the next target rise.
                        #       This will be used to see if the target is close to rising
                        #       in the east in the morning.
                        next_target_set = observer.target_set_time(time,target,
                                            which="next", 
                                            horizon=self.min_alt_degrees )

                        # Step 6: Calculate the time from the next target setting.
                        #        This number will be a small positive number when 
                        #        the target is above the western horizon.
                        td2 = next_target_set - time

                        # if False and verbose:
                        #    print("tx1: {}, sec: {}, tx2: {}, sec: {}".format(tx1, tx1.sec, tx2, tx2.sec))

                        # Note: The next three steps use a "sunset" mode where we want to give
                        #       high priority to targets that will be setting within the next 2-4
                        #       hours.  We use a graded approach for scoring, based on the duration
                        #       of the target/observing block.

                        # Step 8: This is the "2-hour-target-duration" sunset special case.  
                        #       Evaluate for targets/observing blocks that are more the
                        #       _2_ hours in duration, and we are within _4_ hours after sunset,
                        #       and the target will be setting within _4_ hours.
                        #       Only do this for priority 1 targets.
                        #       This is our only chance to observe them.
                        #       The score is set to 1.0 to give it the maximum chance of being observed.

                        if self.mode == "sunset" and \
                                td1 <= dt_4hrs and td1 > dt_0hrs and \
                                td2 <= dt_4hrs and td2 > dt_0hrs and \
                                target.duration >= 2.0 * u.hour and \
                                target.pi_priority == 1.0:
                            score = 1.0
                            if verbose:
                                print("Sunset special case (>= 2-hr duration), score:",score)

                        # Step 9: This is the "1-hour-target-duration" sunset special case.  
                        #       Evaluate for targets/observing blocks that are 1-2 hours
                        #       in duration, and we are within _3_ hours after sunset,
                        #       and the target will be setting within _3_ hours.
                        #       Only do this for priority 1 targets.
                        #       This is our only chance to observe them.
                        #       The score is set to 1.0 to give it the maximum chance of being observed.
                        elif self.mode == "sunset" and \
                                td1 <= dt_3hrs and td1 > dt_0hrs and \
                                td2 <= dt_3hrs and td2 > dt_0hrs and \
                                target.duration >= 1.0 * u.hour and \
                                target.pi_priority == 1.0:
                            score = 1.0
                            if verbose:
                                print("Sunset special case, (>= 1-hour and < 2-hour duration) score:",score)

                        # Step 10: This is the "<1-hour-target-duration" sunset special case.  
                        #       Evaluate for targets/observing blocks that are <1 hour
                        #       in duration, and we are within _2_ hours after sunset,
                        #       and the target will be setting within _2_ hours.
                        #       Only do this for priority 1 targets.
                        #       This is our only chance to observe them.
                        #       The score is set to 1.0 to give it the maximum chance of being observed.
                        elif self.mode == "sunset" and \
                                td1 <= dt_2hrs and td1 > dt_0hrs and \
                                td2 <= dt_2hrs and td2 > dt_0hrs and \
                                target.pi_priority == 1.0:
                            score = 1.0
                            if verbose:
                                print("Sunset special case, (any duration) score:",score)

                        # Step 11:  If all of the other conditions have not been true,
                        #       Determine how far the target is from the meridian in seconds
                        #       and divide by 12 hours (== 43200 seconds) 
                        #       The target can be in either rising towards the meridian or
                        #       setting away from the meridian.
                        else:

                            meridian_time = observer.target_meridian_transit_time(time,target,which='nearest')
                            diff = abs(time.unix - meridian_time.unix)

                            # There are times when the meridian time is 24 hours off.
                            # So, the math here accounts for that.
                            # If the time difference is more than 24 hours (43200 seconds), 
                            # subtract 24 hours.
                            if diff > 43200:
                                diff -= 86400 # 24 hours in seconds.
                                # Recheck that we are using an absolute value.
                                diff = abs(diff)

                            # Here is the meridian scoring algorithm.  
                            # The closer to the meridian the closer the score is to one.  
                            # The range of scores is 1.0 (on the meridian) to 
                            # 0.0 (on the anti-meridian).
                            score = 1.0 - (diff / 43200.0)               
                            # Shouldn't need these, but just checking.
                            # The score should already range from 0.0 to 1.0.
                            if score < 0.0:
                                score = 0.0
                            if score > 1.0:
                                score = 1.0
                    
                        # Add to sqlite database
                        set_score(key,score)
                    # Add the new score to the 1-D list.
                    mask.append(score)
            
        if self.debug:
            print("mask")
            print(repr(mask))
            
        # Turn the mask into a numpy array and reshape.
        mask = np.reshape(np.array(mask),[len(targets), len(times)])
            
        if self.debug:
            print("mask")
            print(repr(mask))
        

#!/usr/bin/python3.6
"""
mmtscheduler.py

Calculates a queue schedule for MMIRS or Binospec at the MMTO.

# ##  
# Note: uses astropy >= 2.0.
# astroplan 0.4 files are included here under mmtqueue.
# External dependencies to astroplan are removed.
# (astropy dependency is still there...)  
# Numpy is used when possible.

# Removing multiprocessing on 2018-04-04.  Using the global daily_events
# for optimization instead.
#
# Going to astropy 3.0, 2018-02-19
# It should be backward compatible with astropy 2.0.
#
# Going to a multiprocessing version, 2018-02-19
# The best performance seems to be with 3 worker processes that calculate
# the constraint scores for the start, middle and end time for an observing
# block.
# This speeds up calculations by about 30%, i.e., 90 seconds goes to 60 seconds.
#

    parser.add_argument('-m','--mode', help='Scheduling mode.  mode 1: "scheduler" mode; mode 2: "dispatcher" mode.', type=int, choices=[1,2], required=False)
    parser.add_argument('-d','--debug', help='Print debugging messages', required=False)
    parser.add_argument('-v','--version', help='Software version.  Release 1.0 is on ops; version 2.0 and above are on scheduler', type=float, required=False)
    parser.add_argument('-s','--start_date', help='Start time for scheduling (UTC time), e.g., "2017-10-01T03:30:00".  Defines the time used for dispatcher mode and the start time for scheduler mode', required=False)
    parser.add_argument('-e','--end_date', help='End time for scheduling (UTC time), e.g., e.g., "2017-10-01T12:30:00".  Defines the end time for scheduler mode.', required=False)
    parser.add_argument('-i','--schedule_id', help='Queue scheduling run id.', type=int, required=False)
    # New 2017-10-16.  Adding binospec option for instrument.  Defaults to type string.
    parser.add_argument('-t','--instrument', help='Instrument  "mmirs" or "binospec"', required=False)

    Command line arguments
    ----------
    -m, --mode : integer
        Mode for computation.  The scheduler mode computes a schedule for the entire queue; the dispatcher mode
        schedules the first observing block
        Example:
            1 (scheduler), e.g., -m=1,
            2 (dispatcher), e.g., -m=2

    -i, --schedule_id : integer
        Schedule ID.  This is the database ID corresponding to the schedule.  A new ID is assigned each time
        a new schedule is generated.  You can get the schedule_id from the ObservatoryManager web interface.
        Example:
            -i=836

    -d, --debug : string
        A debugging flag: 'true' for verbose output.  Defaults to 'false' so only needed if you want
        debugging output.
        Example:
            -d=true
    -s :


Parameter:

Example:
    ./mmtscheduler.py -m=2 -i=836


"""

import argparse
import datetime
import json as json
import numpy as np
from pytz import timezone
import redis as redis
import sqlite3 as sqlite3
import time as time
import subprocess

# Trying this fix for the malformed finals2000A.all at the US Naval Observatory
from astropy.utils import iers
iers.IERS_A_URL = 'http://toshi.nofs.navy.mil/ser7/finals2000A.all'
# Using this to force a new auto_download.  
# iers.auto_max_age = 0.001
iers.conf.auto_download = True

import traceback as traceback

# import warnings

# These are the Classes from astroplan that have been incorporated into the 'mmtqueue' module.
# from astroplan import Observer
# from astroplan import FixedTarget
# from astroplan import Constraint
# from astroplan import AirmassConstraint
# from astroplan import AltitudeConstraint
# from astroplan import AtNightConstraint
# from astroplan import TimeConstraint
# from astroplan import MoonSeparationConstraint
# from astroplan import ObservingBlock
# from astroplan import Transitioner
# from astroplan import Slot

# Only using if we program in short-circuiting for constraint calculations.
# from astroplan.constraints import _get_altaz
# from astropy.coordinates import Angle
from astropy.coordinates import SkyCoord
from astropy.coordinates import AltAz
from astropy.coordinates import EarthLocation
from astropy.time import Time
from astropy.time import TimeDelta
import astropy.units as u

from mmtqueue.constraints import AirmassConstraint
from mmtqueue.constraints import AltitudeConstraint
from mmtqueue.constraints import AtNightConstraint
from mmtqueue.constraints import MoonSeparationConstraint

# New (or modified) constraints.
from mmtqueue.constraints import MaskAngleConstraint
from mmtqueue.constraints import MaskNumberConstraint
from mmtqueue.constraints import MeridianConstraint
from mmtqueue.constraints import PIPriorityConstraint
from mmtqueue.constraints import RotatorConstraint
from mmtqueue.constraints import TimeAllocationConstraint
from mmtqueue.constraints import TimeConstraint

from mmtqueue.constraints import MaskReadyConstraint
from mmtqueue.constraints import TransitConstraint
from mmtqueue.constraints import ObjectTypeConstraint
from mmtqueue.constraints import TargetOfOpportunityConstraint

from mmtqueue.constraints import is_always_observable

from mmtqueue.observer import Observer

from mmtqueue.scheduling import Schedule
from mmtqueue.scheduling import Scheduler
from mmtqueue.scheduling import PriorityScheduler
from mmtqueue.scheduling import Scorer
# from mmtqueue.scheduling import SequentialScheduler
from mmtqueue.scheduling import ObservingBlock
from mmtqueue.scheduling import Transitioner
from mmtqueue.scheduling import Slot

from mmtqueue.target import FixedTarget

# These are my specific utility functions.
from mmtqueue.utilities import add_prefix
from mmtqueue.utilities import dict_factory
from mmtqueue.utilities import get_config_json
from mmtqueue.utilities import get_constraint_name
from mmtqueue.utilities import get_key
from mmtqueue.utilities import get_short_key
from mmtqueue.utilities import get_name
from mmtqueue.utilities import get_mask_id
from mmtqueue.utilities import get_program_id
from mmtqueue.utilities import get_score
from mmtqueue.utilities import set_score
from mmtqueue.utilities import set_scorer
from mmtqueue.utilities import get_now
from mmtqueue.utilities import get_target_name
from mmtqueue.utilities import get_target_key
from mmtqueue.utilities import get_sun_moon_date_key
from mmtqueue.utilities import get_target_date_key
from mmtqueue.utilities import is_numeric
from mmtqueue.utilities import localize_time
from mmtqueue.utilities import to_allocated_time
from mmtqueue.utilities import to_constraint_details
from mmtqueue.utilities import to_completed_time
from mmtqueue.utilities import to_mysql
from mmtqueue.utilities import to_output
from mmtqueue.utilities import to_queue
from mmtqueue.utilities import to_redis
from mmtqueue.utilities import to_status
from mmtqueue.utilities import to_stdout
from mmtqueue.utilities import roundup
from mmtqueue.utilities import roundHour
from mmtqueue.utilities import roundTime
from mmtqueue.utilities import table_to_entries

from mmtqueue.utils import time_grid_from_range

# Global variables.
# (Yes, it's best to avoid global variables, but these variables are used by multiple modules.)

# Global verbosity flags.
verbose = True
verbose2 = False
verbose3 = False
verbose4 = False

doing_line_profiling = False

# Global constraints flags:  whether a global constraint is to be used in calculations.
# Each variable can be set to either True or False
#
# These parameters are used to turn each constraint on or off in this code.
# There are similar parameters loaded in the configuration structure during setup
# from the user.
#
# Both must be True for the constraint to be used.
use_AirmassConstraint = True
use_AltitudeConstraint = True
use_AtNightConstraint = True
use_MaskAngleConstraint = True
use_MaskNumberConstraint = False
use_MaskReadyConstraint = False
use_MeridianConstraint = True
use_MoonSeparationConstraint = True
use_ObjectTypeConstraint = True
use_PIPriorityConstraint = True
use_RotatorConstraint = True
use_TimeAllocationConstraint = True
use_TimeConstraint = True
use_TargetOfOpportunityConstraint = True
use_TransitConstraint = True


# Global constraint debugging flags:  these variables provide detailed outputs for the constraint
# Each variable can be set to either True or False
test_AirmassConstraint = False
test_AltitudeConstraint = False
test_AtNightConstraint = False
test_MoonSeparationConstraint = False
test_MaskAngleConstraint = False
test_MaskNumberConstraint = False
test_MaskReadyConstraint = False
test_MeridianConstraint = False
test_ObjectTypeConstraint = False
test_PIPriorityConstraint = False
test_RotatorConstraint = False
test_TimeAllocationConstraint = False
test_TimeConstraint = False
test_TargetOfOpportunityConstraint = False

# Global sequential scheduling debugging flag:  used for detailed testing of the MMTSequentialScheduler.
test_MMTSequentialScheduler = False

# Global flag on whether we want to post values to the Redis server.
do_redis = True
# The master Redis server.
redis_host = 'redis.mmto.arizona.edu'

# Lambda functions for reuse.
current_milliseconds = lambda: int(round(time.time() * 1000))
current_datetime = lambda: str(timezone('America/Phoenix').localize(datetime.datetime.now()))


# Global schedule id variable:  Typically, reset during configuration
# Notes:
# schedule_id = 804, December binospec run
# schedule_id = 694, Binospec Commissioning Nov Run
#
#  Also, "schedule_id" == "queue_id" in the ObservatoryManager GUI code.
schedule_id = 813

# This is a global dictionary of Astroplan FixedTarget's, indexed by 'objid_id' from the 'fields' structure.
fixed_targets = {}

# Target Of Opportunity info as a mutable object.  This allows new values to be passed to
# constraints for processing.  If it were immutable, we could not programmatically update values
# during the scheduling process.
ToO = {}

# These are elements of the target-of-opportunity structure:
#    "found_ToO" is a boolean on whether there is at least one ToO observable for the scheduling time
#    "objid_ids" are the target (i.e., object) ID's supplied by the user.
#    "all_ToO" are the FixedTarget objects corresponding to the objid_ids.
#    "observable_ToO" are the FixedTarget objects that are observable at a specific scheduling time.

ToO['found_ToO'] = False
ToO['objid_ids'] = []
ToO['all_ToO'] = []
ToO['observable_ToO'] = []

# The sqlite3 database file for caching.
sqlite_file = 'mmtqueue.sqlite'

# Global hash of events such as times for sunrise, sunset, target_set, etc.
use_daily_events = True
daily_events = {}
# Tables within the sqlite3 database:
#
if False:
    print("sqlite_file", sqlite_file)
# The table name for cached constraint scores.
table = "scores"

# Default schema for schedules table.
"""
CREATE TABLE `schedules` (
	`key`	TEXT,
	`value`	TEXT,
	`id`	INTEGER PRIMARY KEY AUTOINCREMENT,
	`timestamp`	TEXT DEFAULT CURRENT_TIMESTAMP
);
"""

# Default schema for scores table
"""
CREATE TABLE `scores` (
	`key`	TEXT,
	`value`	REAL,
	`schedule_MST`	TEXT,
	`schedule_id`	INTEGER,
	`comment`	TEXT,
	PRIMARY KEY(`key`)
);
"""

# To empty the scores table if it gets too big.
"""
DELETE from scores;
vacuum scores;
"""

# Create a shared global database connection
conn = sqlite3.connect(sqlite_file)

# For debugging sqlite.
# cursor = conn.cursor()
# cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
# print(cursor.fetchall())
# if True:
#     print("conn: ", repr(conn))

# Setting query results to be a standard Python dictionary.
conn.row_factory = dict_factory

# This is used for the MaskNumberConstraint
# Always assume that masks 110 and 111 are on Binospec
default_masks = ['110','111']
# This keeps track of which masks have been used in scheduling.
# Used by MaskNumberConstraint
masks_used = default_masks

# Lookup tables.
program_id2program = {}
block_ids = {}
objid_ids = {}


def setup(args):
    """Retrieves configuration parameters from the scheduler server.

    Retrieves rows pertaining to the given keys from the Table instance
    represented by big_table.  Silly things may happen if
    other_silly_variable is not None.

    Args:
        args: An associative array from command line arguments (i.e, from argparse()).

    Returns:
        A boolean: True for success; False for failure

    Raises:
        None implemented yet.
    """

    global schedule_id
    global debug
    global verbose
    global instrument
    global meta
    global allocation
    global stats
    global configuration
    global configs
    global fields
    global start_date
    global end_date
    global mode 
    global version 
    global program_id2program 
    
    global block_ids
    global objid_ids
    
    # Overall approach for input parameters:
    # 1) Command-line arguments get the highest priority,
    # 2) Parameters from the configs dictionary are next highest priority.  These configs values are obtained by
    # a web query.
    # 3) Finally, fall back to a default value.
    # All parameters should be defined in the configs structure with command-line arguments over-ruling those values.
    if 'debug' in args and not args['debug'] is None:
        if args['debug'].lower() == "true":
            txt = "Turning debug on"
            debug = True
        else:
            txt = "Turning debug off"
            debug = False
        print(txt)
    else:
        debug = False

    # Code version.
    version = 2.0 
    # Note: We need to define instrument before printing anything or posting anything to Redis.
    if not args['version'] is None:
        version = float(args['version'])
        if verbose:
            txt = "Setting software version to: {0}".format(version)

    # Note:  schedule_id is also in the "meta" data structure.
    # We need to pass in schedule_id from command line so that
    # the correct "meta" data is returned from the database.
    if 'schedule_id' in args and not args['schedule_id'] is None:
        schedule_id  = int(args['schedule_id'])
        if verbose:
            txt = "Setting schedule_id from args array to: {0}".format(schedule_id)
            print(txt)   

    # Grab all of the "meta" data about this queue schedule.  
    meta = get_config_json(schedule_id)

    # Sending "which_instrument" to dataserver2.
    instrument = meta['instrument']
    cmd = r'echo "set json { \"which_instrument\" : \"' + instrument + r'\" }" | nc -w 1 tcs-server 8660'
    subprocess.call(cmd, shell=True)

    configuration = meta['configuration']

    # An associative array of configuration key:value pairs.
    configs = {}
    for conf in configuration:
        configs[conf] = configuration[conf]['parametervalue']

    # Two different modes to run mmtscheduler:
    # mode == 1, scheduler
    # mode == 2, dispatcher
    if 'mode' in args and not args['mode'] is None:
        mode  = int(args['mode'])
        if verbose:
            txt = "Setting mode to: {0}".format(mode)
            print(txt)      
    # Then, use the configs values.
    elif 'mode' in configs:
        mode = int(configs['mode'])
    # Should not get here...
    else:
        mode = 1
        
    if 'start_date' in args and not args['start_date'] is None:
        time = Time(args['start_date'])
        start_date = time
    elif 'start_date' in configs:
        time = Time(configs['start_date'])
        # Define the start of the scheduling run to the nearest sunrise.
        # self.start_time = self.observer.sun_set_time(time, which='nearest', horizon=-12*u.deg)
        # self.start_time = self.observer.sun_set_time(time, which='nearest', horizon=-6*u.deg)
        start_date = time
    # Should not get here...
    else:
        # Starting now.
        start_date = Time.now()
        
    if 'end_date' in args and not args['end_date'] is None:
        time = Time(args['end_date'])
        end_date = time
    elif 'end_date' in configs:
        time = Time(configs['end_date'])
        # Define the start of the scheduling run to the nearest sunrise.
        # self.start_time = self.observer.sun_set_time(time, which='nearest', horizon=-12*u.deg)
        # self.start_time = self.observer.sun_set_time(time, which='nearest', horizon=-6*u.deg)
        end_date = time
    # Should not get here...
    else:
        # Ending 12 hours from now.
        end_date = Time.now() +  TimeDelta(43200.0, format='sec')      
 
       
    if False:
        print("start_date: ", repr(start_date))
        print("end_date: ", repr(end_date))
    
    # This is the PostgreSQL database (on scheduler) key for the current queue.  All queues
    # have a unique ID.  This schedule ID is used to keep track of different queue runs.
    if 'schedule_id' in meta and not args['schedule_id'] is None:
        if schedule_id  == int(meta['schedule_id']):
            if False:
                print("The schedule_id agrees from both the command line and database: ", \
                      str(schedule_id))
        else:
            # Default to something.  This is the Feb 2018 Binospec run.
            # schedule_id = 812
            if True:
                print("There is a problem.  The command line schedule_id: ", \
                      str(schedule_id), ", does not match the database schedule_id:", \
                      str(meta['schedule_id']))


    # The instrument can currently be either 'binospec' or 'mmirs'
    # The instrument name is cast to all lower characters so 
    # capitalization doesn't matter.
    # Note that the instrument prefix for Redis is all caps:
    # either BINOSPEC or MMIRS.
    if 'instrument' in args and not args['instrument'] is None:
        if args['instrument'].lower() == 'mmirs':
            instrument = "mmirs"
        else:
            instrument = "binospec"           
    # Usually the instrument data is in "meta" unless you want to
    # really confuse things, e.g., running a Binospec queue with
    # MMIRS targets.
    elif 'verbose' in meta:
        if meta['instrument'].lower() == 'mmirs':
            instrument = "mmirs"
        else:
            instrument = "binospec" 
    else:
        instrument = "binospec" 
    
    # Dallan's ObservatoryManager code is providing this value.
    # Just keep track of it here so that we can send it back to
    # the ObservatoryManager after scheduling is done.
    if 'queuerun_id' in meta and not meta["queuerun_id"] is None:
        queuerun_id = meta["queuerun_id"]
    else:
        queuerun_id = "0"

    # Big substructure of "meta" of the allocation of time to different observing programs.
    allocation = meta['allocation']
    if False:
        txt = json.dumps(allocation, indent=2,ensure_ascii=False)
        print(txt)

    # Big substructure of meta of the allocation of stats for programs:
    # How much time has already be used for each program, etc.
    stats = meta['stats']
    if False:
        txt = json.dumps(stats, indent=2,ensure_ascii=False)
        print(txt)


    # Big substructure of meta with the observation fields details.
    # "fields" are really the same as "observing blocks".
    # Using 'fields' here for historical purposes.   
    fields = meta['fields']
    if True:
        txt = json.dumps(fields, indent=2,ensure_ascii=False)
        print(txt)

    # Populate lookup tables.
    for f in fields:
        program_id2program[f['program_id']] = f['program']  
        long_name = get_name(f)
        block_ids[long_name] = f['block_id']
        objid_ids[long_name] = f['objid_id']

    return True

def add_observer():
    """Adds an observer for the MMT Observatory to a scheduling calculation.

    Args:
        None

    Returns:
        Returns an Astroplan Observer instance for the MMT Observatory

    Raises:
        None.
    """
    global mmto
    global mmto_earth_location
    mmto = Observer(longitude=249.11499999999998*u.deg,
                                 latitude=31.688333333333333*u.deg, 
                                 elevation=2608*u.m,
                                 name="mmto",
                                 timezone="America/Phoenix")
    
    # Also, need an EarthLocation version of the MMT location
    mmto_earth_location = EarthLocation(lon=249.11499999999998*u.deg, \
                                lat=31.688333333333333*u.deg,  \
                                height=2608*u.m, )
    
    # times = Time(["2017-08-01 06:00", "2017-08-01 12:00", "2017-08-01 18:00"])

    if False:
        print(get_now(), mmto, mmto_earth_location)

    return mmto

# Read in the table of targets
# from astropy.io import ascii
# target_table = ascii.read('targets.txt')
#targets = [FixedTarget(coord=SkyCoord(ra=ra*u.deg, dec=dec*u.deg), name=name)
#           for name, ra, dec in target_table]

# Define a set a new Astroplan Constraints.
# These include:
#   1) MaskAngleConstraint
#   2) MaskNumberConstraint
#   3) MeridianConstraint
#   4) PIPriorityConstraint
#   5) RotatorConstraint
#   6) TimeAllocationConstraint
#   7) TimeConstraint (modified to handle multiple time interval)
#
#   The Constraints imported from Astroplan include:
#   1) AirmassConstraint
#   2) AltitudeConstraint
#   3) AtNightConstraint

# A total of ten Constraints.


# This is derived from the SequentialScheduler class from astroplan 0.4
# There are additional attributes and logging of results, among other
# changes.
class MMTSequentialScheduler(Scheduler):
    """
    A scheduler that does "stupid simple sequential scheduling".  That is, it
    simply looks at all the blocks, picks the best one, schedules it, and then
    moves on.
    """
    # The list "masks_used" can be local to this class,
    # The structure "stats" needs to be global since it is 
    # updated here, but is used by constraint calculations,
    # based upon those updated values. 
    def __init__(self, masks_used=None, \
                 default_masks=None, \
                 target_of_opportunity_ids='', \
                 *args, **kwargs):
        super(MMTSequentialScheduler, self).__init__(*args, **kwargs)
        # This class attribute keeps track of which masks have been
        # scheduling as the new schedule is computed.
        # Its initial value can either be an empty list or
        # a list with mask_id's [110, 111].  (110 == imaging,
        # 111 = Longslit1)  These two masks are normally on Binospec.
        # There may be a special case where them would be removed.
        # An emply masks_used list would be used in that case.
        self.masks_used = masks_used
        self.default_masks = default_masks
        self.target_of_opportunity_ids = [x.strip() for x in target_of_opportunity_ids.split(',')]


    # Function needs to accept a list.
    def _do_constraint(self, constraint, observer, target, times):
        # Note:  times must be iterable.
        return constraint(observer, target, times)

    def _make_schedule(self, blocks):
        global daily_events, ToO

        # By this time, we know if there are any targets-of-opportunity (ToO's) and the
        # FixedTargets for scheduling have been created.

        # So, loop through the objid_id's for any ToO's and add to ToO['all_ToO']
        for objid in ToO['objid_id']:
            t = fixed_targets[objid]
            # Prevent duplicate entries.
            if t not in ToO['all_ToO']:
                # Note that these are SkyCoord objects, not FixedTarget objects.
                ToO['all_ToO'].append(fixed_targets[objid].coord)

        pre_filled = np.array([[block.start_time, block.end_time] for
                               block in self.schedule.scheduled_blocks])
        if len(pre_filled) == 0:
            a = self.schedule.start_time
            filled_times = Time([a - 1*u.hour, a - 1*u.hour,
                                 a - 1*u.minute, a - 1*u.minute])
            pre_filled = filled_times.reshape((2, 2))
        else:
            filled_times = Time(pre_filled.flatten())
            pre_filled = filled_times.reshape((int(len(filled_times)/2), 2))
        
        for b in blocks:
            
            if b.constraints is None:
                b._all_constraints = self.constraints
            else:
                b._all_constraints = self.constraints + b.constraints
            # to make sure the scheduler has some constraint to work off of
            # and to prevent scheduling of targets below the horizon
            # TODO : change default constraints to [] and switch to append
            if b._all_constraints is None:
                b._all_constraints = [AltitudeConstraint(min=0 * u.deg)]
                b.constraints = [AltitudeConstraint(min=0 * u.deg)]
            elif not any(isinstance(c, AltitudeConstraint) for c in b._all_constraints):
                b._all_constraints.append(AltitudeConstraint(min=0 * u.deg))
                if b.constraints is None:
                    b.constraints = [AltitudeConstraint(min=0 * u.deg)]
                else:
                    b.constraints.append(AltitudeConstraint(min=0 * u.deg))
                    
            if test_MMTSequentialScheduler:
                print("All Constraints: ", repr(b._all_constraints))
                
            b._duration_offsets = u.Quantity([0*u.second, b.duration/2,
                                              b.duration])
            b.observer = self.observer
        current_time = self.schedule.start_time
        
        if test_MMTSequentialScheduler:
            print("Current time: ", repr(current_time))
            # print("Current time: ", repr(current_time), ", blocks:", repr(blocks))
        
        constraint_details = []
        
        while (len(blocks) > 0) and (current_time < self.schedule.end_time):

            t = current_time.datetime
            txt = str(Time.now()) + ", Begin scheduling for: " + str(t)
            to_stdout(txt, mode, redis_client, instrument, do_redis)
            
            # new code to skip through the daylight hours at the MMT...
            # if it's after 7AM MST, jump to the next sunset time (+5 minutes)
            if t.hour >= 14:
                if use_daily_events:
                    comment = 'next_sunset'
                    key = get_sun_moon_date_key(comment, Time(t))

                    if key in daily_events:
                        current_time = daily_events[key]
                        if False:
                            print("current_time: ", repr(current_time))
                    else:
                        current_time = self.observer.sun_set_time(current_time,
                                                                  which='next',
                                                                  horizon=configs['max_solar_altitude'] * u.deg)
                        daily_events[key] = current_time
                else:
                    current_time = self.observer.sun_set_time(current_time,
                                    which='next', 
                                    horizon=configs['max_solar_altitude']*u.deg)
                
                # This is resetting which masks could be used each night.
                # For Binospec, we can change the masks every day.
                # Even changing during the night is possible, but not
                # desirable.
                #
                # See if this works.  default_masks_used is a global.
                self.masks_used = self.default_masks
                ##### This code allows a gap to be inserted into a queue schedule. 

                #####!!!!!! This is a hard-coded hack.  If it's Feb 27, 2018,
                #####!!!!!!  advance scheduling to April 5, 2018
                if t.year == 2018 and t.month == 2 and t.day == 27:
                    # Advancing 39.5 days.  (Feb 26 to April 5)
                    current_time += 60*60*24*39.5*u.second
                    # Get the next sunset

                    if use_daily_events:
                        comment = 'next_sunset'
                        key = get_sun_moon_date_key(comment, Time(t))

                        if key in daily_events:
                            current_time = daily_events[key]
                            if False:
                                print("current_time: ", repr(current_time))
                        else:
                            current_time = self.observer.sun_set_time(current_time,
                                                                      which='next',
                                                                      horizon=configs['max_solar_altitude'] * u.deg)
                            daily_events[key] = current_timg
                    else:
                        current_time = self.observer.sun_set_time(current_time,
                                                                  which='next',
                                                                  horizon=configs['max_solar_altitude'] * u.deg)

                    # current_time = self.observer.sun_set_time(current_time, \
                    #                    which='next', \
                    #                    horizon=configs['max_solar_altitude']*u.deg)
                #####!!!!!! End of hard-coded hack.

                # Adding a five--minute "buffer" after sunset.
                # We need to do this so that we don't fail the "max_solar_altitude" criteria
                # at the beginning of the night.
                current_time += 300*u.second
                current_time.format = 'isot'
                txt = str(Time.now()) + ", Advancing scheduling to: " +                             str(current_time.datetime)
                to_stdout(txt, mode, redis_client, instrument, do_redis)
            
            # first compute the value of all the constraints for each block
            # given the current starting time
            if test_MMTSequentialScheduler:
                print("Current time: ", repr(current_time), ", blocks:", repr(blocks))
            
            block_transitions = []
            block_scores = []
            block_constraint_results = []
            constraint_scores = {}
            target_constraint_scores = {}

            #
            # New:  Target of Opportunity (ToO's) code section
            #
            # Pre-processing observing blocks for any ToO's:
            #
            if True:

                ToO['observable_ToO'] = []

                # Passing in the target_of_opportunity_names through a configuration
                # parameter.  This CSV string is parsed during initialization.
                if False:
                    print("targets of opportunity target ids: ", self.target_of_opportunity_names)
                if len(ToO['all_ToO']) > 0:
                    for b in blocks:
                        name = get_target_name(b)
                        if False:
                            print("Name: ", name)

                        ToO['found_ToO'] = False
                        found = False
                        for skycoord in ToO['all_ToO']:
                            if b.target.coord == skycoord:
                                target = b.target

                                # Calculate the time duration of the block.
                                #
                                # First figure out the transition
                                #
                                # We are not using transition blocks at the moment so "trans" should
                                # be None and the transition_time will be 0 seconds.
                                # Adding the code in case we start using transition blocks.
                                if len(self.schedule.observing_blocks) > 0:
                                    trans = self.transitioner(
                                        self.schedule.observing_blocks[-1], b, current_time, self.observer)
                                else:
                                    trans = None
                                transition_time = 0*u.second if trans is None else trans.duration

                                times = current_time + transition_time + b._duration_offsets
                                # If the block is a ToO and is always observable for the time duration,
                                # then we want to schedule it.
                                #
                                # A zero will be set to any non-ToO blocks in the next section of code.
                                #
                                if np.all(is_always_observable(b._all_constraints, self.observer, [target], times)):
                                    ToO['found_ToO'] = True
                                    # Note that which_ToO is a SkyCoord to be compatible with constraints,
                                    # not a FixedTarget or an observing block.
                                    if skycoord not in ToO['observable_ToO']:
                                        ToO['observable_ToO'].append(skycoord)
                                    if False:
                                        print("1421: ToO", repr(ToO))
                                        print("break")

            if False:
                print("3643 ToO: ", repr(ToO))

            # This is where we would parallelize for multiple blocks.
            for b in blocks:
                # first figure out the transition
                if len(self.schedule.observing_blocks) > 0:
                    trans = self.transitioner(
                        self.schedule.observing_blocks[-1], b, current_time, self.observer)
                else:
                    trans = None
                block_transitions.append(trans)
                transition_time = 0*u.second if trans is None else trans.duration

                times = current_time + transition_time + b._duration_offsets

                # make sure it isn't in a pre-filled slot
                if (any((current_time < filled_times) & (filled_times < times[2])) or
                        any(abs(pre_filled.T[0]-current_time) < 1*u.second)):
                    block_constraint_results.append(0)

                else:
                    constraint_res = []

                    #
                    # Target of Opportunity (ToO) calculations:
                    #
                    # If we found a ToO and if the target for the block is 
                    # the ToO, append a constraint score result of 1.  
                    # Otherwise, append a constraint score of 0.  
                    # This latter case will cause the overall_score to be zero for all non-ToO blocks.
                    #
                    # Code needs testing.
                    if False:
                        print("15352, ToO: ", repr(ToO))
                        if ToO['found_ToO']:
                            # See if the target is a ToO.
                            if b.target.coord is ToO['observable_ToO']:
                                # Appending a 1, which has no overall effect when the
                                # constraint results are multiplied for the overall score.
                                # Note that several OB may have a score of 1 here if they
                                # are for the same ToO target.
                                # Other constraints will rank blocks with the same ToO target.
                                if False:
                                    print("constraint_res append 1")
                                constraint_res.append(1.0)
                            else:
                                # Appending a zero to the constraint results list, this will
                                # cause the overall score to be zero when all results are
                                # multiplied together.
                                if False:
                                    print("constraint_res append 0")
                                constraint_res.append(0.1)
                        else:
                            # Doing this for consistency so that each block will have the same
                            # number of constraint results to be multiplied for the overall score.
                            if False:
                                print("constraint_res append 1")
                            constraint_res.append(1.0)

                    for constraint in b._all_constraints:
                        key = get_key(constraint, b.target, times[0], time2=times[-1])
                        if False:
                            # This "if" branch:
                            #    1) does not use the sqlite database for caching
                            c = constraint(self.observer, b.target, times)
                            score = np.prod(c)
                        else:
                            obj = get_score(conn, key)
                            # These astronomical constraints will not change with programmatic changes,
                            # so we can use cached values.
                            if ("RotatorConstraint" in key or \
                                "TargetOfOpportunityConstraint" in key or \
                                "AirmassConstraint" in key or \
                                "AltitudeConstraint" in key or \
                                "MeridianConstraint" in key or \
                                "AtNightConstraint" in key or \
                                "MoonSeparationConstraint" in key) and \
                                obj is not None:
                              
                                    # Not using stored scores for now.  We may be able to
                                    # use some scores, such as MeridianConstraint, but
                                    # probably not programmatic constraints, such as
                                    # TimeAllocationConstraint.  Programmatic constraints
                                    # should be recalculated in case something has changed.
                                    if False:
                                        print("obj: ", repr(obj))
                                    score = float(obj['value'])
                                    if schedule_id != obj['schedule_id']:
                                        # Updating this score.  Commonly, the current 
                                        # schedule_id will be newer than the saved score.
                                        # We've calculated this target/time/constraint
                                        # combination before for an older schedule.
                                        # We want to keep the schedule_id with
                                        # the most recent schedule.
                                        schedule_MST = obj['schedule_MST']
                                        set_score(conn,key,score,time=schedule_MST,schedule_id=schedule_id)
                                    if False:
                                        print("key:", key, ", using saved score:", score)
                            else:
                                if True:
                                    try:
                                        c = constraint(self.observer, b.target, times)
                                        score = np.prod(c)
                                    except:
                                        print("traceback: ", traceback.format_exc())
                                        # We shouldn't get here.  Something is wrong with the input values.
                                        print("Cannot calculate constraint score for ", get_constraint_name(constraint))

                                schedule_time = localize_time(times[0])
                                set_score(conn,key,score,time=schedule_time,schedule_id=schedule_id)
                            if verbose3:
                                print(key, ": ", score)    

                        constraint_res.append(score)
                            
                        target_name = get_target_name(b)
                        constraint_name = get_constraint_name(constraint)
                        if target_name not in target_constraint_scores:
                            target_constraint_scores[target_name] = {}

                        if False:
                            print("target_name: ", target_name, ", constraint_name: ", constraint_name, ", score: ", str(score))

                        target_constraint_scores[target_name][constraint_name] = score

                        schedule_time = localize_time(times[0])
                        set_score(conn,key,score,time=schedule_time,schedule_id=schedule_id)
                                
                    # take the product over all the constraints *and* times
                    overall_score = np.prod(constraint_res)

                    # Saving overall score for this target to sqlite.
                    target_key = get_target_key(target_name, times[0], time2=times[-1])
                    schedule_time = localize_time(times[0])
                    set_score(conn,target_key,overall_score,time=schedule_time,schedule_id=schedule_id)

                    block_constraint_results.append(overall_score)
                    # Need the block to calculated RA, Dec, etc. later.
                    block_scores.append({"overall_score":overall_score,"name":b.target.name,"block":b})
                    constraint_scores[b.target.name] = overall_score
                    if verbose:
                        txt = "***** Target: " + b.target.name.ljust(25) +  " overall_score: " + \
                                str(round(overall_score,8)).ljust(12) + " *****"
                        to_stdout(txt, mode, redis_client, instrument, do_redis)

            # now identify the block that's the best
            bestblock_idx = np.argmax(block_constraint_results)
            txt = "***** Best result: " + str(block_constraint_results[bestblock_idx]) + " *****"
            to_stdout(txt, mode, redis_client, instrument, do_redis)

                        
            if block_constraint_results[bestblock_idx] == 0.:
                # if even the best is unobservable, we need a gap
                current_time += self.gap_time
                txt = "Adding gap_time", str(self.gap_time)
                to_stdout(txt, mode, redis_client, instrument, do_redis)
            else:
                # If there's a best one that's observable, first get its transition
                trans = block_transitions.pop(bestblock_idx)
                if trans is not None:
                    self.schedule.insert_slot(trans.start_time, trans)
                    current_time += trans.duration

                # now assign the block itself times and add it to the schedule
                newb = blocks.pop(bestblock_idx)
                newb.start_time = current_time
                current_time += newb.duration
                newb.end_time = current_time
                newb.constraints_value = block_constraint_results[bestblock_idx]

                # This is a modification from the previous version of the code.
                
                if current_time > self.schedule.end_time:
                    # We need to append the block back onto the blocks list
                    # since it was popped off and we didn't schedule it.
                    blocks.append(newb)
                    # Do we need to set the current time back as well???
                    # Will this put us into an endless loop of popping and appending the same block???
                    # current_time -= newb.duration
                    # We can't insert this block, so bail out.
                    break               
                
                else:
                    self.schedule.insert_slot(newb.start_time, newb)

                    # This is where new code starts.  We need to publish
                    # results of scheduling while a new schedule is being
                    # computed.

                    # Inserting print information on block being inserted.  JDG 2016-11-30
                    txt = str(Time.now())  +", Scheduling " + newb.target.name + \
                            " from " + str(newb.start_time) + " to " + str(newb.end_time)
                    to_stdout(txt, mode, redis_client, instrument, do_redis)   
                    txt = ""
                    to_stdout(txt, mode, redis_client, instrument, do_redis)
            
                    # Adding this newblock into time used for the correct program.
                    """
                    "program_hours_allocated": 5.99,
                    "total_hours_requested": 7.5,
                    "total_hours_used": 0,
                    """
                    # name = newb.target.name
                    name = get_target_name(newb)
                    program_id = get_program_id(name)
                    program = program_id2program[program_id]
                    dt = (newb.end_time.unix - newb.start_time.unix)/3600.0
                    if test_MMTSequentialScheduler:
                        print("name: ", name)
                        print("program_id: ", program_id)
                        print("dt: ", dt)
                        # This is before (i.e., "pre") adding in the time for the block

                        print("total_hours_used (pre): ", stats[program]['total_hours_used'])

                    # Add time to stats for program (n hours)
                    stats[program]['total_hours_used'] += dt
                    if test_MMTSequentialScheduler:
                        # This is after (i.e., "post") adding in the time for the block
                        print("total_hours_used (post): ", stats[program]['total_hours_used'])

                    # Add the mask_id associated this block into "masks_used"
                    # The global masks_used list is used by MaskNumberConstraint
                    mask_id = get_mask_id(name)
                    if mask_id not in self.masks_used:
                        if test_MMTSequentialScheduler:
                                print("Appending to masks_used: ", mask_id, repr(masks_used))
                        self.masks_used.append(mask_id)

            
            sorted_block_scores = sorted(block_scores, key=lambda x: x['overall_score'], reverse=True)
            dispatcher_list = []
            rank = 0
            for bs in sorted_block_scores:
                rank += 1
                name = bs['name']
                block = bs['block']
                target = block.target
                
                # Search through the fields to find which field is associated
                # with this observing block.
                field = None
                for f in fields:
                    # HERE2 !!!
                    # if f['block_id'] + '_' + f['objid'] == name:
                    # my_name = f['objid'] + '_' +f['block_id'] + "_P" + str(f['pi_priority'])
                    my_name = get_name(f)
                    if my_name == name:
                        # print ("my_name: ", my_name)
                        field = f
                        break
                  
                
                altaz = target.coord.transform_to(AltAz(obstime=Time(t),location=mmto_earth_location))
                ra_dec = target.coord.to_string("hmsdms").split()
                ra = ra_dec[0]
                dec = ra_dec[1]
             
                # With the short-circuit of calculations, some values may not be defined.
                # Force them to be None.
                my_keys = ["TimeAllocationConstraint",
                           "MeridianConstraint",
                           "PIPriorityConstraint",
                           "AirmassConstraint",
                           "AltitudeConstraint",
                           "AtNightConstraint",
                           'MaskAngleConstraint',
                           "MoonSeparationConstraint",
                           "RotatorConstraint",
                           "TimeConstraint"]

                for k in my_keys:
                    if k not in target_constraint_scores[name]:
                        target_constraint_scores[name][k] = '' 
                

                debug_json = False
                # Build up a nested structure for JSON and Redis.
                data = {"timestamp":int(time.time()),                   
                        "update":time.strftime('%Y-%m-%d %H:%M:%S') + " UTC"}
                
                if rank is None:
                    print('rank failed!!!')
                    data['rank'] = 'Undefined'
                else:
                    if debug_json:
                        print("rank: ", rank)
                    data['rank'] = rank
                    
                if name is None:
                    print('name failed!!!')
                    data['name'] = 'Undefined'
                else:
                    if debug_json:
                        print("name: ", name)
                    data['name'] = name
                
                if t is None:
                    print('start_time failed!!!')
                    data['start_time'] = 'Undefined'
                else:
                    x = str(t) + " UTC"
                    if debug_json:
                        print("t: ", x)
                    data['start_time'] = x
                
                if ra is None:
                    print("ra failed!!!") 
                    data['ra'] = 'Undefined'
                else:
                    if debug_json:
                        print("ra: ", ra)
                    data['ra'] = ra
                    
                if dec is None:
                    print("dec failed!!!") 
                    data['dec'] = 'Undefined'
                else:
                    if debug_json:
                        print("dec: ", dec)
                    data['dec'] = dec
                    
                if 'pi_priority' in field:
                    if debug_json:
                        print("pi_priority: ", field['pi_priority'])
                    data['pi_priority'] = field['pi_priority']
                else:
                    print("pi_priority failed!!!") 
                    data['pi_priority'] = 'Undefined'
                
                if 'tac_priority' in field:
                    if debug_json:
                        print("tac_priority: ", field['tac_priority'])
                    data['tac_priority'] = field['tac_priority']
                else:
                    print("tac_priority failed!!!")     
                    data['tac_priority'] = 'Undefined'
                
                 # duration is in minutes
                if 'duration' in field:
                    x = round(float(field['duration'])/60.0,2)
                    if debug_json:
                        print("duration: ", x)
                    data['duration'] = x
                else:
                    print("duration failed!!!")     
                    data['duration'] = 'Undefined'
                    
                if 'block_id' in field:
                    x = field['block_id']
                    if debug_json:
                        print("block_id: ", x)
                    data['block_id'] = x
                else:
                    print("block_id failed!!!")     
                    data['block_id'] = 'Undefined'
                
                if 'objid_id' in field:
                    x = field['objid_id']
                    if debug_json:
                        print("objid_id: ", x)
                    data['objid_id'] = x
                else:
                    print("objid_id failed!!!")     
                    data['objid_id'] = 'Undefined'
                
                if altaz is None:
                    print("alt failed!!!") 
                    data['alt'] = 'Undefined'
                else:
                    x = str(altaz.alt)
                    if debug_json:
                        print("alt: ", x)
                    data['alt'] = x
             
                if altaz is None:
                    print("az failed!!!") 
                    data['az'] = 'Undefined'
                else:
                    x = str(altaz.az)
                    if debug_json:
                        print("az: ", x)
                    data['az'] = x       
                
                for key in my_keys:
                    if name != "" and \
                            name in target_constraint_scores and \
                            key in target_constraint_scores[name]:
                        x = target_constraint_scores[name][key]
                        if len(str(x)) == 0:
                            x = 'Undefined'
                        # Ugh!  numpy.int64 are not serializable. Go figure...
                        # Casting it as a string.
                        data[key] = str(x)
                        if debug_json:
                            print("{}: {}".format(key, x))
                    else:
                        print("{} failed!!!".format(key)) 
                        data[key] = 'Undefined'
                # "overall_score":round(bs['overall_score'],8) }
                
                key = 'overall_score'
                if key in bs:
                    x = round(bs[key],8)
                    data[key] = x
                    if debug_json:
                        print("{}: {}".format(key, x))
                else:
                    print("{} failed!!!".format(key)) 
                    data[key] = 'Undefined'
                    
                """
                data = {"rank":rank,                   
                        "timestamp":int(time.time()),                   
                        "update":time.strftime('%Y-%m-%d %H:%M:%S') + " UTC",                   
                        "name":name,                   
                        "start_time":str(t) + " UTC",                   
                        "ra":ra,                   
                        "dec":dec,                   
                        "pi_priority":field['pi_priority'],                   
                        "tac_priority":field['tac_priority'],  
                        # duration is in minutes
                        "duration":round(float(field['duration'])/60.0,2),
                        # "overhead":field['overhead'],
                        "block_id":field['block_id'],
                        "objid_id":field['objid_id'],
                        "alt":str(altaz.alt), \
                        "az":str(altaz.az), \
                        "TimeAllocationConstraint": target_constraint_scores[name]["TimeAllocationConstraint"], \
                        "MeridianConstraint": target_constraint_scores[name]["MeridianConstraint"], \
                        "PIPriorityConstraint": target_constraint_scores[name]["PIPriorityConstraint"], \
                        # "TACPriorityConstraint": scheduler.targets[name].constraint_scores["TACPriorityConstraint"], \
                        "AirmassConstraint": target_constraint_scores[name]["AirmassConstraint"], \
                        "AltitudeConstraint": target_constraint_scores[name]["AltitudeConstraint"], \
                        "AtNightConstraint": target_constraint_scores[name]["AtNightConstraint"], \
                        "MoonSeparationConstraint": target_constraint_scores[name]["MoonSeparationConstraint"], \
                        # "MaskAngleConstraint": target_constraint_scores[name]["MaskAngleConstraint"], \
                        # "MaskNumberConstraint": target_constraint_scores[name]["MaskNumberConstraint"], \
                        "RotatorConstraint": target_constraint_scores[name]["RotatorConstraint"], \
                        "TimeConstraint": target_constraint_scores[name]["TimeConstraint"], \
                        "overall_score":round(bs['overall_score'],8) }
                """
                dispatcher_list.append(data)
            
            # Add this "dispatcher" constraint structure to the growing list.
            constraint_details.append(dispatcher_list)

            # Done if Dispatcher mode
            if mode == 2:
                txt = get_now() + " Pushing dispatcher results to Redis"
                to_stdout(txt, mode, redis_client, instrument, do_redis)
                # to_output(conn, tbl, mode, redis_client, instrument, schedule_id, do_redis)
                to_status('stopped', mode, redis_client, instrument, do_redis)
                txt = get_now() + " Dispatcher stopped"
                to_stdout(txt, mode, redis_client, instrument, do_redis)
                # Trigger logging the dispatcher output to MySQL
                to_mysql(schedule_id, mode, redis_client, instrument, do_redis)
                to_output(conn, dispatcher_list, mode, redis_client, \
                         instrument, schedule_id, do_redis)
                break
                
        # Scheduler mode
        if mode == 1:
            txt = "Pushing constraint time series results to Redis"
            to_stdout(txt, mode, redis_client, instrument, do_redis)
            to_constraint_details(constraint_details, mode, redis_client, instrument, do_redis)
        
        return self.schedule



def add_global_constraints():
    global global_constraints
    #
    # Set up the three global Constraints:
    # 1) AirmassConstraint 
    # 2) AltitudeConstraint
    # 3) AtNightConstraint
    # 4) MoonSeparationConstraint
    #
    # All are used as booleans (0 or 1) in that
    # they either completely pass or fail the constraint.
    # Note that constraints are evaluated at the beginning,
    # middle, and end of each time block.  If a boolean
    # constraint fails any one of these time periods, 
    # it fails overall.
    #

    # Create the list of constraints that all targets must satisfy
    global_constraints = []
    
    #
    # Evaluating AirmassConstraint
    if use_AirmassConstraint and "use_airmass_constraint" in configs and \
                configs['use_airmass_constraint'].lower() == "true" and \
                is_numeric(configs['max_airmass']):
            # If True, return the score as either 1.0 or 0.0, 
            # not as as float.
            #
            # Just assume this is True for now.  We could add this
            # as a configs parameter in the future if we want.
            boolean_constraint = True
            c = AirmassConstraint(max = float(configs['max_airmass']), 
                                    boolean_constraint = boolean_constraint)
            if test_AirmassConstraint:
                print("Adding global constraint: ", repr(c))
            global_constraints.append(c)
    
   
    # 
    # Evaluating AltitudeConstraint
    if use_AltitudeConstraint and "use_altitude_constraint" in configs and \
                configs['use_altitude_constraint'].lower() == "true" and \
                is_numeric(configs['max_alt_degrees']) and \
                is_numeric(configs['min_alt_degrees']) and  \
                'use_altitude_boolean' in configs:
            # Forcing it to be True since we are using the MeridianConstarint
            # as well.
            if True or configs['use_altitude_boolean'].lower() == 'true':
                boolean_constraint = True
            else:
                boolean_constraint = False
            c = AltitudeConstraint(min = configs['min_alt_degrees']*u.deg,
                                    max = configs['max_alt_degrees']*u.deg,
                                    boolean_constraint = boolean_constraint)
            if test_AltitudeConstraint:
                print("Adding global constraint: ", repr(c))
            global_constraints.append(c)
    

    #
    # Evaluating AtNightConstraint
    if use_AtNightConstraint and "use_at_night_constraint" in configs and \
                configs['use_at_night_constraint'].lower() == "true" and \
                'max_solar_altitude' in configs and \
                is_numeric(configs['max_solar_altitude']):
            if float(configs['max_solar_altitude']) == -6:
                c = AtNightConstraint.twilight_civil()
            elif float(configs['max_solar_altitude']) == -12:
                c = AtNightConstraint.twilight_civil()
            else:
                c = AtNightConstraint.twilight_astronomical()
            """
            'max_solar_altitude' in configs and \
            is_numeric(configs['max_solar_altitude']):
                pass
    
                # Only allowing -6 and -12, otherwise use astronomical twilight.
                # We could use max_solar_altitude directly
                if 'max_solar_altitude' in configs and \
                is_numeric(configs['max_solar_altitude']) and \
                float(configs['max_solar_altitude']) == -6:
                    c = AtNightConstraint.twilight_civil()
                elif 'max_solar_altitude' in configs and \ 
                    is_numeric(configs['max_solar_altitude']) and \
                    float(configs['max_solar_altitude']) == -12:
                    c = AtNightConstraint.twilight_nautical()
                else:
                    c = AtNightConstraint.twilight_astronomical()
            """
            if test_AtNightConstraint:
                print("Adding global constraint: ", repr(c))
            global_constraints.append(c)
    

    # Evaluating MoonSeparationConstraint
    if use_MoonSeparationConstraint and "use_moon_separation_constraint" in configs and \
                configs['use_moon_separation_constraint'].lower() == "true" and \
                is_numeric(configs['moon_separation_degrees']):
            # If True, return the score as either 1.0 or 0.0, 
            # not as as float.
            #
            # Just assume this is True for now.  We could add this
            # as a configs parameter in the future if we want.
            boolean_constraint = True
            
            c = MoonSeparationConstraint(min=float(configs['moon_separation_degrees'])*u.deg)
            if test_MoonSeparationConstraint:
                print("Adding global constraint: ", repr(c))
            global_constraints.append(c)


# In[ ]:


def add_block_constraints():
    #
    # Now we want to process each of the "fields" (an
    # alias for observing block).  The rest of the 
    # constraints will be block-specific.  We use this
    # approach to pass in block-specific parameters
    # to the constraint that are used when the score for
    # the constraint is calculated.
    #
    global blocks
    global fixed_targets
    global target_names
    global daily_events
    global ToO

    fixed_targets = {}
    blocks = []
    target_names = []
    for field in fields:
        # Create a name for this observing block.  
        # The name will contain information about the 
        # mask_id, the block_id, program, etc.
        # Encapsulating this information in the name
        # makes the name unique and also allows the
        # information to be extracted later.
        # Astroplan FixedTargets are commonly cast as
        # Astropy SkyCoord.  When this happens, we lose
        # any target-specific information, such as a 
        # mask_id.  This target-specific information
        # can be extracted from the field/observing
        # block name as needed.
        name = get_name(field)
        
        target_names.append(name)
        
        if False:
            print("name: ", name, "field: ", repr(field))
        r = field['ra_hms']
        d = field['dec_dms']
        id = field['objid_id']
        if False:
            print("ra: ", repr(r))
            print("dec: ", repr(d))
        #
        # Create a standard Astroplan FixedTarget instance.
        t = FixedTarget(coord=SkyCoord(ra=r, dec=d), name=name)
        fixed_targets[id] = t
        if 'time_resolution_seconds' in configs:
            time_resolution_seconds = int(configs['time_resolution_seconds'])
        else:
            time_resolution_seconds = '20'
        # duration

        # The total duration of a field equals the "duration"
        # plus the "overhead".  We may use Astroplan's 
        # TransitionBlocks at a future time.  At the moment,
        # we are using these fixed overheads.
        d = float(field['duration'])
        o = float(field['overhead'])
        # Total duration, "df":  duration plus overhead
        df = roundup(d + o, time_resolution_seconds)*u.second
        if False:
            print("Duration+overhead: ", repr(df))

        # A default priority for this block We're not 
        # using this for PriorityScheduling so it doesn't really
        # matter what the value is.
        p = 1.0

        # Now, add all of the block constraints as needed for the 
        # observing block/field.
        block_constraints = []

        # Evaluating MaskAngleConstraint
        if use_MaskAngleConstraint and 'design_parang' in field and 'max_mask_angle' in configs:

            if test_MaskAngleConstraint:
                print("Doing MaskAngleConstraint")

            try:
                design_parang = field['design_parang']*u.deg
            except:
                print("design_parang invalid: ", repr(field['design_parang']))
                # Defaults to 0.0 degrees.
                design_parang = 0.0*u.deg

            try:
                max_mask_angle = configs['max_mask_angle']*u.deg
            except:
                print("max_mask_angle invalid: ", repr(configs['max_mask_angle']))
                # Defaults to 30.0 degrees
                max_mask_angle = 30.0*u.deg

            if test_MaskAngleConstraint:
                print("Defining MaskAngleConstraint, p: ", design_parang, ", m: ", design_parang)
            c = MaskAngleConstraint(mmto,
                     design_parang=design_parang, 
                     max_mask_angle=max_mask_angle)

            if test_MaskAngleConstraint:
                print("Adding constraint: ", repr(c))

            block_constraints.append(c) 

        # Evaluating MaskNumberConstraint
        if use_MaskNumberConstraint and "mask_id" in field:
            c = MaskNumberConstraint(field['mask_id'], 
                    masks_used)
            if test_MaskNumberConstraint:
                print("Adding constraint: ", repr(c))
            block_constraints.append(c)

        #  Evaluating MaskReadyConstraint
        if use_MaskReadyConstraint and \
                "use_mask_ready_constraint" in configs and \
                configs['use_mask_ready_constraint'].lower() == "true" and \
                'mask_ready' in field:

            try:
                mask_ready = field['mask_ready'].lower()
                if mask_ready == 'null':
                    mr = None
                elif mask_ready == 'true':
                    mr = True
                else:
                    mr = False
            except:
                mr = False

            if mr is not None:
                c = MaskReadyConstraint(mask_ready=mr)
                if test_MaskReadyConstraint:
                    print("Adding constraint: ", repr(c))
                block_constraints.append(c)

        #  Evaluating MeridianConstraint    
        if use_MeridianConstraint and 'pi_priority' in field:
            try:
                pp = float(field['pi_priority'])
            except:
                print("pi_priority invalid: ", repr(field['pi_priority']))
                # Defaults to a PIPriority of 1.0, the highest value.
                # We're currently using priorities 1 (highest), 2 and 3 (lowest).
                pp = 1.0

            c = MeridianConstraint(df, pp, daily_events)
            if test_MeridianConstraint:
                print("Adding constraint: ", repr(c))
            block_constraints.append(c) 

        #  Evaluating PIPriorityConstraint
        if use_PIPriorityConstraint and \
                "use_pi_priority_constraint" in configs and \
                configs['use_pi_priority_constraint'].lower() == "true" and \
                'pi_priority' in field:

            try:
                pp = float(field['pi_priority'])
            except:
                print("pi_priority invalid: ", repr(field['pi_priority']))
                pp = 1.0

            c = PIPriorityConstraint(pi_priority=pp)
            if test_PIPriorityConstraint:
                print("Adding constraint: ", repr(c))
            block_constraints.append(c)

        #  Evaluating RotatorConstraint
        if False:
            print("Trying RotatorConstraint")
            print("use_RotatorConstraint", repr(use_RotatorConstraint))
            if  "use_rotator_constraint" in configs:
                print("use_rotator_constraint is in configs")
            else:
                print("use_rotator_constraint is NOT in configs")
            
            print("use_rotator_constraint", repr(configs['use_rotator_constraint'].lower()))
            
            if "use_rotator_constraint" in configs:
                "use_rotator_constraint is in configs"
            else:
                "use_rotator_constraint is NOT in configs"
            
            if "max_rot_degrees" in configs:
                print("max_rot_degrees is in configs")
            else:
                print("max_rot_degrees is NOT in configs")
            
            if "min_rot_degrees" in configs:
                print("min_rot_degrees is in configs")
            else:
                print("min_rot_degrees is NOT in configs")
                 
            if "posang" in field:
                print("posang is in field")
            else:
                print("posang is NOT in field")
            
        
        if use_RotatorConstraint and \
            "use_rotator_constraint" in configs and \
            configs['use_rotator_constraint'].lower() == "true" and \
            'max_rot_degrees' in configs and \
            'min_rot_degrees' in configs and  \
            'posang' in field:

            if False:
                print("Doing RotatorConstraint")

            try:
                max_rot_degrees = float(configs['max_rot_degrees'])*u.deg
            except:
                print("max_rot_degrees invalid: ", repr(configs['max_rot_degrees']))
                # Defaults to 160 degrees. This should be safe for all
                # instruments.
                max_rot_degrees = 160.0*u.deg

            try:
                min_rot_degrees = float(configs['min_rot_degrees'])*u.deg
            except:
                print("min_rot_degrees invalid: ", repr(configs['min_rot_degrees']))
                # Defaults to -160 degrees. This should be safe for all
                # instruments.
                min_rot_degrees = -160.0*u.deg

            try:
                posang = float(field['posang'])*u.deg
            except:
                print("posang invalid: ", repr(field['posang']))
                # Position angle defaults to 0.0(degrees)
                posang = 0.0*u.deg

            c = RotatorConstraint( max=max_rot_degrees,
                     min=min_rot_degrees,
                     posang=posang)

            if test_RotatorConstraint:
                print("Adding constraint: ", repr(c))
            block_constraints.append(c)

        #  Evaluating ObjectTypeConstraint
        if use_ObjectTypeConstraint and \
                "use_object_type_constraint" in configs and \
                configs['use_object_type_constraint'].lower() == "true" and \
                'objtype' in field:

            try:
                object_type = field['object_type'].lower()
                c = ObjectTypeConstraint(object_type=object_type)
                if test_object_type_contraint:
                    print("Adding constraint: ", repr(c))
                block_constraints.append(c)
            except:
                print("ObjectTypeConstraint definition error")

        #  Evaluating TargetOfOpportunityConstraint
        if use_TargetOfOpportunityConstraint and \
                "use_target_of_opportunity_constraint" in configs and \
                configs['use_target_of_opportunity_constraint'].lower() == "true":


            if False:
                print("+++++++In TargetOfOpportunityConstraint+++++++")
                print("found_ToO: ", ToO['found_ToO'])
                print("which_ToO: ", ToO['which_ToO'])

            try:
                c = TargetOfOpportunityConstraint(ToO=ToO)
                if test_TargetOfOpportunityConstraint:
                    print("Adding constraint: ", repr(c))
                block_constraints.append(c)
            except:
                print("TargetOfOpportunityConstraint definition error")

              
        #  Evaluating TimeAllocationConstraint
        if use_TimeAllocationConstraint and  \
            "use_time_allocation_constraint" in configs and \
            configs['use_time_allocation_constraint'].lower() == "true" and                 'program' in field: 

            p = field['program']
            c = TimeAllocationConstraint(p,
                 stats)

            if test_TimeAllocationConstraint:
                print("constraint:", repr(c))
            block_constraints.append(c)

        # Evaluating TimeConstraint
        if use_TimeConstraint and \
            "use_time_constraint" in configs and \
            configs['use_time_constraint'].lower() == "true":

            # Handle case where a list of time_constraints are supplied.
            if 'time_constraints' in field:
                block_times = []
                for [t1,t2] in field['time_constraints']:
                    arr = [Time(t1),Time(t2)]
                    block_times.append(arr)
                    c = (TimeConstraint(None,None,block_times))
                    block_constraints.append(c)
            # Else, case where we have just one time constraint with
            # a single start and end.
            elif 'time_constraint_start' in field and 'time_constraint_end' in field:
                c = TimeConstraint(Time(field['time_constraint_start']),
                                        Time(field['time_constraint_end']),
                                        None)
                block_constraints.append(c)

        b = ObservingBlock(t,df,p,constraints=block_constraints)
        if False:
            print("Appending block: ", repr(b))
            print("block_constraints: ", repr(block_constraints))
        blocks.append(b)

def create_transitioner():
    global transitioner
    if 'slew_rate' in configs:
        try:
            slew_rate = float(configs['slew_rate'])* u.deg / u.second
        except:
            print("slew_rate invalid: ", crepr(configs['slew_rate']))
            slew_rate = 1.0* u.deg / u.second
    else:
        slew_rate = 1.0* u.deg / u.second

    transitioner = Transitioner(slew_rate,
                                 {'filter':{('B','G'): 10*u.second,
                                            ('te_schG','R'): 10*u.second,
                                            'default': 30*u.second}})




def create_scheduler():
    global seq_scheduler
    global masks_used
    global configs
    names = ''
    if 'target_of_opportunity_names' in configs:
        names = configs['target_of_opportunity_names']
    seq_scheduler = MMTSequentialScheduler(constraints = global_constraints,
                        observer = mmto,
                        transitioner = transitioner,
                        masks_used=masks_used,
                        default_masks=default_masks,
                        target_of_opportunity_names=names)
    if False: 
        print("seq_scheduler: ", repr(seq_scheduler))


# Testing a PriorityScheduler
def create_pscheduler():
    global pri_scheduler
    global masks_used
    pri_scheduler = PriorityScheduler(constraints = global_constraints,
                                           observer = mmto,
                                           transitioner = transitioner)

    if True:
        print("Here we go doing PriorityScheduler!!!")
        print("pri_scheduler: ", repr(pri_scheduler))

def run_pscheduler():
    global priority_schedule
    global target_names
    global pri_scheduler

    st = datetime.datetime.now()
    print("PriorityScheduler start time: ", st)

    priority_schedule = Schedule(start_date, end_date)
    if True:
        print(get_now(), "priority_schedule: ", repr(priority_schedule))

    s = pri_scheduler(blocks, priority_schedule)
        
    # print("Priority table: ", priority_schedule.to_table())
    print("Observing blocks:")
    for b in priority_schedule.observing_blocks:
        print("    block: ", repr(b))
    print("Scheduled blocks:")
    for b in priority_schedule.scheduled_blocks:
        print("    block: ", repr(b))
    print("Open slots:")
    for b in priority_schedule.open_slots:
        print("    slot: ", repr(b))

    et = datetime.datetime.now()
    print("PriorityScheduler end time: ", et)

    # Elapsed time to run the scheduler.
    time_diff = et - st
    txt = "Elapsed time: ", str(time_diff)
    print(txt)

def run_scheduler():
    global sequential_schedule
    global target_names
        
    st = datetime.datetime.now()
    print("start time: ", st)
  
    sequential_schedule = Schedule(start_date, end_date)
    if True: 
        print(get_now(), "sequential_schedule: ", repr(sequential_schedule))
  
    s = seq_scheduler(blocks, sequential_schedule)
    if False: 
        print(get_now(), "seq_scheduler: ", repr(s))

     
    # Porting from mmtscheduler.py
    txt = ""

    # Print out which targets are scheduled and not scheduled if in scheduler mode (== 1)
    if int(mode) == 1:
        try:
            tbl = table_to_entries(s.to_table(), block_ids, objid_ids)
        except Exception as e:
            print("table_to_entries error")
            print("traceback: " , traceback.format_exc())
            
        try:
            txt = json.dumps(tbl,ensure_ascii=False)
        except Exception as e:
            print("JSON error")
            print("traceback: " , traceback.format_exc())
        
        try:
            #
            # Do we only need to_output() to pass in the object and not txt???
            # 
            # to_output(txt, mode, redis_client, instrument, schedule_id, do_redis)
            to_output(conn, tbl, mode, redis_client, instrument, schedule_id, do_redis)

        except Exception as e:
            print("to_output error")
            print("traceback: " , traceback.format_exc())
            
            # to_output(scheduler.table_to_python(scheduler.queue))
            # Convert the queue to JSON and save to Redis.
        
        try:
            to_queue(txt, mode, redis_client, instrument, do_redis)
        except Exception as e:
            print("to_queue error")
            print("traceback: " , traceback.format_exc())            
        
        try:     
            to_stdout(txt, mode, redis_client, instrument, do_redis)
        except Exception as e:
            print("to_queue error")
            print("traceback: " , traceback.format_exc())  
            
        
        # Checking if all targets have been scheduled.
        # Perhaps make this a class method at some point.
        targets_scheduled = []
        targets_unscheduled = []

        for name in target_names:
            found = False
            for block in s.scheduled_blocks:
                try:
                    if block.target.name == name:
                        found = True
                except Exception as e:
                    pass
    
            if found:
                txt = "Target {0} is scheduled.".format(name)
                targets_scheduled.append(name)
            else:
                txt = "Target {0} is NOT scheduled.".format(name)
                targets_unscheduled.append(name)
            to_stdout(txt, mode, redis_client, instrument, do_redis)

        # Keeping Redis backward compatible for existing MMIRS redis parameters.
        if instrument == 'mmirs':
            to_redis('SCHEDULER.scheduled', targets_scheduled, redis_client, do_redis)
            to_redis('SCHEDULER.unscheduled', targets_unscheduled, redis_client, do_redis)
        # Setting either MMIRS or BINOSPEC.
        to_redis(add_prefix(instrument, 'SCHEDULER.scheduled'), targets_scheduled, \
                redis_client, do_redis)
        to_redis(add_prefix(instrument, 'SCHEDULER.unscheduled'), targets_unscheduled, 
                redis_client, do_redis)

    # Extracting details of allocated and completed times.
    program_ids = {}
    for field in fields:
        # program_id is the "name" of the program, e.g, SAO-4
        # program is the database number of the program, e.g., 459
        program_ids[field['program']] = field['program_id']
        if True:
            print("Adding program_id: ", field['program_id'], \
                  ", program: ", field['program'])

    allocated_time = {}
    completed_time = {}
    for program in stats:
        # "stats" is indexed on "program", e.t., 459.
        allocated_time[program] = stats[program]['program_hours_allocated']
        completed_time[program] = stats[program]['total_hours_used']
        if True:
            print("Working on program: ", program)
    
    to_allocated_time(allocated_time,  mode, redis_client, instrument, do_redis)
    to_completed_time(completed_time,  mode, redis_client, instrument, do_redis)
    
    # Print out the remaining allocated time for each program.
    txt = "Program".rjust(20) +  \
              "Scheduled Time (hrs)".rjust(25) +  \
              " Allocated Time (hrs)".rjust(25) +  \
              " Percent Scheduled (%)".rjust(25)
    to_stdout(txt, mode, redis_client, instrument, do_redis)
    txt = '-------------------------------------------------------------------------------------'
    to_stdout(txt, mode, redis_client, instrument, do_redis)
    for program in completed_time:
        # Checking for division-by-zero error.  JDG 2018-06-01
        if allocated_time[program] > 0.0:
            # Checking that there is at least a target in the program.  JDG 2018-06-01
            if program in program_ids:
                program_id = program_ids[program] 
                txt = program_id.rjust(20) + \
                    str(round(completed_time[program],2)).rjust(25) + \
                    str(round(allocated_time[program],2)).rjust(25) + \
                    str(round(100 * completed_time[program]/allocated_time[program], 2)).rjust(25)
                to_stdout(txt, mode, redis_client, instrument, do_redis)
    
    # Trigger logging the scheduler output to MySQL
    to_mysql(schedule_id, mode, redis_client, instrument, do_redis)

    # Scheduling is done.  Set SCHEDULER.status to 'stopped'.
    to_status('stopped', mode, redis_client, instrument, do_redis)
    to_stdout("Calculations complete", mode, redis_client, instrument, do_redis)

    et = datetime.datetime.now()
    print("end time: ", et)
    
    # Elapsed time to run the scheduler.
    time_diff = et - st
    txt = "Elapsed time: ", str(time_diff)
    print(txt)
    to_stdout(txt, mode, redis_client, instrument, do_redis)

# time_res in seconds, probably just let it default to an hour.
def show_scores(time_res=3600):
    global sequential_schedule
    global target_names

    # start_date and end_date should already be defined and available globally
    # my_schedule = Schedule(start_date, end_date)
    start_time = Time(roundHour(start_date.datetime))
    if False:
        print("Start time: ", repr(start_time))
    end_time = Time(roundHour(end_date.datetime))
    if False:
        print("End time: ", repr(end_time))
    my_schedule = Schedule(start_time, end_time)
    if False:
        print(get_now(), "my_schedule: ", repr(my_schedule))

    times = time_grid_from_range((start_time, end_time), time_res*u.second)
    
    # blocks should already be defined
    # global_constraints should already be defined and available globally
    for block in blocks:
        if True:
            print("Processing block: ", block.target.name)
        my_scorer = Scorer([block], mmto, my_schedule, global_constraints)
        if False:
            print(get_now(), "my_scorer: ", repr(my_scorer))

        # print("Creating score array:")
        score_array = my_scorer.create_score_array(time_resolution=time_res*u.second)
        if False:
            print(repr(score_array))
            
        
            
        if True:
            target = block.target
            for i in range(1,len(times)):
                try:
                    time = Time(roundTime(times[i].datetime))
                    local_time = localize_time(time)
                    value = score_array[0][i]
                    print("time: ", local_time, " MST, score: ", value)
                    # if value > 0.0:
                    key = get_short_key(target,time)
                    set_scorer(conn, key, value, time=localize_time(time), schedule_id=schedule_id, comment=get_program_id(key))
                except:
                    print("show_scores() traceback: " , traceback.format_exc())
  


"""  
This is the schedule for the two following time evaluations:
2018-01-27 10:07:34.649577 : scheduling  <astroplan.scheduling.ObservingBlock (stream_tri_03_7982_P1, 2017-12-15 00:51:00.000 to 2017-12-15 02:36:00.000) at 0x111ccdcc0>
2018-01-27 10:07:38.869270 : scheduling  <astroplan.scheduling.ObservingBlock (jj0151-2_7905_P1, 2017-12-15 02:36:19.166 to 2017-12-15 03:26:19.166) at 0x115574b38>
2018-01-27 10:07:43.259979 : scheduling  <astroplan.scheduling.ObservingBlock (PS16fgt (copy)_7753_P1, 2017-12-15 03:26:31.772 to 2017-12-15 03:46:31.772) at 0x114fe9eb8>
2018-01-27 10:07:47.742714 : scheduling  <astroplan.scheduling.ObservingBlock (jj0151-1_7904_P1, 2017-12-15 03:46:44.378 to 2017-12-15 06:16:44.378) at 0x115574940>
2018-01-27 10:07:51.891498 : scheduling  <astroplana.scheduling.ObservingBlock (SAO-9_L1527IRS (copy)_7874_P1, 2017-12-15 06:17:24.019 to 2017-12-15 07:07:24.019) at 0x111cb5198>
seq_scheduler:  Schedule containing 5 observing blocks between 2017-12-14 22:46:00.000 and 2017-12-15 07:07:24.019
"""

"""
2018-01-27a
Constraint calculations with no stored values in scores.sqlite took:
start: '2018-01-27 08:58:34.161705'
end: 2018-01-27 09:35:11.433527'
Around 36 minutes for 12720 entries in the scores table.
"""

"""
2018-01-27b (rerunning same calculations with saved values)
Constraint calculations with all stored values in scores.sqlite took:
start: '2018-01-27 10:07:26.775868'
end: '2018-01-27 10:07:51.919122'
Around 25 seconds for 12720 entries in the scores table.
"""


def main(args):
    """ 
        Entry point for overall execution.
    """
    global mode
    global instrument
    
    # Setting default values here.  These will typically be
    # command line arguments.
    # mode = 2
    # instrument = "binospec"
    setup(args)
    add_observer()
    add_global_constraints()
    add_block_constraints()    
    create_transitioner()

    case = 1
    if case == 1:
        # Normal sequential scheduling
        create_scheduler()
        run_scheduler()
    # Try running a PriorityScheduler
    elif case == 2:
        # Experimental priority scheduling.
        # Don't use for now. 2018-02-26
        print("Running experimental PriorityScheduler()")
        create_pscheduler()
        run_pscheduler()
    elif case == 3:
        # Create score array from schedule.
        print("Create a score array for the schedule...")
        create_scheduler()
        # time_resolution ins in seconds.
        # We add the u.seconds within show_scores()
        # 3600 seconds = 1 hour
        # Times are round to this resolution, e.g., every hour.
        show_scores()
        pass

    return

if __name__ == "__main__":
    redis_client = redis.StrictRedis( host=redis_host, decode_responses=True ) # without that, strings have a 'b' in front
    
    parser = argparse.ArgumentParser(description='mmtscheduler.py')
    parser.add_argument('-m','--mode', help='Scheduling mode.  mode 1: "scheduler" mode; mode 2: "dispatcher" mode.', type=int, choices=[1,2], required=False)
    parser.add_argument('-d','--debug', help='Print debugging messages', required=False)
    parser.add_argument('-v','--version', help='Software version.  Release 1.0 is on ops; version 2.0 and above are on scheduler', type=float, required=False)
    parser.add_argument('-s','--start_date', help='Start time for scheduling (UTC time), e.g., "2017-10-01T03:30:00".  Defines the time used for dispatcher mode and the start time for scheduler mode', required=False)
    parser.add_argument('-e','--end_date', help='End time for scheduling (UTC time), e.g., e.g., "2017-10-01T12:30:00".  Defines the end time for scheduler mode.', required=False)
    parser.add_argument('-i','--schedule_id', help='Queue scheduling run id.', type=int, required=False)
    # New 2017-10-16.  Adding binospec option for instrument.  Defaults to type string.
    parser.add_argument('-t','--instrument', help='Instrument  "mmirs" or "binospec"', required=False)

    args = vars(parser.parse_args())
    if False:
        print("args: ", repr(args))

    main(args)


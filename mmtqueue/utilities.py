import datetime
import json
import math
import requests
import time
from pytz import timezone
import traceback

import sqlite3
# from sqlite3 import Error

__all__ = ['get_score', 'set_score', 'set_scorer', 'set_schedule', 'get_constraint_name', 'dict_factory', \
           'is_numeric', 'localize_time', 'get_config_json', 'get_name',\
           'get_key', 'get_short_key', 'get_now', 'get_target_key', 'get_target_name', \
           'to_allocated_time', 'to_completed_time', 'to_constraint_details', \
           'to_exception', 'to_mysql', 'to_output', 'to_queue', 'to_redis', \
           'to_status', 'to_stdout', 'roundup', 'roundTime', 'roundHour', 'add_prefix', 'get_mask_id', \
           'get_block_id', 'get_program_id', 'table_to_python', 'table_to_entries',
           'get_target_date_key','get_sun_moon_date_key']

debug = True

# Lambda funcstions for reuse.
current_milliseconds = lambda: int(round(time.time() * 1000))
""" Utility lambda function to return the current timestamp in milliseconds
"""

current_datetime = lambda: str(timezone('America/Phoenix').localize(datetime.datetime.now()))



def add_prefix(instrument, key):
    if instrument.lower() == 'binospec':
        prefix = "BINOSPEC."
    else:
        prefix = "MMIRS."
    if False:
        print("from add_prefix: " , prefix + key)
    return prefix + key


# Used by sqlite to return the result as a dictionary.
def dict_factory(cursor, row):
    d = {}
    for idx,col in enumerate(cursor.description):
        d[col[0]] = row[idx]
    return d


# 
# Sample key name in sqlite scores database table.
# Eisenstein-rset56-mask4_305-2*305*8161*459*P1.AirmassConstraint.2018-02-07T02:01:49.559.2018-02-07T03:31:49.559
#
    # Order in key name from member of a "field" structure:
    # 1) objid
    # 2) mask_is
    # 3) block_id
    # 4) program
    # 5) pi_priority
    #
def get_name(f):
    try:
        # We don't want any "*" in the objid or program_id.
        # There shouldn't be, but just to be sure.
        # They would mess up the parsing of the name.
        objid = f['objid'].replace("*","_")
        objid = objid.replace(" ","")
        program_id = f['program_id'].replace("*","_")
        program_id = program_id.replace(" ","")
        # This version has the mask_id in the key name.
        if 'mask_id' in f:
            my_name = objid + '*' + f['mask_id'] + "*" + f['block_id'] + "*" + program_id +  "*P" + str(f['pi_priority'])
        # While this version does not have the mask_id (for older queues)
        else:
            # my_name = f['objid'] + "*" + f['block_id'] + "*" + f['program'] +  "*P" + str(f['pi_priority'])
            my_name = objid + "*" + f['block_id'] + "*" + program_id +  "*P" + str(f['pi_priority'])
    except:
        # Really shouldn't get here...
        my_name = ""

    if False:
        print("my_name: ", my_name)
    return my_name

# Extracts the block_id out of the target name

def get_block_id(name):
    # Assumes that there are no "*" in the user's target name.
    arr = name.split('*')
    # Case of having a start_date and end_date
    if len(arr) == 5:
        block_id = arr[2]
    # Case of having only a start_date (no longer used)
    elif len(arr) == 4:
        block_id = arr[1]
    else:
        print("Key name does not include block_id.  Length should be 4 or 5")
        block_id = 0
        
    if False:
        print("block_id: ", block_id)
    return block_id

# Extracts the mask_id out of the target name
def get_mask_id(name):
    arr = name.split('*')
    # Need to have two versions: one that includes the mask_id 
    # within the name and the other that doesn't (and fails)
    if len(arr) == 5:
        mask_id = arr[1]
    else:
        print("Key name does not include mask_id.  Length should be 5")
        mask_id = 0
        
    if False:
        print("mask_id: ", mask_id)
    return mask_id


# Extracts the program out of the target name
def get_program_id(name):
    arr = name.split('*')
    program = arr[3] 
    
    if len(arr) == 5:
        program = arr[3]
    elif len(arr) == 4:
        program = arr[2]
    else:
        print("Key name does not include program.  Length should be 4 or 5")
        program = 0
        
    if False:
        print("program: ", program)
    return program


# For a target, returns a key with the time rounded to a integer modified
# Julian day.
# mode can be 'set', 'rise' or 'transit' 
def get_target_date_key(target, time, mode='set'):
    """
    
    """
    key = ''
    try:
        
        if target.name == 'icrs':
            # It's a astropy SkyCoord.  Get the "name" as the coordinates.
            target_name = target.to_string('hmsdms')
        else:
            target_name = target.name
        # Getting the integer value of the modified Julian day.
        time.format = 'mjd'
        key = str(target_name) + '.' + mode + '.' + str(int(time.mjd))
    except:
        print("traceback: ", traceback.format_exc())

    if False:
        print("get_target_date_key: ", key)
    return key

# mode could be: sunrise, sunset, moonrise, moonset, sun_transit, moon_transit
def get_sun_moon_date_key(mode, time):
    """
    
    """
    key = ''
    try:
        # Getting the integer value of the modified Julian day.
        time.format = 'mjd'
        key = mode + "." + str(int(time.mjd))
    except:
        print("traceback: ", traceback.format_exc())

    if False:
        print("get_sun_moon_date_key: ", key)
    return key


# This one has no constraint
def get_short_key(target, time, time2=None):
    """
    
    """
    key = ''
    try:
        time.format = 'isot'
        # Typically going with both a start and end time.
        # So, there would be both "time" and "time2" parameters
        if target.name == 'icrs':
            # It's a astropy SkyCoord.  Get the "name" as the coordinates.
            target_name = target.to_string('hmsdms')
        else:
            target_name = target.name

        # Remove any whitespace in target name.
        target_name = target_name.replace(" ", "")
        if time2 == None:
            key = "{}.{}".format(target_name,
                                        str(time))
        else:
            time2.format = 'isot'
            key = "{}.{}.{}".format(target_name,
                            str(time),
                            str(time2))
    except:
        print("traceback: ", traceback.format_exc())

    if False:
        print("key:", key)
    return key

# This one has a constraint.
def get_key(constraint, target, time, time2=None):
    """
    
    """
    key = ''
    try:
        constraint_name = get_constraint_name(constraint)
        time.format = 'isot'
        # Typically going with both a start and end time.
        # So, there would be both "time" and "time2" parameters
        if target.name == 'icrs':
            # It's a astropy SkyCoord.  Get the "name" as the coordinates.
            target_name = target.to_string('hmsdms')
        else:
            target_name = target.name

        # Remove any whitespace in target name.
        target_name = target_name.replace(" ", "")
        if time2 == None:
            key = "{}.{}.{}".format(target_name,
                                        constraint_name,
                                        str(time))
        else:
            time2.format = 'isot'
            key = "{}.{}.{}.{}".format(target_name,
                            constraint_name,
                            str(time),
                            str(time2))

    except:
        print("traceback: ", traceback.format_exc())

    if False:
        print("key:", key)
    return key


def get_score(conn, key):
    """Gets a score value from the mmtqueue.sqlite:scores database table

    Keyword arguments:
    :param conn: sqlite3 Connection
    :param key: 'key' for lookup in the schedules database
    :type conn: sqlite3.connection
    :type key: String
    """
    if False:
        print("Entering get_score: ", repr(conn), ", ", repr(key))
    row = None
    if conn is not None:
        try:
            c = conn.cursor()
            c.execute("SELECT * FROM scores WHERE key='{}'".format(key))
            row = c.fetchone()
        except sqlite3.DatabaseError:
            print('DatabaseError: Cannot execute select for {}'.format(key))
        except sqlite3.Error:
            print('ERROR: Cannot execute select for {}'.format(key))
    return row


def set_schedule(conn, key, value):
    """Saves a schedule output to the mmtqueue.sqlite:schedules database table

    Keyword arguments:
    :param conn: sqlite3 Connection
    :param key: 'key' for lookup in the schedules database
    :param value: 'value' to insert/replace in the schedules database
    :type conn: sqlite3.connection
    :type key: String
    :type value: String
    """
    # Using replace rather then insert so that fields can be updated.
    if conn is not None:
        try:
            sql = "REPLACE INTO schedules (key, value) VALUES ('{}', '{}')".format(key, value)
            if True:
                print("sql: ", sql)
            
            # Can do this with new versions of sqlite.  It does commits automatically.
            with conn:
                conn.execute(sql)
        except sqlite3.IntegrityError:
            print('setschedule IntegrityError: ID already exists in PRIMARY KEY column {}'.format("key"))
        except sqlite3.DatabaseError:
            _print('setschedule DatabaseError: Cannot execute replace for {}'.format(key))
        except sqlite3.Error:
            print('set_scheduleERROR: Cannot execute replacd for {}'.format(key))


def set_scorer(conn, key, value, time='', schedule_id=0, comment=''):
    """Sets a score in the mmtqueue.sqlite.scorer database table.



    """

    if conn is not None:
        try:
            sql = "REPLACE INTO scorer VALUES ('{}', {}, '{}', {}, '{}')".format(key, value, time, schedule_id, comment)
            if False:
                print("sql: ", sql)
            # Can do this with new versions of sqlite.  It does commits automatically.
            with conn:
                conn.execute(sql)
    
        except sqlite3.IntegrityError:
            print('set_score IntegrityError: ID already exists in PRIMARY KEY column {}'.format(key))
        except sqlite3.DatabaseError:
            print('set_score DatabaseError: Cannot execute replace for {}'.format(key))
        except sqlite3.Error:
            print('set_score ERROR: Cannot execute replace for {}'.format("key"))


def set_score(conn, key, value, time='', schedule_id=0, comment=''):
    """Sets a score in the mmtqueue.sqlite.scores database table.



    """

    if conn is not None:
        try:
            sql = "REPLACE INTO scores VALUES ('{}', {}, '{}', {}, '{}')".format(key, value, time, schedule_id, comment)
            if False:
                print("sql: ", sql)
            # Can do this with new versions of sqlite.  It does commits automatically.
            with conn:
                conn.execute(sql)
    
        except sqlite3.IntegrityError:
            print('set_score IntegrityError: ID already exists in PRIMARY KEY column {}'.format(key))
        except sqlite3.DatabaseError:
            print('set_score DatabaseError: Cannot execute replace for {}'.format(key))
        except sqlite3.Error:
            print('set_score ERROR: Cannot execute replace for {}'.format("key"))


def get_config_json(schedule_id: object) -> object:
    url = 'https://scheduler.mmto.arizona.edu/QueueSchedules/config_json.php?formatted=0&schedule_id='
    url += str(schedule_id)
    retries = 0
    txt = ''
    # Making this request more robust by trying several times to get the git
    # configuration json string.
    # The text should easily be over 100 characters, probably > 10000.
    while retries < 5 and len(txt) < 100:
        r = requests.get(url)
        txt = r.text
        retries += 1
    txt = txt.replace("<pre>\n","")
    txt = txt.replace("</pre>\n","")
    txt = txt.replace("\s+","")
    txt = txt.replace("\n","")
    txt = txt.replace("&amp;&amp;","")
    if False:
        print(txt)
    obj = json.loads(txt)
    return obj

# Returns the name of a constraint.
def get_constraint_name(constraint):
    if constraint is None:
        name = ''
    else:
        name = type(constraint).__name__
    if False:
        print("Constraint name: ", name)
    return name





def get_now():
    return str(datetime.datetime.now())



# Returns teh target name from an observing block
def get_target_name(block):
    return block.target.name


def get_target_key(target_name, time, time2=None):
    """
    
    """
    constraint_name = "OverallScore"
    time.format = 'isot'
    # Typically going with both a start and end time.
    # So, there would be both "time" and "time2" parameters
    if time2 == None:
        key = "{}.{}.{}".format(target_name, 
                        constraint_name,
                        str(time))
    else:
        time2.format = 'isot'
        key = "{}.{}.{}.{}".format(target_name, 
                        constraint_name,
                        str(time),
                        str(time2))
        
    if False:
        print("key:", key)
    return key

def is_numeric(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def localize_time(t):
    # return str(timezone('America/Phoenix').localize(t.datetime - datetime.timedelta(hours=7)))
    return timezone('America/Phoenix').localize(t.datetime - datetime.timedelta(hours=7)).strftime('%Y-%m-%d %H:%M:%S')
    
def to_queue(txt, mode, redis_client, instrument, do_redis):
    if debug:
        print(txt)

    if do_redis:
        if mode == 2:
            key = 'DISPATCHER.queue'
        else:
            key = 'SCHEDULER.queue'
        
        # Old version without instrument prefix.
        if instrument == 'mmirs':
            redis_client.set(key, txt)
            redis_client.publish(key, txt)

        # New version with instrument prefix
        redis_client.set(add_prefix(instrument, key), txt)
        redis_client.publish(add_prefix(instrument, key), txt)
        
def to_redis(key, value, redis_client, do_redis):
    if debug:
        print(value)
    # do_redis is a global Boolean on whether we
    # want to post to the Redis server.
    # It is usually True, but if you don't
    # want to interfere with other calculations,
    # you can turn if off.
    try:
        if do_redis:
            d_json = to_json(key, value)
            """
            #
            # This has moved to to_json()
            #
            d = {'key':key,
                 'value': value,
                 'datetime': current_datetime(),
                 'timestamp': current_milliseconds()}
            d_json = json.dumps(d,sort_keys=True,ensure_ascii=False)
            d_json = d_json.replace('\\\"','\"')
            # Not sure that we need this.
            if True:
                d_json = d_json.replace('\"','"')
            
            # Removing the double quotes from an array.
            if True:
                if '"[' in d_json:
                    print("d_json1: ", d_json)
                    d_json = d_json.replace('"[','[')
                    print("d_json2: ", d_json)
                if ']"' in d_json:
                    print("d_json3: ", d_json)
                    d_json = d_json.replace(']"',']')
                    print("d_json4: ", d_json)
            else:
                d_json = d_json.replace('"[','[')
                d_json = d_json.replace(']"',']')
            """ 
            redis_client.set(key, d_json)
            redis_client.publish(key, d_json) 
    except Exception as e:
        print("to_redis JSON error")
        print("traceback: " , traceback.format_exc())
        
# This is a helper function to log the outputs from the dispatcher and scheduler
# to a mysql table.  This is done by calling a URL with the correct schedule_id.
# Different CGI par are used for the dispatcher and the scheduler.
def to_mysql(sid, mode, redis_client, instrument, do_redis):
    txt = ""
    try:
        # Using scheduler for everything now:
        # url = 'https://scheduler.mmto.arizona.edu/QueueSchedules/notify.php?message=complete&type=SCHEDULER&schedule_id=' + str(sid)
        if mode == 2:
            url = 'https://scheduler.mmto.arizona.edu/QueueSchedules/notify.php?message=complete&type=DISPATCHER&schedule_id=' + str(sid)
        else:
            url = 'https://scheduler.mmto.arizona.edu/QueueSchedules/notify.php?message=complete&type=SCHEDULER&schedule_id=' + str(sid)
      
        """
        if mode == 2:
            url = 'https://ops.mmto.arizona.edu/ObservatoryManager/QueueSchedules/notify.php?message=complete&type=DISPATCHER&schedule_id=' + str(sid)
        else:
            url = 'https://ops.mmto.arizona.edu/ObservatoryManager/QueueSchedules/notify.php?message=complete&type=SCHEDULER&schedule_id=' + str(sid)
        """
        txt = "url = " + url
        # Don't need the response for now.
        requests.get(url)
    except:
        txt = "Error logging to MySQL via URL."
    to_stdout("Logging to mysql via URL:" + txt, mode, \
              redis_client, instrument, do_redis)

def to_stdout(txt, mode, redis_client, instrument, do_redis):
    # Old version without instrument prefix.
    if mode == 2:
        key = 'DISPATCHER.stdout'
    else:
        key = 'SCHEDULER.stdout'

    # For now, pushing output from both instrument to the same Redis parameter.
    if True or instrument.lower() == 'mmirs':
        to_redis(key, txt, redis_client, do_redis)
    # New version with instrument prefix.
    to_redis(add_prefix(instrument, key),txt, redis_client, do_redis)

def to_status(txt, mode, redis_client, instrument, do_redis):
    # Old version without instrument prefix.
    if mode == 2:
        key = 'DISPATCHER.status'
    else:
        key = 'SCHEDULER.status'
    
    if True or instrument.lower() == 'mmirs':
        to_redis(key, txt, redis_client, do_redis)
    # New version with instrument prefix.
    to_redis(add_prefix(instrument, key), txt, redis_client, do_redis)

def to_json(key, value):
    d_json = ''
    if False:
        print("In to_json: key: ", key, ", value: ", repr(value))
    try:
        j =json.dumps(value,sort_keys=True,ensure_ascii=False)
        data = value
    except Exception as e:
        print("traceback1: " , traceback.format_exc())   

        try:
            d = value[0]
            j =json.dumps(d,sort_keys=True,ensure_ascii=False)
            data = d
        except Exception as e:
            print("traceback2: " , traceback.format_exc())   
    
    try:
        d = {'key':key,
                 'value': data,
                 'datetime': current_datetime(),
                 'timestamp': current_milliseconds()}
        d_json = json.dumps(d,sort_keys=True,ensure_ascii=False)
        d_json = d_json.replace('\\\"','\"')
        # Not sure that we need this.
        if True:
            # This is a "raw" string replacement
            d_json = d_json.replace(r'\"','"')
            # This would also work.
            d_json = d_json.replace('\\"','"')
            # This doesn't do anything???
            d_json = d_json.replace('\"','"')
        
        # Removing the double quotes from an array.
        if False:
            if '"[' in d_json:
                print("d_json1: ", d_json)
                d_json = d_json.replace('"[','[')
                print("d_json2: ", d_json)
            if ']"' in d_json:
                print("d_json3: ", d_json)
                d_json = d_json.replace(']"',']')
                print("d_json4: ", d_json)
        else:
            d_json = d_json.replace('"[','[')
            d_json = d_json.replace(']"',']')
    except Exception as e:
        print("to_redis JSON error")
        print("traceback: " , traceback.format_exc())


    return d_json

def to_output(conn, tbl, mode, redis_client, instrument, schedule_id, do_redis):
    """ Utility function to set/publish to Redis the completed time by program.

    The scheduler details is a large JSON structure that contains all of the
    computational subproducts during scheduling. This includes each constraint score
    for each observing block for each time slot. It is similar to the dispatcher
    output, but for the entire scheduling run.

    """
    if mode == 2:
        if schedule_id is not None:
            key = "DISPATCHER.{}.output".format(schedule_id)
            to_redis(add_prefix(instrument, key), tbl, redis_client, do_redis)
            # Putting into constraints.schedules.
            d_json = to_json(add_prefix(instrument, key), tbl)
            if True:
                print("Calling set_schedule.")
            set_schedule(conn, add_prefix(instrument, key), d_json)
        key = 'DISPATCHER.output'
    else:
        if schedule_id is not None:
            key = "SCHEDULER.{}.output".format(schedule_id)
            to_redis(add_prefix(instrument, key), tbl, redis_client, do_redis)
            # Putting into constraints.schedules.
            d_json = to_json(add_prefix(instrument, key), tbl)
            if True:
                print("Calling set_schedule.")
            set_schedule(conn, add_prefix(instrument, key), d_json)
        key = 'SCHEDULER.output'

    to_redis(add_prefix(instrument,key), tbl, redis_client, do_redis)
    

# =============================================================================
#     # New version with instrument prefix.
#     if instrument.lower() == 'binospec':
#         prefix = "BINOSPEC."
#     else:
#         prefix = "MMIRS."
# 
#     if mode == 2:
#         if schedule_id is not None:
#             key = prefix + "DISPATCHER.{}.output".format(schedule_id)
#             to_redis(key,txt, redis_client)
#         key = prefix + 'DISPATCHER.output'
#     else:
#         if schedule_id is not None:
#             key = prefix + "SCHEDULER.{}.output".format(schedule_id)
#             to_redis(key,txt, redis_client)
#         key = prefix + 'SCHEDULER.output'
#     to_redis(add_prefix(instrument, key), txt, redis_client)
# =============================================================================
        
def to_constraint_details(txt, mode, redis_client, instrument, do_redis):

    """ Utility function to set/publish to Redis the constraint details for an entire schedule.

    The SCHEDULER.details Redis parameter is a large JSON structure that contains all of the
    computational subproducts during scheduling. This includes each constraint score
    for each observing block for each time slot. It is similar to the dispatcher
    output, but for the entire scheduling run.

    """
    # Old version without instrument prefix.
    if mode == 2:
        key = 'SCHEDULER.details'
    else:
        key = 'DISPATCHER.details'

    if True or instrument == 'mmirs':
        to_redis(key,txt,redis_client, do_redis)

    # New version with instrument prefix.
    to_redis(add_prefix(instrument, key),txt,redis_client, do_redis)


def to_allocated_time(txt, mode, redis_client, instrument, do_redis):
    """ Utility function to set/publish to Redis the allocated time by program.

    The DISPATCHER.allocated_time and SCHEDULER.allocated_time Redis parameters are
    JSON structures of the allocated time for each program.  The allocated time is 
    used to compute a priority for each observing block, based upon how much time
    each program has already used.
    

    """
    # Old version without instrument prefix.
    if mode == 2:
        key = 'DISPATCHER.allocated_time'
    else:
        key = 'SCHEDULER.allocated_time'
    
    if True or instrument == 'mmirs':
        to_redis(key, txt, redis_client, do_redis)

    # New version with instrument prefix.
    to_redis(add_prefix(instrument, key),txt, redis_client, do_redis)

def to_completed_time(txt, mode, redis_client, instrument, do_redis):
    """ Utility function to set/publish to Redis the completed time by program.

    The DISPATCHER.allocated_time and SCHEDULER.allocated_time Redis parameters are
    JSON structures of the allocated time for each program.  The allocated time is 
    used to compute a priority for each observing block, based upon how much time
    each program has already used.
    
    """

    # Old version without instrument prefix.
    if mode == 2:
        key = 'DISPATCHER.completed_time'
    else:
        key = 'SCHEDULER.completed_time'
    
    if True or instrument == 'mmirs':
        to_redis(key,txt, redis_client, do_redis)
    # New version with instrument prefix.
    to_redis(add_prefix(instrument, key),txt, redis_client, do_redis)
    
def to_exception(txt, mode, redis_client, instrument, do_redis):
    """ Utility function to set/publish to Redis any run-time exceptions.
    """

    # Old version without instrument prefix.
    if mode == 2:
        key = 'DISPATCHER.exception'
    else:
        key = 'SCHEDULER.exception'

    if True or instrument == 'mmirs':
        to_redis(key,txt,redis_client, do_redis)
    # New version with instrument prefix.
    to_redis(add_prefix(instrument, key),txt, redis_client, do_redis)



# General function to round a number up to a multiple of "interval".
def roundup(x, interval):
    """ Utility function to round a float number up to the next higher integer.
    """
    return int(math.ceil(float(x) / float(interval)) * interval)

#From:  https://stackoverflow.com/questions/3463930/how-to-round-the-minute-of-a-datetime-object-python/10854034#10854034
def roundTime(dt=None, roundTo=60):
   """Round a datetime object to any time laps in seconds
   dt : datetime.datetime object, default now.
   roundTo : Closest number of seconds to round to, default 1 minute.
   Author: Thierry Husson 2012 - Use it as you want but don't blame me.
   """
   if dt == None : dt = datetime.datetime.now()
   seconds = (dt.replace(tzinfo=None) - dt.min).seconds
   rounding = (seconds+roundTo/2) // roundTo * roundTo
   return dt + datetime.timedelta(0,int(rounding-seconds),-dt.microsecond)

   """
   # Rounding down to the minute, e.g., 10 minutes.
   # From https://stackoverflow.com/questions/3463930/how-to-round-the-minute-of-a-datetime-object-python
   """
   # dt = dt - datetime.timedelta(minutes=dt.minute % roundTo,
   #                           seconds=dt.second,
   #                           microseconds=dt.microsecond)
   # return dt

def roundHour(dt=None):
   """
   Rounding down to the hour.
   From https://stackoverflow.com/questions/3463930/how-to-round-the-minute-of-a-datetime-object-python
   """
   if dt == None : dt = datetime.datetime.now()
   dt = dt - datetime.timedelta(minutes=dt.minute,
                             seconds=dt.second,
                             microseconds=dt.microsecond)
   return dt

# convert the internal astroplan table structure to a Python structure (usually a dictioary.)
# From:  https://github.com/astropy/astropy/issues/4604
def table_to_python(table):
    """Convert Astropy Table to Python dict.

    Numpy arrays are converted to lists, so that
    the output is JSON serialisable.

    Can work with multi-dimensional array columns,
    by representing them as list of list.
    """
    total_data = {}
    for name in table.colnames:
        data = table[name].tolist()
        total_data[name] = data
    return total_data

# convert the internal astroplan table structure to a Python structure that can
# be used as "entries" for the javascript fullcalendar.
# From:  https://github.com/astropy/astropy/issues/4604
def table_to_entries(table, block_ids, objid_ids):
    """
    Convert Astropy Table to Python list, suitable for JSON.
    Resulting fields are inputs to fullcalendar.io entries as JSON.
    """
    entries = []
    for row in table:
        # Not including TransitionBlocks
        if row['target'] != "TransitionBlock":
            name = row['target']
            
            block_id = ""
            objid_id = ""
            if name in block_ids:
                block_id = block_ids[name]
            if name in objid_ids:
                objid_id = objid_ids[name]
            
            data =  {"title": name,  \
                    "block_id": block_id, \
                    "objid_id": objid_id, \
                    "start": row['start time (UTC)'] + " UTC", \
                    "end": row['end time (UTC)'] + " UTC"}
            entries.append(data)
    return entries

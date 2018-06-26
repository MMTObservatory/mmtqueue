#!/usr/bin/python3.6
#
"""
This program listens on a Redis parameter for new commands to recalculated the dispatcher or scheduler
queues.  The same Redis parameters handles both dispatcher and scheduler commands as well as both
the MMIRS and BinoSpec instruments.  This program calls another program, mmtsheduler.py, that does the
actual calcuations of a new dispatcher output or queue schedule.  

Much of the complexity of this program is keeping track of the type of calculation (dispatcher or scheduler)
and the instrument (MMIRS or BinoSpec).  Additional complexity is introduced because of the migration
of the related web-based API's from ops to the vSphere-based dbshare/scheduler virtual machine.

There are two instances of this program running on the cruncher computer:  dispatcher.service and 
scheduler.service.  When a new command is posted to Redis, any existing calculation is cancelled and
a new one is started.  Also of importance is that the same process computes dispatcher or scheduler
outputs for both the MMIRS and BinoSpec instruments.  Only one dispatcher or scheduler calculation
can be run for either MMIRS or BinoSpec.

A "version" parameter is being used for the transition from "ops" to "scheduler" as well as the
transition from astroplan 0.2.1 to more recent version, currently 0.4.  Astropy and Numpy have
also moved to newer versions.

"""
import argparse
import redis
import time
import sys
import subprocess
import os
import datetime
import signal
import json
from pytz import timezone

# For reuse.
current_milliseconds = lambda: int(round(time.time() * 1000))
current_datetime = lambda: str(timezone('America/Phoenix').localize(datetime.datetime.now()))
# current_datetime = lambda: timezone('America/Phoenix').localize(datetime.datetime.now()).strftime("%Y-%m-%d %H:%M:%S")

# Global parameters
debug = True 
# Redis clients.
redis_client = None  # get/set
redis_pubsub = None  # pub/sub
thread = None
# Defaults to mode==2 (dispatcher)
mode = 2
# PID structures of either the scheduler or dispatcher subprocess.
pid = None

cwd = None

# Redis callback when receiving a change in the DISPATCHER.command or SCHEDULER.command redis parameter
# See a sample JSON structure for this parameter above.
# msg['command'] values can be "stort" or "stop"
def on_command_message( message ):
    global mode, pid, debug, version, instrument, cwd
    if debug:
        print ("message", message)
    try:
        # Convert the Redis data string to JSON and extract the "value" parameter from that structure.
        m = json.loads(message['data'])
        msg = m['value']
        if debug:
            print ("msg", msg)
       
        # Default to version 1.0.  This assumes that the API will be queried on ops.mmto.arizona.edu.
        version = 2.0

        # The "version" parameter of the message can set a new version number, typically 2.0.
        # Version 2.0 uses the API on scheduler.mmto.arizone.edu.
        # 
        # Note that if no version is given, the API uses ops.mmto.arizona.edu
        # This makes the code backwards-compatible.
        if 'version' in msg:
            try: 
                float(msg['version'])
                version = msg['version']
            except:
                print('Invalid version number')

        # Received a "start" command.  This can be for either the dispatcher or scheduler.
        if msg['command'] == 'start':
            # Opening a file handle to devnull.
            # FNULL = open(os.devnull, 'w')
            # Trying to open a log file.
            FNULL = open('log.txt', 'w')
            
            if debug:
                txt = "%s, Received start request: %s\n" % (time.strftime('%Y-%m-%d %H:%M:%S'), message)
                to_stdout(txt)
                
            # Make sure that we have everything that we need in the JSON structure before
            # attempting to do any calculations.
            if 'scheduler_id' in msg and  \
                int(msg['scheduler_id']) > 0 and \
                'start_date' in msg and \
                'end_date' in msg:
                
                # Note: Make sure that mmtscheduler is in the path and using the correct version of python.
                # 
                # I'm currently using a symbolic link from /usr/local/bin/mmtscheder.py to this file
                # in the mmtqueue subdirectory on cruncher.
                if instrument == 'binospec':
                    version = 2.0

                # This is a hack to force binospec to use the scheduler.mmto.arizona.edu API.
                if str(version) == '2.0':
                    instrument = 'binospec'

                # cmd = "cd /home/mmtqueue/src/QueueScheduler/;./mmtscheduler.py -m=" + str(mode) \
                # + " -t=" + instrument \
                cmd = "./mmtscheduler.py -m=" + str(mode) \
                        + " -i=" + str(msg['scheduler_id']) \
                        + " -s=" + msg['start_date'] \
                        + " -e=" + msg['end_date'] \
                        + " -v=" + str(version)


                
                if debug:
                    txt = "%s, cmd: %s" % (datetime.datetime.now(), cmd)
                    print("cmd: ", txt)
                    to_stdout(txt)
                
                # Mode 1 is the scheduler
                # Mode 2 is the dispatcher
                mode = int(mode)
                if mode == 1 or mode == 2:
                    # Kill the subprocess if it exists.
                    if pid is not None:
                        stop_computing(mode)
                            
                    # Create a new subprocess
                    #
                    # Start a process and keep track of the pid.
                    # The os.setsid() is passed in the argument preexec_fn so
                    # it's run after the fork() and before exec() to run the shell.
                    #
                    # Redirecting subprocess output to DEVNULL
                    if True:
                        print("Before pid: ", repr(pid))
                    # pid = subprocess.Popen(cmd, stdout=FNULL, stderr=subprocess.STDOUT, 
                    #            stdin=subprocess.PIPE, shell=True, preexec_fn=os.setsid)
                    #
                    # pid = subprocess.Popen(cmd, stdout=FNULL, stderr=subprocess.STDOUT, 
                    #       stdin=subprocess.PIPE, shell=True, cwd=cwd,
                    #       preexec_fn=os.setsid)
                    #
                    pid = subprocess.Popen(cmd, shell=True, cwd=cwd)
                    if True:
                        print("After pid: ", repr(pid))
                    
                    # Compose the JSON for the DISPLATCHER.running/SCHEDULER.running redis parameter.
                    stuff = {'scheduler_id': msg['scheduler_id'],
                             'start_date': msg['start_date'],
                             'end_date': msg['end_date'],
                             'pid': pid.pid,
                             'timestamp': current_milliseconds(),
                             'datetime': current_datetime()}
                    
                    # Post to the DISPATCHER/SCHEDULER.running Redis parameter
                    # This is shared between Binospec and MMIRS.
                    to_running(stuff)
                                                          
                    # Post an output line to DISPATCHER/SCHEDULER.stdout Redis parameter.
                    if mode == 1:
                        txt = "Starting a scheduler at %s with pid %s" % (time.strftime('%Y-%m-%d %H:%M:%S'), pid.pid)
                    else:
                        txt = "Starting a dispatcher at %s with pid %s" % (time.strftime('%Y-%m-%d %H:%M:%S'), pid.pid)
                    print(txt)
                    to_stdout(txt)
                  
                else:
                    txt = "Invalid mode"
                    print (txt)
                    to_stdout(txt)
          
        # Received a "stop" Redis command.  The "mode" will sort out which PID to kill.
        elif msg['command'] == 'stop':
            if debug:
                txt = "%s, Received stop request: %s\n" % (time.strftime('%Y-%m-%d %H:%M:%S'), message)
                to_stdout(txt)
            
            if mode == 1 or mode == 2:
                to_status('stopped')
                stop_computing(mode)
                to_stdout("Stopped computation")
            else:
                to_stdout("Invalid stop mode: only 1 or 2 are allowed.")
        
        # Shouldn't get here.
        else:
            pass
    except Exception as e: 
       to_stdout("on_command_message error, exception: %s" % repr(e))
    
    return         

# This function kills a running dispatcher or scheduler process, if one exists.  
# We don't want multiple instances of dispatchers and schedulers running, only
# one of each.
def stop_computing(m):
    global pid, debug
    try:
        if m == 1 or m == 2:
            if pid is not None:
                try:
                    if mode == 1:
                        txt = "Sending a scheduler kill command on pid %s" % (pid.pid)
                    else:
                        txt = "Sending a dispatcher kill command on pid %s" % (pid.pid)
                    to_stdout(txt)
                    # os.killpg(os.getpgid(pid.pid), signal.SIGTERM)  # Send the signal to all the process groups
                    pid.kill()
                    pid = None
                    
                    # Compose the JSON for the DISPATCHER.running or SCHEDULER.running redis parameter.
                    stuff = {'scheduler_id': 'NA',
                             'start_date': 'NA',
                             'end_date': 'NA',
                             'pid': 'NA',
                             'timestamp': current_milliseconds(),
                             'datetime': current_datetime()}
                    to_running(stuff)
                    
                    # Change DISPATCHER.status or SCHEDULER.status to "stopped"
                    to_status('stopped')
                except Exception as e: 
                    to_stdout("Process probably done, exception: %s" % repr(e))
            else:
                to_stdout("Process is not running.")
        
        elif m == "NA":
            pass
        
        else:
            to_stdout("Invalid stop data: only 1, 2 or NA are allowed.")
    except Exception as e: 
        to_stdout("stop_computing, exception: %s" % repr(e))

def to_redis(key,value):
    global mode, redis_client, debug
    if debug:
        print(value)
    d = {'key':key,
         'value': value,
         'datetime': current_datetime(),
         'timestamp': current_milliseconds()}
    d_json = json.dumps(d,sort_keys=True)
    d_json = d_json.replace('\\\"','\"')
    redis_client.set(key, d_json)
    redis_client.publish(key, d_json) 

def to_running(txt):
    if mode == 2:
        key = 'DISPATCHER.running'
    else:
        key = 'SCHEDULER.running'
    to_redis(key,txt)
        
def to_status(txt):
    if mode == 2:
        key = 'DISPATCHER.status'
    else:
        key = 'SCHEDULER.status'
    to_redis(key,txt)

def to_stdout(txt):
    if mode == 2:
        key = 'DISPATCHER.stdout'
    else:
        key = 'SCHEDULER.stdout'
    to_redis(key,txt)

def to_exception(txt):
    global mode
    if mode == 2:
        key = 'DISPATCHER.exception'
    else:
        key = 'SCHEDULER.exception'
    to_redis(key,txt)


if __name__ == "__main__":
    # Initialize Redis clients.
    # Redis client "redis_client" is for get/set commands.
    redis_client = redis.StrictRedis( host='redis.mmto.arizona.edu', decode_responses=True ) # without that, strings have a 'b' in front
    
    # Create a separate Redis client, "redis_client2", for publish/subscribe commmands.
    redis_client2 = redis.StrictRedis( host='redis.mmto.arizona.edu', decode_responses=True ) # without that, strings have a 'b' in front
    redis_pubsub = redis_client2.pubsub( ignore_subscribe_messages=True ) 
    
    # Running the Redis publish/subscribe cline in a separate thread.
    thread = redis_pubsub.run_in_thread( sleep_time=0.001 )
    
    parser = argparse.ArgumentParser(description='mmtscheduler.py')
    
    # Mode defaults to 2 (dispatcher)
    parser.add_argument('-m','--mode', help='Scheduling mode.  mode 1: "scheduler" mode; mode 2: "dispatcher" mode.', type=int, choices=[1,2], required=False)
  
    # New 2017-10-16.  Adding binospec option for instrument.  Defaults to type string.
    parser.add_argument('-t','--instrument', help='Instrument  "mmirs" or "binospec"', required=False)
   
 
    # New 2018-04-05.  Adding an optional current working directory,.
    parser.add_argument('-c','--cwd', help='Current working direction', required=False)

    # Parse the command line argument.
    args = vars(parser.parse_args())
    
    if not args['mode'] is None:
        mode  = int(args['mode'])
        txt = "Setting mode to: {0}".format(mode)
        to_stdout(txt)
    else:
        txt = "Default mode is 2 (dispatcher)"
        to_stdout(txt)
   
    instrument = 'mmirs'
    if not args['instrument'] is None:
        if args['instrument'] == 'binospec':
            instrument  = args['instrument']
            print("Setting instrument to: ", instrument)

    # cwd = "/home/mmtqueue/src/QueueScheduler/   "
    # cwd = "/home/mmtqueue/git/QueueScheduler2.0/"
    cwd = "/Users/jdgibson/git/QueueScheduler2.0/"
    if not args['cwd'] is None:
        cwd  = args['cwd']
        print("Setting current working directory to: ", cwd)
    else:
        print("Current working directory is: ", cwd)

    # Use Redis to subscribe to the appropriate Redis command:
    #   Mode 1 is the scheduler,
    #   Mode 2 is the dispatcher.
    if mode == 1:
        redis_pubsub.subscribe( **{ 'SCHEDULER.command': on_command_message } )
        redis_pubsub.subscribe( **{ 'BINOSPEC.SCHEDULER.command': on_command_message } )
        redis_pubsub.subscribe( **{ 'MMIRS.SCHEDULER.command': on_command_message } )
        # if instrument == 'binospec':
        #     redis_pubsub.subscribe( **{ 'BINOSPEC.SCHEDULER.command': on_command_message } )
        # else:
        #     redis_pubsub.subscribe( **{ 'MMIRS.SCHEDULER.command': on_command_message } )

    else:
        redis_pubsub.subscribe( **{ 'DISPATCHER.command': on_command_message } )
        redis_pubsub.subscribe( **{ 'BINOSPEC.DISPATCHER.command': on_command_message } )
        redis_pubsub.subscribe( **{ 'MMIRS.DISPATCHER.command': on_command_message } )
        # if instrument == 'binospec':
        #     redis_pubsub.subscribe( **{ 'BINOSPEC.DISPATCHER.command': on_command_message } )
        # else:
        #    redis_pubsub.subscribe( **{ 'MMIRS.DISPATCHER.command': on_command_message } )

    
    # os.setpgrp() # create new process group, become its leader
    try:
        while True:
            for m in redis_pubsub.listen():
                if m['type'] == 'subscribe':
                    if m['data'] == 1:
                        to_stdout('subscribed to: %s' % (m['channel']))
                # This is what most pubsub submissions will match
                elif m['type'] == 'message':
                    # Our payload is in the data field of the response
                    to_stdout(json.loads(m['data']))
                else:
                    to_stdout(m)
    except KeyboardInterrupt:
        try:
            # Check if subprocesses need to be terminated.
            if pid is not None:
                    pid.terminate()
        except OSError:
            pass
    except:
        to_exception("Unexpected error: " + sys.exc_info()[0])
    finally:
        to_status('stopped')
        os.killpg(0, signal.SIGTERM) # kill all processes in my group
        os.killpg(0, signal.SIGKILL) # kill all processes in my group

# Sample systemd services.
# 
# queue@queue:/lib/systemd/system$ pwd
# /lib/systemd/system
# queue@queue:/lib/systemd/system$ cat mmtdispatcher.service
# [Unit]
# Description=Server control for MMIRS/Binospec queue dispatcher
# After=syslog.target network.target
# 
# [Service]
# Type=simple
# Restart=always
# ExecStart=/usr/bin/python /usr/local/bin/mmtscheduler_server.py -m=2
# 
# [Install]
# WantedBy=multi-user.target
# queue@queue:/lib/systemd/system$ cat mmtscheduler.service
# [Unit]
# Description=Server control for MMIRS/Binospec queue scheduler
# After=syslog.target network.target
# 
# [Service]
# Type=simple
# Restart=always
# ExecStart=/usr/bin/python /usr/local/bin/mmtscheduler_server.py -m=1
# 
# [Install]
# WantedBy=multi-user.target
# queue@queue:/lib/systemd/system$

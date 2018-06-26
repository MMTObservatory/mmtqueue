#!/usr/bin/env python
#
# This script uses Redis to control the mmtscheduler.py script for the MMIRS/Binospec queue scheduler and dispatcher.
# The scheduler and dispatcher are run as separate processes and started (or restarted) by subscripting to Redis parameters.
# If a new request to run a scheduler or dispatcher is received, any currently running
# subprocess is killed and a new subprocess started.
# 
# Redis parameters:
#    SCHEDULER.command  (JSON structure of scheduler_id, mode, start_date, end_date, and command.  Possible command states: "start", "stop", "NA".)
#    SCHEDULER.running (JSON structure of PID, scheduler_id, mode, start_date, end_date if running.  "NA" if not running.)
#    SCHEDULER.stdout    (line-by-line of output from mmtscheduler.py)
#    SCHEDULER.queue   (JSON-format schedule)
#    SCHEDULER.status   (Single word for status:  'running'|stopped)
#    SCHEDULER.dispatcher   (JSON-format of dispatcher output)
#
# Current plans are to run it as a systemd service on a new, fast ("number crunching") computer.
#  command = {}
#
import redis
import time
import sys
import subprocess
import os
import datetime
import signal
import json

stuff = {'scheduler_id': 29,
            'command': 'start',
            'mode': 2,
            'start_date': '2017-01-10T19:00:00',
            'end_date': '2017-01-18T19:00:00',
            'scheduler_pid': -1,
            'computation_start': time.strftime('%Y-%m-%d %H:%M:%S')}

# print json.dumps(stuff)

print "publish SCHEDULER.command '" + json.dumps(stuff) + "'"

output = """
publish SCHEDULER.command '{"scheduler_id": 29, "end_date": "2017-01-10T19:00:00", "scheduler_pid": -1, "computation_start": "2017-02-15 11:27:07", "command": "start", "mode": 2, "start_date": "2017-01-10T19:00:00"}'
"""
2018-04-03:  Temporarily running the scheduler.service and dispatcher.service on ops2.  Using "/usr/bin/python3.6" for the version of python.
Used "/usr/bin/pip3.6" to install astropy 3, numpy 1.14, mkl, redis, pytz, and requests.  The default version of python3 on CentOS 7 is python3.4.  Astropy 3.0 requires python 3.5 or higher.

The two services here, scheduler.service and dispatcher.service, are copied to /lib/systemd/system/ rather than using symbolic links.  "systemctl enable *.service" will not work with symlinks.  Make sure that the versions of the services in this directory are copied to /lib/systemd/system/.


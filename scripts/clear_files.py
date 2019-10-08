""" A daemon to handle file cleanup on various locations."""

# a useful thing is to use find with -mmin and -delete, e.g.

# find /mnt/ssd/fildata -type f -mmin +59 -delete -name "*.fil"

# (tested this command, works like a charm.)

# should delete all files older than an hour.  This could be run with a cron
# job or else a persistent python script that runs every now and then.

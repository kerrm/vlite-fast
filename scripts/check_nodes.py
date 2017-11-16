#!/usr/bin/python
""" Periodically check on all of the hosts in the provided configuration
    file.  Send email if one of the servers drops out."""

import time
import os

# poll nodes every interval seconds
interval = 60

# alert these email addresses
emails = ['matthew.kerr@gmail.com','matthew.kerr@nrl.navy.mil',
    'aerickso@nrao.edu']

conf = '/home/vlite-master/mtk/config/hosts'
import sys
if len(sys.argv) > 1:
  conf = sys.argv[1]
  print 'Using configuration file %s.'%conf

def good_line(l):
  l = l.strip()
  if l[0] == '#':
    return False
  return True

def host_sort_key(h):
  toks = h.split('.')
  return int(toks[0][10:])

def send_email(host):
  print 'node %s is unresponsive, sending email'%(host)
  cmd = 'mailx -s "%s is unresponsive" %s < /dev/null > /dev/null'%(
      host,' '.join(emails))
  os.system(cmd)

good_lines = filter(good_line,file(conf).readlines())
hosts = sorted(list(set([x.split()[0] for x in good_lines])),
    key=host_sort_key)

server_ok = [True]*len(hosts)
# NB this is redundant...
message_sent = [False]*len(hosts)

while(True):
  for ihost,host in enumerate(hosts):
    #print 'checking host %s'%host
    cmd = 'ping -c 3 %s > /dev/null 2>&1'%host
    ok = os.system(cmd) == 0
    if (not ok):
      if server_ok[ihost]:
        # server has gone down!  send a message
        send_email(host)
        server_ok[ihost] = False
        message_sent[ihost] = True
      else:
        # server is still down -- do nothing
        pass
    else:
      if not server_ok[ihost]:
        # server has come back up
        print 'node %s is responsive again'%host
        server_ok[ihost] = True
        message_sent[ihost] = False
      else:
        # server was up and is up, do nothing
        pass
  time.sleep(interval)


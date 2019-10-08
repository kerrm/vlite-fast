#!/home/vlite-master/mtk/bin/python
""" Get a report on available disk space.
"""

import subprocess
hosts = [1,2,3,4,5,6,7,8,9,10,11,12]
hosts = ['vlite-difx%d'%i for i in hosts]

def parse_ssd_result(output):
    """ Return a dictionary with key usage entries."""
    results = dict()
    results['fildata'] = 0
    results['dumps'] = 0
    results['total'] = 0
    lines = output.split('\n')
    for line in lines:
        toks = line.split()
        if len(toks) == 2:
            if toks[1] == '/mnt/ssd/fildata':
                results['fildata'] = float(toks[0])
            elif toks[1] == '/mnt/ssd/dumps':
                results['dumps'] = float(toks[0])
            elif toks[1] == '/mnt/ssd':
                results['total'] = float(toks[0])
    return results

def print_results(results):
    hosts = sorted(results.keys(),key=lambda h: int(h[10:]))
    for host in hosts:
        hresults = results[host]
        print '\n'.join((
            '==========================',
            'Host: %s'%host,
            '  Filterbanks: %03.2f GB'%(hresults['fildata']/1e6),
            '  Voltages:    %03.2f GB'%(hresults['dumps']/1e6),
            '  Total:       %03.2f GB'%(hresults['total']/1e6)))

results = dict()
for host in hosts:
    cmd = 'ssh mkerr@%s "du -c /mnt/ssd"'%host
    try:
        output = subprocess.check_output(cmd,shell=True,
                stderr=subprocess.STDOUT)
    # permissions errors on lost+found etc. cause retcode!=0
    except subprocess.CalledProcessError as e:
        output = e.output
    results[host] =  parse_ssd_result(output)

print_results(results)

# TODO -- make a histogram plot using cumulative breakdowns.  This can be
# pushed to the monitoring website.

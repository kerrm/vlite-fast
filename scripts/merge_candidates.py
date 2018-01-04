#!/usr/bin/python
import glob
import os
from collections import deque,defaultdict


# read in all non-merged candidates
#fnames = sorted(glob.glob('/home/vlite-master/mtk/data/*.cand.000'))
#stems = [fname.rstrip('.000') for fname in fnames]

# set up a dictionary with the matching filenames
all_fnames = sorted(glob.glob('/home/vlite-master/mtk/data/*.cand.[0-9][0-9][0-9]'))
fname_dict = defaultdict(deque)
for fname in all_fnames:
  stem = fname[:len(fname)-4]
  fname_dict[stem].append(fname)

remove_original = True

numfiles = len(fname_dict.keys())

for istem,stem in enumerate(sorted(fname_dict.keys())):
  if istem % 100 == 0:
    print 'Working on %s (%04d/%04d).'%(stem,istem+1,numfiles)
  # avoid this I/O
  #fnames = sorted(glob.glob('%s.[0-9][0-9][0-9]'%stem))
  fnames = fname_dict[stem]
  all_lines = deque()
  for fname in fnames:
    lines = map(str.strip,file(fname).readlines())
    lines = [x for x in lines if len(x) > 0]
    if len(lines) > 0:
      all_lines.extend(lines)
  f = file(stem,'w')
  f.write('\n'.join(all_lines))
  f.write('\n')
  f.close()
  if remove_original:
    cmd = 'rm -f %s'%(' '.join(fnames))
    os.system(cmd)

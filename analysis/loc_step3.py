""" Using delays measured in Step3, along with metadata about the VLITE
    antenna positions and delays, use TDOA algorithm to localize the
    pulse source.
"""
# TODO -- restore astropy units

import filecmp
import glob
import os
import pickle
import xml.dom.minidom
from xml.etree import ElementTree

from astropy.coordinates import SkyCoord
from astropy.coordinates.builtin_frames import ICRS,ITRS
from astropy.time import Time
import astropy.units as u
import healpy
import numpy as np
import pylab as pl

C = 299792458 
tsamp = 1./128e6

def get_vla_center():
    #define VLA_CENTER_X                    -1601185.4      /* [m] */
    #define VLA_CENTER_Y                    -5041977.5      /* [m] */
    #define VLA_CENTER_Z                    3554875.9       /* [m] */
    x0 = -1601185.4
    y0 = -5041977.5
    z0 =  3554875.9
    return np.asarray([x0,y0,z0])

def parse_antprop_xml(fname):
    """
    Parse antenna names and positions out of saved XML data.

    Returns
    -------
    Dict with antenna names as keys and X,Y,Z position tuples as values.
    """

    rvals = dict()
    e = ElementTree.fromstring(open(fname).read())
    for elem in e:
        if elem.tag == 'AntennaProperties':
            name = elem.attrib['name']
            X = float(elem.find('X').text)
            Y = float(elem.find('Y').text)
            Z = float(elem.find('Z').text)
            rvals[name] = np.asarray([X,Y,Z])
    return rvals

def cmp_antprop(ap1,ap2):
    for key in ap1.keys():
        if not key in ap2:
            return False
        if not np.all(ap1[key]==ap2[key]):
            return False
    return True

def get_antenna_positions(t0):
    """ Find the antenna prop file closest in time before t0."""
    antprops = glob.glob('localization_data/antprop/antprop_1*.xml')
    times = np.asarray([int(os.path.split(x)[-1].split('_')[-1][:-4]) for x in antprops])
    a = np.argmin(np.abs(times-t0))
    if times[a] > t0:
        if a == 0:
            print('Warning, t0 precedes earliest antprop file.')
        else:
            a -= 1
    fname = antprops[a]
    epoch = fname[-14:-4]
    rvals_follow = None
    print('Using %s for antenna positions (%s).'%(
        os.path.basename(fname),Time(epoch,format='unix').datetime))
    rvals = parse_antprop_xml(fname)
    if a == (len(times)-1):
        print('Warning, t0 follows latest antrop file.')
    else:
        follow_fname = antprops[a+1]
        rvals_follow = parse_antprop_xml(follow_fname)
        if not cmp_antprop(rvals_follow,rvals):
            follow_epoch = follow_fname[-14:-4]
            for key in rvals.keys():
                if key in rvals_follow.keys():
                    diff = (np.asarray(rvals[key])-np.asarray(rvals_follow[key]))
                    if not np.all(diff==0):
                        print('Possible position shift for %s : '%(key),diff)
            #print('Warning, antenna positions in %s (%s) differ.'%(
                #os.path.basename(follow_fname),
                #Time(follow_epoch,format='unix').datetime))
    x0 = get_vla_center()
    for key in rvals.keys():
        rvals[key] += x0
    return rvals

def parse_correlator_delays(t0):
    """ Parse the "vliteantennas.in" file for delays from the correlator.

    Returns
    ------
    Dict keyed on antennas with delay values (in nanoseconds).
    """

    # find nearest preceding meta data file.  For concreteness, compare to
    # the next one to see if there are any changes and alert to them.

    metas = sorted(glob.glob('localization_data/meta/META*'))
    t0s = [meta[-10:] for meta in metas]
    idx = np.searchsorted(t0s,t0)
    first = last = False
    if idx == 0:
        print('Warning, t0 occurs before first META file date.')
        first = True
    elif idx == len(metas):
        print('Warning, t0 occurs after last META file date.')
        last = True
    if not first:
        idx -= 1
    fname = metas[idx]
    print('Selecting META file %s (%s).'%(
        os.path.basename(fname),Time(t0,format='unix').datetime))
    if not (first or last):
        # compare bracketing files and warn if they are different
        if not filecmp.cmp(fname,metas[idx+1]):
            print('Warning, subsequent META file %s differs.'%metas[idx+1])

    # format of a line
    # 0 10 vlite-difx7 p1p2 5021.000000 E08 5021.000000 1
    lines = [x.strip() for x in open(fname).readlines()]
    lines = [x for x in lines if not x.startswith('#')]
    delays = dict()
    antenna_map = dict()
    for line in lines:
        try:
            toks = line.split()
            vant = 'V%d'%(int(toks[0]))
            antenna = 'ea%02d'%(int(toks[1]))
            delay = float(toks[4])
            delays[antenna] = delay
            antenna_map[vant] = antenna
        except KeyError:
            pass
    return delays,antenna_map

def load_pipeline_delays(t0,antenna_map,pol=0):
    """ Parse tabulated delays from VLITE pipeline to find best-matching
        delays for epoch.

    Returns
    -------
    Dict keyed on antennas with delay values (nanoseconds).
    """

    in_delays = pickle.load(
            open('localization_data/vlite_delays/vlite_delays.pickle','rb'))
    delays = dict() 
    mjd = Time(t0,format='unix').mjd
    for vant in antenna_map.keys():
        ant = antenna_map[vant] 
        if vant == 'V0':
            delays[ant] = [0,0]
        else:
            fix_key = 'V%02d'%(int(vant[1:]))
            try:
                d = in_delays[fix_key]
                idx = np.searchsorted(d['edges_mjd'],mjd)
                if (idx==len(d['edges_mjd'])):
                    print('Warning, extrapolating delays for antenna %s.'%(ant))
                idx -= 1
                pxd = d['polx_delays'][idx]
                pyd = d['poly_delays'][idx]
                #delays[ant] = [pxd,pyd]
                # ACHTUNG -- NB that thread=0 == pol Y (I think...)
                delays[ant] = [pyd,pxd]
            except KeyError:
                print('Warning, no VLITE pipeline delay for %s/%s.'%(ant,vant))

    return delays

def load_delays(t0,valid_antennas=None,select_pol=None):
    measured_antennas = set()
    outdir = 'localization_output/%.3f'%(t0)
    measured_delays = pickle.load(open(os.path.join(outdir,'step2_test.pickle'),'rb'))
    rvals = dict()
    mask = dict()
    for baseline in measured_delays.keys():
        a1,a2,this_pol = baseline.split('-')
        if select_pol is not None:
            if int(this_pol[-1]) != select_pol:
                continue

        idx0,x0,err,phase,sn,coh1,coh2,viz = measured_delays[baseline] 

        if (valid_antennas is not None):
            if (a1 not in valid_antennas) or (a2 not in valid_antennas):
                continue
        measured_antennas.add(a1)
        measured_antennas.add(a2)

        rvals[baseline] = (idx0+x0)*tsamp
        mask[baseline] = sn > 10

    return measured_antennas,rvals,mask

#t0 = 1584034753.161
#t0 = 1584408723.204
t0 = 1581646177.340

outdir = 'localization_output/%.3f'%(t0)

# load in antenna positions -- generally a superset of vlite antennas
positions = get_antenna_positions(t0)

# load in "correlator" delay offsets; rough estimates of delays
correlator_delays,antenna_map = parse_correlator_delays(t0)

# load in additional delays from the VLITE pipeline
pipeline_delays = load_pipeline_delays(t0,antenna_map,pol=0)

# form a set of antennas with good delay measurements
valid_antennas = set(positions.keys())
valid_antennas = valid_antennas.intersection(correlator_delays.keys())
valid_antennas = valid_antennas.intersection(pipeline_delays.keys())

# manually excise any antennas
bad_antennas = []
bad_antennas = ['ea03']
valid_antennas = valid_antennas.difference(bad_antennas)

# manually add in delay for ea05
#pipeline_delays['ea05'][0] += 3.01
#pipeline_delays['ea05'][1] += 3.01



# now load in previously computed delays
measured_antennas,measured_delays,baseline_mask = load_delays(
        t0,valid_antennas=valid_antennas)
valid_antennas = list(valid_antennas.intersection(measured_antennas))

# make total set of "systematic" delays for each antenna
clock_delays = dict()
for ant in valid_antennas:
    clock_delays[ant] = [correlator_delays[ant] + pipeline_delays[ant][0],
                         correlator_delays[ant] + pipeline_delays[ant][1]]


# observation specific information
obstime = Time(t0,format='unix',scale='utc')
# Crab
#pos_trans = SkyCoord('05:34:31.93830','+22:00:52.1758',frame=ICRS,unit=[u.hourangle,u.deg])
#pos_trans = SkyCoord.from_name('PSR B0329+54')
pos_trans = SkyCoord.from_name('PSR J0534+2200')
# 3C 147
#pos_point = SkyCoord('05:42:36.13789843','+49:51:07.2337251',unit=[u.hourangle,u.deg])
XYZ_trans= pos_trans.cartesian.get_xyz().value

# form ICRS vectors for each antenna using this hack
antennas_icrs = dict()
XYZ_icrs = dict()
for ant in valid_antennas:
    itrs_ant = ITRS(positions[ant],obstime=obstime)
    icrs_ant = itrs_ant.transform_to(ICRS)
    antennas_icrs[ant] = icrs_ant
    distance = np.sum(positions[ant]**2)**0.5*u.m
    XYZ_icrs[ant] = icrs_ant.cartesian.get_xyz().value*distance

# project baselines
baselines = list(measured_delays.keys())
projected_baselines = dict()
delays_baselines = dict()
for baseline in baselines:
    a1,a2,pol = baseline.split('-')
    pol = int(pol[-1])
    projected_baselines[baseline] = XYZ_icrs[a1]-XYZ_icrs[a2]
    delays_baselines[baseline] = clock_delays[a1][pol]-clock_delays[a2][pol]

# form arrays for rapid computation of delays / chi^2 later
projected_baselines_array = np.asarray([projected_baselines[baseline] for baseline in baselines])/C
clock_delays_array = np.asarray([delays_baselines[baseline] for baseline in baselines])*1e-9
measured_delays_array = np.asarray([measured_delays[baseline] for baseline in baselines])
baseline_mask_array = np.asarray([baseline_mask[baseline] for baseline in baselines])


# make a sky grid
nside = 2**6
ras_hp = np.empty(healpy.nside2npix(nside))
decs_hp = np.empty(healpy.nside2npix(nside))
for i in range(healpy.nside2npix(nside)):
    ang = healpy.pix2ang(nside,i,lonlat=True)
    ras_hp[i] = ang[0]
    decs_hp[i] = ang[1]
sc_hp = SkyCoord(ras_hp,decs_hp,frame=ICRS,unit=[u.deg,u.deg])
XYZ_hp = sc_hp.cartesian.get_xyz().value.transpose()

# make a zoomed in sky grid
ngrid = 200
if ngrid % 2 == 1:
    ngrid += 1
X,Y = np.meshgrid(np.arange(ngrid+1),np.arange(ngrid+1))
rarcmin = 10
zoom_ra_grid = pos_trans.ra.value + np.linspace(-rarcmin/60,rarcmin/60,ngrid+1)*0.5/np.cos(np.radians(pos_trans.dec.value))
zoom_dec_grid = pos_trans.dec.value + np.linspace(-rarcmin/60,rarcmin/60,ngrid+1)*0.5/np.cos(np.radians(pos_trans.dec.value))
zoom_ras = zoom_ra_grid[X]
zoom_decs = zoom_dec_grid[Y]
zoom_sc = SkyCoord(zoom_ras,zoom_decs,frame=ICRS,unit=[u.deg,u.deg])
XYZ_zoom = zoom_sc.cartesian.get_xyz().value.transpose()

# compute delays on the grids, and for the exact source position
tau_hp = np.einsum('ij,kj->ik',XYZ_hp,projected_baselines_array)
tau_zoom = np.einsum('ijk,lk->ijl',XYZ_zoom,projected_baselines_array)
tau_trans = (XYZ_trans*projected_baselines_array).sum(axis=1)

# calculate residuals -- TODO -- fix signs to be more sensible
total_delays_ns = measured_delays_array.copy()*1e9
total_delays_ns += clock_delays_array*1e9
chi_hp = total_delays_ns[None,:]+tau_hp*1e9
# we do *not* want to subtract mean, because this is a direct prediction
# between the two antennas, does not need to be referenced to anything
#chi -= chi.mean(axis=1)[:,None]

chi_zoom = total_delays_ns[None,None,:] + tau_zoom*1e9
#chi_zoom -= chi_zoom.mean(axis=-1)[:,:,None]

"""
# what we can do, though, is an antenna-dependent delay to attempt to remove
# any offsets
#antenna_indices = dict()
for antenna in valid_antennas:
    if antenna != 'ea05':
        continue
    qpos = np.asarray([antenna==baseline.split('-')[0] for baseline in baselines])
    qneg = np.asarray([antenna==baseline.split('-')[1] for baseline in baselines])
    sign = np.zeros(len(qpos))
    sign[qpos] = 1
    sign[qneg] = -1
    mask = sign != 0
    #antenna_indices[antenna] = sign
    q = (chi[...,mask]*sign[mask]).mean(axis=-1)
    chi[...,mask] += q[...,None]
    q = (chi_zoom[...,mask]*sign[mask]).mean(axis=-1)
    chi_zoom[...,mask] += q[...,None]
"""

rms_hp = (chi_hp[:,baseline_mask_array]**2).mean(axis=1)**0.5
rms_zoom = (chi_zoom[:,:,baseline_mask_array]**2).mean(axis=-1)**0.5

# plot residuals
pl.close('all')
healpy.mollview(np.log10(rms_hp),coord='C',flip='astro')
healpy.graticule()
default = dict(marker='x',color='red',markersize=8)
default['lonlat'] = True
pl.gca().projplot(pos_trans.ra,pos_trans.dec,**default)
pl.savefig(os.path.join(outdir,'step3_skymap.png'))
#healpy.mollview(np.log10(1./chi2),coord='C',flip='astro')
#healpy.mollview(rms,coord='C',flip='astro')
#pl.axis([-0.9557,-0.7784,0.1827,0.3945])
#healpy.graticule()
#default['marker'] = '*'
#pl.gca().projplot(crab.ra.value,crab.dec.value,**default)
pl.figure(2); pl.clf()
#pl.pcolor(zoom_ras,zoom_decs,rms_zoom)#,edgecolor='k')
#pl.imshow(rms_zoom.transpose(),extent=[,zoom_decs,rms_zoom)#,edgecolor='k')
dra = zoom_ra_grid[1]-zoom_ra_grid[0]
ra0 = zoom_ra_grid[0] - 0.5*dra
ra1 = zoom_ra_grid[-1] + 0.5*dra
ddec = zoom_dec_grid[1]-zoom_dec_grid[0]
dec0 = zoom_dec_grid[0] - 0.5*ddec
dec1 = zoom_dec_grid[-1] + 0.5*ddec
pl.imshow(rms_zoom.transpose(),origin='lower',extent=[ra0,ra1,dec0,dec1],aspect='auto')
pl.colorbar()
idx = np.argmin(np.ravel(rms_zoom))
xidx = idx//len(zoom_ras)
yidx = idx-len(zoom_ras)*xidx
best_loc = SkyCoord(zoom_ra_grid[xidx],zoom_dec_grid[yidx],unit=[u.deg,u.deg])
#resids = chi_zoom[xidx,yidx][baseline_mask_array]
# TMP -- replace with known position to estimate error
resids = chi_zoom[ngrid//2,ngrid//2][baseline_mask_array]
#std = np.median(np.abs(resids-np.median(resids)))*1.4826
# don't correct for mean
std = np.median(np.abs(resids))*1.4826
resids = chi_zoom[xidx,yidx][baseline_mask_array]
# set contours based on chi^2 drop
yvals = rms_zoom**2/std**2
# make it such that minimum value has chi^2/dof = 1
chisq_dof = yvals.min()
dof = baseline_mask_array.sum()
#yvals *= dof/chisq_dof
yvals *= dof
# draw contours for likelihood decreases by given amounts -- calibrate with
# the Rayleigh distribution once I remember the numbers...
yvals -= yvals.min()
pl.contour(zoom_ras,zoom_decs,yvals.transpose(),levels=[1,2.5,6],colors='white',linestyles='--')
# don't use std to compute trans_rms, since it subtracts the mean, which
# we don't wanna do
#trans_rms = (measured_delays_array+clock_delays_array+tau_trans)[baseline_mask_array]
#trans_rms = (trans_rms**2).mean()**0.5*1e9
# this should also be
trans_rms = rms_zoom[ngrid//2,ngrid//2]
# 
print('Separation=%.2f arcmin, fit_rms = %.2f ns, trans_rms = %.2f ns'%(best_loc.separation(pos_trans).value*60,rms_zoom.min(),trans_rms))
pl.plot(pos_trans.ra.value,pos_trans.dec.value,marker='*',markersize=8,color='red',label='Known Pos r.m.s = %.2f ns'%trans_rms)
pl.plot(best_loc.ra.value,best_loc.dec.value,marker='^',label='Best Fit r.m.s. = %.2f ns, $\chi^2_{%d}$ = %.1f'%(rms_zoom.min(),dof,round(chisq_dof*dof)))
pl.legend(loc='lower left')
pl.xlabel('Right Ascension')
pl.ylabel('Declination')
pl.savefig(os.path.join(outdir,'step3_zoom.png'))
#pl.title('Global Best r.m.s. = %.2f ns'%(rms_zoom.min()))
#pl.figtext(0.05,0.95,'Best r.m.s. = %.2f ns'%(rms_zoom.min()))
#pl.figtext(0.05,0.90,'Source r.m.s. = %.2f ns'%(trans_rms))

# make a nice grid plot of baselines
va = sorted(valid_antennas)
resids = (measured_delays_array+clock_delays_array+tau_trans)*1e9
#resids = chi_zoom[xidx,yidx]
grid_resids = np.empty((2*len(va),len(va)))
grid_resids[:] = np.nan
for ibaseline,baseline in enumerate(baselines):
    ant1,ant2,pol = baseline.split('-')
    resid = resids[ibaseline]
    idx1 = va.index(ant1)
    idx2 = va.index(ant2)
    if int(pol[-1]) != 0:
        grid_resids[idx1+len(va),idx2] = resid
        grid_resids[idx2+len(va),idx1] = resid
    else:
        grid_resids[idx1,idx2] = resid
        grid_resids[idx2,idx1] = resid

pl.close(3);
pl.figure(3,(8,4))
vmax = np.abs(resids[baseline_mask_array]).max()
vmin = -vmax
pl.imshow(grid_resids.transpose(),interpolation='nearest',aspect='auto',vmin=vmin,vmax=vmax,origin='lower',cmap='PiYG')
pl.colorbar()
pl.xticks(np.arange(2*len(va)),labels=va*2,rotation='vertical')
pl.yticks(np.arange(len(va)),labels=va,rotation='horizontal')
pl.axvline(len(va)-0.5,color='k')
pl.tight_layout()
pl.subplots_adjust(right=1.05,top=0.92)
pl.figtext(0.25,0.94,'pol=0',size='large')
pl.figtext(0.65,0.94,'pol=1',size='large')
pl.savefig(os.path.join(outdir,'step3_baseline_diagnostic.png'))

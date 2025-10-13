"""
FILE
    flexce.py

DESCRIPTION

    Flexible Galactic Chemical evolution code.
    
"""

version = 0.0

import os
import sys
home = os.path.expanduser('~')
stem_flexce = home + '/FlexCE_v' + str(version) + '/'
stem_scripts = stem_flexce + 'scripts/'
sys.path.append(stem_scripts)

import numpy as np
import pickle

import flexce_main
from flexce_main import ChemEvolZone

import import_yields
from import_yields import Yields

import calc_abundances
from calc_abundances import Abundances

from astropy.io import fits

try:
    file_in = sys.argv[1]
except:
    print('Usage: %s parameter_file_name' % sys.argv[0])
    sys.exit(1)



##### Set path #####
stem_runs = stem_flexce + 'runs/'
stem_param = stem_flexce + 'param_files/'
sim_id = file_in.strip('param').strip('.txt')
####################


##### Read in parameter file #####
dir_args = {}
yld_args = {}
initialize_args = {'sim_id': sim_id}
mass_bins_args = {}
snia_dtd_args = {}
inflows_args = {}
outflows_args = {}
warmgasres_args = {}
sf_args = {}

f = open(stem_param + file_in, 'r')
for line in f:
    if (line[0] != '#') and (line != '\n'):
        if ' = ' in line:
            k = line.strip().split(' = ')[0]
            v = line.strip().split(' = ')[1]
            try:
                v = float(v)
            except:
                pass
            try:
                if v[0] == '[':
                    v = v.strip('[').strip(']').split(', ')
                    try:
                        v = [float(item) for item in v]
                    except ValueError:
                        if v ==['']:
                            v = []
            except:
                pass
            if v == 'True':
                v = True
            elif v == 'False':
                v = False
        if 'yields' in k:
            if (',' in v) and (v[-1] != ','):
                yld_args[k.split('yields_')[1]] = v.split(', ')
            elif (',' in v) and (v[-1] == ','):
                 yld_args[k.split('yields_')[1]] = [v.strip(',')]
            else:
                yld_args[k.split('yields_')[1]] = v
        elif 'initialize' in k:
            initialize_args[k.split('initialize_')[1]] = v
        elif 'mass_bins' in k:
            mass_bins_args[k.split('mass_bins_')[1]] = v
        elif 'snia_dtd' in k:
            if k.split('snia_dtd_')[1] == 'func':
                snia_dtd_args[k.split('snia_dtd_')[1]] = v
            else:
                if 'kwargs' not in snia_dtd_args:
                    snia_dtd_args['kwargs'] = {}
                snia_dtd_args['kwargs'][k.split('snia_dtd_')[1]] = v
        elif 'inflows' in k:
            if k.split('inflows_')[1] in ['M1', 'M2', 'b1', 'b2']:
                if 'k' not in inflows_args:
                    inflows_args['k'] = {}
                inflows_args['k'][k.split('inflows_')[1]] = v
            else:
                inflows_args[k.split('inflows_')[1]] = v
        elif 'outflows' in k:
            outflows_args[k.split('outflows_')[1]] = v
        elif 'warmgasres' in k:
            warmgasres_args[k.split('warmgasres_')[1]] = v
        elif 'sf_' in k:
            sf_args[k.split('sf_')[1]] = v
        elif 'sfh_file' in k:
            sfh_file = v
            
#inflows_args
#print inflows_args

f.close()

# Load the inflow file
if inflows_args['func'] == 'custom':
    print("Reading inflow rates from = " + inflows_args['file'])
    f = open(inflows_args['file'],'r')
    inflow_rate = []
    for line in f:
        inflow_rate.append(float(line))
    f.close()
    #print inflow_rate
    inflows_args['inflow_rate'] = inflow_rate
    #print inflows_args

# Load the SFE file
if 'sfe_file' in sf_args:
    print("Reading SFE from = " + sf_args['sfe_file'])
    f = open(sf_args['sfe_file'],'r')
    sfe = []
    for line in f:
        sfe.append(float(line))
    f.close()
    #print sfe
    sf_args['sfe'] = sfe
    #print sf_args
else:
    n_steps = int(float(initialize_args['time_tot']) / float(initialize_args['dt']) + 1.)
    sf_args['sfe'] = np.zeros(n_steps) + float(sf_args['nu_kslaw'])

# Load the SFH file
if 'sfh_file' in locals():
    print("Reading SFH from = " + sfh_file)
    f = open(sfh_file,'r')
    sfh = []
    for line in f:
        sfh.append(float(line))
    f.close()
    #print sfh
    
#print sf_args

####################################

##### Stellar Mass Bins #####
mass_bins = np.concatenate(
    (np.arange(mass_bins_args['low'], 8., mass_bins_args['dm_low']),
     np.arange(8., mass_bins_args['high'] + 0.1, mass_bins_args['dm_high'])))
#############################

##### Yields #####
yld = Yields(stem_flexce, snii_dir=yld_args['snii_dir'],
             agb_dir=yld_args['agb_dir'], snia_dir=yld_args['snia_dir'],
             rprocess_dir=yld_args['rprocess_dir'],
             sprocess_dir=yld_args['sprocess_dir'], mbins=mass_bins)
yld.load_rprocess_yields()
yld.load_sprocess_yields()
yld.load_snii_yields()
yld.load_agb_yields()
yld.load_snia_yields(model=yld_args['snia_model'])
yld.concat_ncapture_yields(r_elements=yld_args['r_elements'],
                           s_elements=yld_args['s_elements'])
yld.load_solar_abund()
#################

##### Evolve Box #####
box = ChemEvolZone(yld.mass_bins, **initialize_args)
box.snia_dtd(**snia_dtd_args)
box.inflow_rx(**inflows_args)
box.outflow_rx(**outflows_args)
box.warmgasres_rx(**warmgasres_args)
box.star_formation(**sf_args)
if 'sfh' in locals():
    box.evolve_box(yields=yld,sfh=sfh)
else:
    box.evolve_box(yields=yld)
######################

##### Abundances #####
ab = Abundances(stem_flexce, yld.sym, box.mgas_iso, box.survivors, box.t,
                box.param, box.sim_id)
ab.load_solar_abund()
ab.calc_abundances()

apogee_el = np.array(['C', 'N', 'O', 'Na', 'Mg', 'Al', 'Si', 'S',
                      'K', 'Ca', 'Ti', 'V', 'Cr', 'Mn', 'Co', 'Ni'])
ab.select_elements(apogee_el)
#####################
 
##### Save simulation runs #####
fout = open(stem_runs + 'box' + sim_id + '.pck', 'w')
pickle.dump(box, fout, 2)
fout.close()

# same some info
#dt = np.dtype([('agb',float),('dm_sfr',float),('inflow',float),('inflow_rate',float),
#               ('mass_ave',float),('mass_ave2',float),('mass_bins',float),('mass_bins2',float),
#               ('mass_frac',float),('mass_frac2',float),('mass_int',float),('mass_int2',float),
#               ('metallicity',float),('mfrac',float),('mgas_iso',float),('mremnant',float),('mstar',float),
#               ('mstar_left',float),('mstar_stat',float),('mwarmfrac',float),('mwarmgas_iso',float),
#               ('outflow',float),('sf',float),('sfe',float),('sfr',float),('snia',float),('snii',float),
#               ('survivors',int),('t',float),('tau_m',float)])
#data = np.zeros(len(box.t),dtype=dt)
#data['agb'] = box.agb
#data['dm_sfr'] = box.dm_sfr
#data['inflow'] = box.inflow
#data['inflow_rate'] = box.inflow_rate
#data['mass_ave'] = box.mass_ave
#data['mass_ave2'] = box.mass_ave2
#data['mass_bins'] = box.mass_bins
#data['mass_bins2'] = box.mass_bins2
#data['mass_frac'] = box.mass_frac
#data['mass_frac2'] = box.mass_frac2
#data['mass_int'] = box.mass_int
#data['mass_int2'] = box.mass_int2
#data['metallicity'] = box.metallicity
#data['mfrac'] = box.mfrac
#data['mgas_iso'] = box.mgas_iso
#data['mremnant'] = box.mremnant
#data['mstar'] = box.mstar
#data['mstar_left'] = box.mstar_left
#data['mstar_stat'] = box.mstar_stat
#data['mwarmfrac'] = box.mwarmfrac
#data['mwarmgas_iso'] = box.mwarmgas_iso
#data['outflow'] = box.outflow
#data['sf'] = box.sf
#data['sfe'] = box.sfe
#data['sfr'] = box.sfr
#data['snia'] = box.snia
#data['snii'] = box.snii
#data['survivors'] = box.survivors
#data['t'] = box.t
#data['tau_m'] = box.tau_m

# some as 401 while others are 401x300, not sure why

fout = open(stem_runs + 'ab' + sim_id + '.pck', 'w')
pickle.dump(ab, fout, 2)
fout.close()

# Output to FITS file
dt = np.dtype([('t',float),('FeH',float),('CFe',float),('NFe',float),('OFe',float),('NaFe',float),
               ('MgFe',float),('AlFe',float),('SiFe',float),('PFe',float),('SFe',float),('KFe',float),
               ('CaFe',float),('TiFe',float),('VFe',float),('CrFe',float),('MnFe',float),('CoFe',float),
               ('NiFe',float)])
n = len(ab.xfe[0])
data = np.zeros(n,dtype=dt)
# xfe or xfe_all??
elem = ['C','N','O','Na','Mg','Al','P','S','K','Ca','Ti','V','Cr','Mn','Co','Ni']
for e in elem:
    ind = np.where(ab.elements == e)[0][0]
    data[e+'Fe'] = ab.xfe_all[ind]
data['t'] = ab.t[0:n]
data['FeH'] = ab.feh
outfile = stem_runs + 'ab' + sim_id + '.fits'
fits.writeto(outfile,data,overwrite=True)

################################

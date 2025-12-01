"""Main script.

Command line args:
    config file that starts with "sim".

Usage:
    python flexce.py ../config/sim0.cfg
"""

from __future__ import print_function, division, absolute_import

import os
from os.path import join
import sys
import numpy as np
from astropy.io import fits

from .fileio.cfg_io import read_sim_cfg
from .fileio.pickle_io import pickle_write
from .fileio.txt_io import txt_write
from . import utils
from .chemevol import ChemEvol
from .abundances import Abundances


path_flexce = join(os.path.abspath(os.path.dirname(__file__)), '')
path_data = join(path_flexce, 'data')
default_config_path = join(path_flexce, 'config')
default_config = read_sim_cfg(os.path.join(default_config_path,'sim0.cfg'))

def evolve(yld, initialize_kws, snia_dtd_kws, inflows_kws, outflows_kws,
           warmgasres_kws, sf_kws):
    """Evolve the galaxy.

    Args:
        yld: Yields instance
        initialize_kws (dict): args to initialize instance of ChemEvol class.
        mass_bins_args (dict): args to define stellar mass bins.
        snia_dtd_kws (dict): args to set SNIa delay time distribution of
            ChemEvol instance.
        inflows_kws (dict): args to set inflow rate and composition of
            ChemEvol instance.
        outflows_kws (dict): args to set outflow rate and composition of
            ChemEvol instance
        warmgasres_kws (dict): turn on warm ISM reservoir in ChemEvol
            instance.
        sf_kws (dict): args to set star formation rate in ChemEvol instance.

    Returns:
        ChemEvol instance
    """
    gal = ChemEvol(yld.mass_bins, **initialize_kws)
    gal.snia_dtd(**snia_dtd_kws)
    gal.inflow_rx(**inflows_kws)
    gal.outflow_rx(**outflows_kws)
    gal.warmgasres_rx(**warmgasres_kws)
    gal.star_formation(**sf_kws)
    gal.evolve_box(yields=yld)
    return gal


def calc_abundances(path, sym, mgas, survivors, time, parameters, sim_id):
    """Calculate abundances of box.

    Args:
        path (str): data directory.
        sym (array): Isotope abbreviations.
        mgas (array): Mass of each isotope in gas-phase at each timestep.
        survivors (array): Number of stars from each timestep that survive to
            the end of the simulation.
        time (array): time in Myr.
        parameters (dict): parameters of the simulation.
        sim_id (str): simulation ID number.

    Returns:
        Abundances instance
    """
    abund = Abundances(path, sym, mgas, survivors, time, parameters, sim_id)
    abund.load_solar_abund()
    abund.calc_abundances()
    #apogee_el = np.array(['C', 'N', 'O', 'Na', 'Mg', 'Al', 'Si', 'S',
    #                      'K', 'Ca', 'Ti', 'V', 'Cr', 'Mn', 'Co', 'Ni'])
    #abund.select_elements(apogee_el)
    return abund


def output(path, sim_id, gal, abund):
    """Write simulation results to pickle and txt files.

    Args:
        path (str): output directory.
        sim_id (str): simulation ID number.
        gal: ChemEvol instance.
        abund: Abundances instance.
    """
    path_sim = join(path, ''.join(['sim', sim_id]))
    if not os.path.isdir(path_sim):
        os.makedirs(path_sim)

    pickle_write(gal, join(path_sim, ''.join(('box', sim_id, '.pck'))))
    pickle_write(abund, join(path_sim, ''.join(('ab', sim_id, '.pck'))))

    txt_write(path_out, sim_id, gal, abund)

    # Output to FITS file
    dt = [('t',float),('FeH',float)]
    for i in range(len(abund.elements)):
        if abund.elements[i] != 'Fe':
        dt.append((abund.elements[i]+'Fe',float))
    n = len(abund.xfe[0])
    data = np.zeros(n,dtype=np.dtype(dt))
    for e in abund.elements:
        if e != 'Fe':
            ind = np.where(abund.elements == e)[0][0]
            data[e+'Fe'] = abund.xfe_all[ind] 
    data['t'] = abund.t[0:n]
    data['FeH'] = abund.feh
    outfile = join(path_sim, ''.join(('ab' + sim_id + '.fits')))
    fits.writeto(outfile,data,overwrite=True)
    

def run(config,path_out=None):
    
    mass_bins = utils.define_mass_bins(**config['mass_bins'])
    ylds = utils.load_yields(path_data, config['yields'], mass_bins)
    box = evolve(ylds, config['init'], config['dtd'], config['inflows'], config['outflows'],
                 config['warmgas'], config['sf'])
    ab = calc_abundances(path_data, ylds.sym, box.mgas_iso, box.survivors,
                         box.t, box.param, box.sim_id)
    return box,ab


    
# TODO (specify elements_out in config file)

if __name__ == '__main__':

    argv = None
    if argv is None:
        argv = sys.argv

    try:
        default_config_path = join(path_flexce_root, 'config')
        fname, path_config = utils.set_path(argv[1], default_config_path)
    except IndexError:
        path_config = join(os.getenv('HOME'), 'flexce', 'examples')
        fname = 'sim0.cfg'
        print('\nUsing default parameters in \n{}'.format(argv[1]))

    file_in = join(path_config, fname)

    # Load the config file
    config = read_sim_cfg(file_in)
    simulation_id = config['id']
    
    # TODO Add try...except to handle user-defined output path
    path_out = utils.substitute_dir_in_path(path_config, 'config', 'output')
    
    # run it
    box,ab = run(config,path_out)

    output(path_out, simulation_id, box, ab)

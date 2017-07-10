# Retrieval of species data from indices in binned coordinates
##### IMPORTS #####
import pygeode as pyg
import numpy as np
import rpnpy.librmn.all as rmn
import time, glob, math, operator, itertools
import matplotlib.pyplot as plt
import rpnpy.vgd.all as vgd
from mpl_toolkits.basemap import Basemap
from levels import *

##### CONSTANTS #####
pv_file_dir = '/home/ords/aq/alh002/pyscripts/workdir/pv_files'

##### FXNS #####
def build_data(species, ni, nj, nk, timesteps):
    '''
    builds data for the species identified based on binned coordinates

    ni: number of lon grid points
    nj: number of lat grid points
    species: species you wish to build the data for
    '''
    global pv_file_dir
    start_time = time.time()
    tropos = np.zeros((timesteps, nk, ni, nj))
    i = 0
    
    # reads nc file to retrieve data.
    # you may need to call splice_month in pv_calc.py if you're not dealing with january.
    species_data = pyg.open('/space/hall1/sitestore/eccc/aq/r1/alh002/NCDF/SPECIES/GO3/2009.nc')
    # flips levels so it's compatible with the data in tropo file
    species_data = species_data.go3[:timesteps,::-1]

    with open(pv_file_dir + '/tropo_coords.txt', 'r') as tropo_file:
        for line in tropo_file:
            # if the line read is a timestep barrier, retrieve the timestep for the following lines
            if line[:-1].split(',')[0] == 'TIMESTEP':
                timestep = line[:-1].split(',')[1]
                print "WORKING ON TIMESTEP: {0}".format(timestep)
                continue
            # otherwise, this will run
            coordinate = [int(x) for x in line[:-1].split(',')]  # creates a list with 3 ints [nk,ni,nj]
            for ind in coordinate:
                tropos[timestep, coordinate[0],coordinate[1],coordinate[2]] = species_data[timestep, coordinate[0],coordinate[2],coordinate[1]]
    print "That took {0} seconds.".format(str(time.time() - start_time))
    return tropos

##### MAIN #####
species = 'go3'  # the species data will be retrieved from is GEMS Ozone (from MACC Reanalysis)
tropo_go3 = build_data(species, 320, 73, 60, 28)


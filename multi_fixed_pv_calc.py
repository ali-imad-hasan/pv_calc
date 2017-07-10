# Calculating potential vorticity from fst file
##### IMPORTS #####
import pygeode as pyg
import numpy as np
import rpnpy.librmn.all as rmn
import time, glob, math, operator, itertools
import matplotlib.pyplot as plt
import rpnpy.vgd.all as vgd
from mpl_toolkits.basemap import Basemap
from levels import *
#const_pressure = nointerp_pressure
##### CONSTANTS #####
a       = 6.371229e6        #earth radius in m
g0      = 9.81
kn2ms   = 0.51444444     #conversion from knots to m/s
kappa   = 0.287514
dlatlon = 1.125
l       = 0             # used to name files
#filename = '/home/ords/aq/alh002/pyscripts/workdir/momentum_test.fst'
strato_file = '/home/ords/aq/alh002/pyscripts/workdir/pv_files/strato_coords.txt'
tropo_file = '/home/ords/aq/alh002/pyscripts/workdir/pv_files/tropo_coords.txt'
filenames  = ['/home/ords/aq/alh002/NCDF/testmom.nc']
lnsp_files = ["/home/ords/aq/alh002/NCDF/LNSP/2008.nc", 
              "/home/ords/aq/alh002/NCDF/LNSP/2009.nc", 
              "/home/ords/aq/alh002/NCDF/LNSP/2010.nc", 
              "/home/ords/aq/alh002/NCDF/LNSP/2011.nc", 
              "/home/ords/aq/alh002/NCDF/LNSP/2012.nc"]  # in this list, place the directory of the year
lnsp_files = [lnsp_files[1]]  # this is temporary, just to work with one file

##### FXN #####
def calc_qq(uu, vv, dx, dy, cosphi):
    '''
    calculates qq from u component of wind, v component of wind, dx, dy, and cosphi.

    uu/vv should be 4d arrays, dx/dy/cosphi should be 2d in the form nj, ni.

    returns 2d array in shape nj,ni
    '''

    # dvdx/dudy are calculated through leapfrogging technique
    dvdx = np.zeros((nj,ni))
    dudy = np.zeros((nj,ni))
    qq   = np.zeros((nj,ni))

    # calculation of qq from dvdx and dudy
    dvdx[:,0]  = (vv[t,k,:,1] - vv[t,k,:,0]) / dx[:,0]
    for i in xrange(1, ni-1):
        dvdx[:,i] = ((vv[t,k,:,i+1] - vv[t,k,:,i-1]) / dx[:,i]) / 2.
    dvdx[:,-1] = (vv[t,k,:,-1] - vv[t,k,:,-2]) / dx[:,-1]

    try:
        dudy[0,:]   = ((uu[t,k,1,:] * cosphi[1]) - (uu[t,k,0,:] * cosphi[0])) / dy[0,:]
        for j in xrange(1, nj-1):
            dudy[j,:] = (((uu[t,k,j+1,:] * cosphi[j+1]) - (uu[t,k,j-1,:] * cosphi[j-1])) / dy[j,:])/2.
        dudy[-1,:]  = ((uu[t,k,-2,:] * cosphi[-2]) - (uu[t,k,-1,:] * cosphi[-1])) / dy[-1,:]
    except:
        raise

    for i in xrange(ni):
        qq[:, i] = f0_p + dvdx[:,i] - dudy[:,i]/cosphi 

    return qq


def bin_coords(PV_array, timestep):
    '''
    bins coords based on if they are stratospheric or tropospheric relative to 
    if PV value is > 2

    writes the index coords in the form nk,ni,nj (comma seperated) in its respective file
    depending on if it's stratospheric or tropospheric
    '''
    global ni, nj, nk_mom, strato_file, tropo_file
    start_time = time.time()
    # reverses the levels so k will go from bottom to top (1000 - 0.1)
    PV_array = PV_array[::-1]

    # initializes lists
    strato_coords   = []
    tropo_coords    = []

    # if you wish to use a numpy solution (much faster), uncomment these two
    #strato_coords   = np.zeros((nk_mom,ni,nj))
    #tropo_coords    = np.zeros((nk_mom,ni,nj))

    # creates file
    strato  = open(strato_file,'a')
    tropo   = open(tropo_file,'a')
    strato.write("TIMESTEP,{0}\n".format(timestep))
    tropo.write("TIMESTEP,{0}\n".format(timestep))

    for i in xrange(ni):
        for j in xrange(nj):
            for k in xrange(nk_mom):
                if PV[k,i,j] > 2:
                    # if you wish to use a numpy solution (much faster), uncomment these two
                    #tropo_coords[:k,i,j] = 1
                    #strato_coords[k:,i,j] = 1
                    tropo_coords    += [(z, i, j) for z in xrange(k)]
                    strato_coords   += [(z, i, j) for z in xrange(k,nk_mom)]
                    break
    # conducts writing to file
    print 'WRITING BINNED INDEX POSITIONS'
    for coord in tropo_coords:
        tropo.write('{0},{1},{2}\n'.format(coord[0], coord[1], coord[2]))
    for coord in strato_coords:
        strato.write('{0},{1},{2}\n'.format(coord[0], coord[1], coord[2]))
    print 'DONE WRITING'
    strato.close()
    tropo.close()

    print "Binning time: {0}".format(str(time.time() - start_time))
    #return strato_coords, tropo_coords


def build_fst(params):
    '''
    (dict) -> (int, dict)
    builds the file as per the parameters defined in the params dict
    
    returns the file_id
    '''

    # makes an empty .fst file
    day = time.gmtime()
    temp = '' 
    for x in day: temp += str(x)
    #new_nc = '/home/ords/aq/alh002/pyscripts/workdir/pv_files/TEST5.fst'
    new_nc = '/home/ords/aq/alh002/pyscripts/workdir/pv_files/POTVOR_file_{0}.fst'.format(temp)
    tmp = open(new_nc, 'w+'); tmp.close()
    output_file = new_nc

    try:
        file_id = rmn.fnom(output_file)
        open_fst = rmn.fstouv(file_id, rmn.FST_RW)
        print(file_id, open_fst)

        MACC_grid = rmn.encodeGrid(params)
        print("Grids created.")
        print 'Grid Shape:' + str(MACC_grid['shape'])

        rmn.writeGrid(file_id, MACC_grid)
        toc_record = vgd.vgd_new_pres(const_pressure, ip1=MACC_grid['ig1'], ip2=MACC_grid['ig2'])
        vgd.vgd_write(toc_record, file_id)

        return file_id, MACC_grid

    except:
        rmn.fstfrm(file_id)
        rmn.fclos(file_id)
        raise


def get_grid_descriptors(file_id):
    '''
    gets grid descriptors from a fst file
    '''

    tic_record = rmn.FST_RDE_META_DEFAULT.copy()
    tac_record = rmn.FST_RDE_META_DEFAULT.copy()
    tac = rmn.fstinl(file_id, nomvar='>>')[0]
    tic = rmn.fstinl(file_id, nomvar='^^')[0]
    tic_record.update(rmn.fstprm(tic))
    tac_record.update(rmn.fstprm(tac))

    return tic_record, tac_record


def get_zonal_mean(array):
    '''
    array must be 3d of size nk, nlon, nlat
    '''
    
    zonal_mean = np.zeros(array.shape)
    for k in xrange(len(array)):
        for j in xrange(len(array[0,0])):
            zonal_mean[k,:,j] = array[k,:,j].mean()

    return zonal_mean


def shift_lon(array):
    '''
    shifts the longitude from -180 - 180 to 0 - 360
    '''
    new_array = np.zeros(array.shape)
    lonlen = len(new_array[0])
    new_array[:, :lonlen/2] = array[:, lonlen/2:]
    new_array[:, lonlen/2:] = array[:, :lonlen/2]
    return new_array


def vert_interp(pressure, org):
    '''
    org should be a 3d array, pressure should be a 3d array
    of pressures.
    '''
    global const_pressure    
    start_time = time.time()
    print "Conducting vertical interpolation..."

    lon         = len(org[0])
    lat         = len(org[0, 0])
    lev         = len(org)
    y_interp    = np.zeros((len(const_pressure), lon, lat))
    x_interp    = const_pressure

    for i in xrange(lon):
        for j in xrange(lat):
            try:
                x_initial = pressure[:, j, i]
                y_initial = org[:, i, j]
                y_interp[:, i, j] = np.interp(x_interp, x_initial, y_initial)
            except:
                print pressure.shape
                raise
                
    print "That took " + str(time.time() - start_time) + " seconds"
    return y_interp


def splice_month(open_var, m_int):
    '''
    takes an open var, splices out the timesteps that dont regard
    to that month (with m_int)
    '''

    if len(open_var.time) == 1464:
        leap_year = True
    elif len(open_var.time) == 1460:
        leap_year = False
    else:
        print("File must contain all 4 time values for each day of the year")
        exit()
    
    if leap_year:
        stored_array = [
                        open_var[:124], open_var[124:240], open_var[240:364], open_var[364:484],
                        open_var[484:608], open_var[608:728], open_var[728:852], open_var[852:976],
                        open_var[976:1096], open_var[1096:1220], open_var[1220:1340], open_var[1340:]
                        ]
    else:
        stored_array = [
                        open_var[:124], open_var[124:236], open_var[236:360], open_var[360:480],
                        open_var[480:604], open_var[604:724], open_var[724:848], open_var[848:972],
                        open_var[972:1092], open_var[1092:1216], open_var[1216:1336], open_var[1336:]
                        ]
    return stored_array[m_int]


def build_hhmmss(timestep):
    '''
    calculates seconds from timestep

    given seconds (0 being 00 hours 00 minutes 00 seconds in the day, 
    86399 being 23 hours 59 minutes 59 seconds in the day), returns an
    int in the form hhmmss00 for use in record definition with rmn.newdate 
    '''
    seconds = int((((timestep) % 4.)/4) * 86400.)

    hh   = int(math.floor(seconds/3600))
    mm   = int(math.floor((seconds - (hh * 3600)) / 60))
    ss   = int(math.floor((seconds - (hh * 3600) - (mm * 60)) / 60))

    time_int = int(hh * 1e6 + mm * 1e4 + ss * 1e2)
    
    print 'TIMEINT: {0}'.format(str(time_int))
    return time_int

def get_pressures(open_nc, m_int):
    '''
    (open .nc, int) -> np.ndarray
    gets pressures from open_file, based on m_int as an integer
    representing the month you wish to obtain.
    '''
    global Ak, Bk, const_pressure
                                 
    start_time = time.time()
    print "Getting pressures..."
    
    lnsp = open_nc.lnsp  # open file containing every lnsp for that year
    lnsp_array = splice_month(lnsp, m_int)
    print lnsp_array.shape

    shape = list(lnsp_array.shape); shape.insert(1, 60)  # creates a dimension for levels
    print "PRESSURE SHAPE IS {0}".format(str(shape))
    pressure = np.zeros(shape)
    full_pressure = np.zeros(shape)
    try:
        for lev in xrange(len(pressure[0])):
            pressure[:, lev] = ((Ak[lev] + Bk[lev]*np.exp(lnsp_array)) / 100.)

    except:
        debug_tuple = (lev,)
        print "LEV: {0}".format(debug_tuple)
        print pressure.shape, lnsp_array.shape
        raise

    for lev in xrange(1, len(pressure[0])):
        try:
            full_pressure[:, lev-1] = (pressure[:,lev] + pressure[:, lev-1]) / 2.
        except IndexError:
            if k+1 != 0:
                raise
            else:
                pass
    full_pressure[:,-1] = 1000.0

#    print pressure[5, :, 12, 33]
#    print pressure[15, :, 12, 23]
#    print pressure[45, :, 12, 53]
#    print pressure[75, :, 2, 39]
#    print full_pressure[23, :, 12, 23]
    print "That took " + str(time.time() - start_time) + " seconds"
    return full_pressure 


##### MAIN #####
# this portion of the code handles parsing values from the nc file
for filename in filenames:
    pv_list  = []
    zon_list = []
    nc = pyg.open(filename)
    lnsp_file = pyg.open(lnsp_files[0])

    # all instances of '[:4]' are to limit to 4 timesteps, or 1 day (in this case 01012012)
    uu = nc.u.get()[:28,:,::-1]
    vv = nc.v.get()[:28,:,::-1]
    qq = nc.vo.get()[:28,:,::-1]
    th = nc.t.get()[:28,:,::-1]
    pressures = get_pressures(lnsp_file, 0)
    pressures = pressures[:28,:,:len(uu[0,0])]
    pressures = pressures[:,:,::-1]

    uu.setflags(write=True)
    vv.setflags(write=True)
    qq.setflags(write=True)
    th.setflags(write=True)

    # for some reason the nc file gets th[t,k,-1,:] as all the same temperature.
    # in other words, at longitude index (West 180 degrees) 0, 
    # all the temperatures are recorded as the same in the ncdf
    # i will conduct mean of th[t,k,-2,:] and th[t,k,0,:] to estimate these values
    #th[:,:,-1,:] = (th[:,:,-2,:] + th[:,:,-3,:]) / 2. 

    lon = nc.longitude.values
    lat = nc.latitude.values[:len(uu[0,0])]
    lat = lat[::-1]
    ni = len(lon)
    nj = len(lat)
    nk_mom = len(pressures[0])#len(const_pressure)
    nk_thm = nk_mom
    qq2 = np.zeros((nk_mom, nj, ni))

    # stuff on level k - 1
    uulast = np.zeros([nj, ni])
    vvlast = np.zeros([nj, ni])

    # stuff on level k + 0.5 (th will be on this level)
    dthdx = np.zeros([nj, ni])
    dthdy = np.zeros([nj, ni])

    # stuff on level k - 0.5
    thlast = np.zeros([nj, ni])

    # calculation must be done with actual lon-lat, not rotated
    cosphi = np.cos(lat * (np.pi / 180.))
    cosphi[-1] = cosphi[-2]
    f0_p   = 2 * np.pi/86400. * np.sin(lat * np.pi/180.)
    for t in xrange(len(qq)): 
        for x in xrange(ni):
            for k in xrange(nk_mom):
                qq[t,k,:,x] = f0_p + qq[t,k,:,x]
    dx     = np.zeros([nj, ni])  # shape is lon, lat
    dy     = np.zeros([nj, ni])

    # setting dx and dy partial derivative arrays through leapfrogging
    try:
        print "Setting dx..."
        for i in xrange(1, ni-1):
            # in absolute vorticity, this is part of du*cosphi / dx
            dx[:,i] = a * cosphi * ((lon[i+1]-lon[i-1]) * np.pi/360.)
            # bottom level init. 
        dx[:,0] = a * cosphi * (lon[1] - lon[0]) * np.pi/180. 
        dx[:,-1] = a * cosphi * (lon[-1] - lon[-2]) * np.pi/180.  # top level init
    except:
        print i 
        raise
    try:
        print "Setting dy..."
        for j in xrange(1, nj-1):
            dy[j,:] = a * (lat[j+1] - lat[j-1]) * np.pi/360.
        dy[0,:] = a * (lat[1] - lat[0]) * np.pi/180.  # bottom level init. NOTE: why do we not divide by 2 here?
        dy[-1,:] = a * (lat[-1] - lat[-2]) * np.pi/180.  # top level init
    except:
        print j
        raise

    # params for partial grid descriptors and grid definition for fst files
    params0 = {
            'grtyp' : 'Z',
            'grref' : 'L',
            'nj'    : nj,
            'ni'    : ni,
            'lat0'  : 9,#lat[0],
            'lon0'  : -180,#lon[0],
            'dlat'  : dlatlon,
            'dlon'  : dlatlon
            }

    tempq = np.zeros((nk_mom, nj, ni))
    tempv = np.zeros((nj, ni))
    tempu = np.zeros((nj, ni))

    # builds the fst file as per the params in params0
    file_id, MACC_grid = build_fst(params0)

    # gets the grid descriptors for later use with definition of data records
    tic_record, tac_record = get_grid_descriptors(file_id)
    
    # creates surf/pressure array and turns temp into theta
    for o in xrange(len(uu)):
        for k in xrange(len(uu[0])):
            thet1 = pressures[o, -1] / pressures[o, k]
            th[o,k] = th[o,k] * thet1 ** kappa

    # holds all the bin lists for each timestep
    #strato_timed_bin_list   = np.zeros((len(uu), nk_mom, ni, nj))
    #tropo_timed_bin_list    = np.zeros((len(uu), nk_mom, ni, nj))
    start_time = time.time()
    for t in xrange(len(uu)):
    #for t in xrange(1):
        # Potential vorticity field
        PV          = np.zeros([nk_mom, nj, ni])
        pres_levs   = np.zeros(PV.shape)  # really lazy
        dp          = np.zeros(PV.shape)  # really lazy

        pres_levs   = pressures[t]
        surf_pres   = pres_levs[-1]
        tempq       = qq[t]

        try:
            dp[-1] = 100 * (pres_levs[-1] - pres_levs[-2])
            dp[0] = 100 * pres_levs[0]
            for k in xrange(1, nk_mom-1):
                dp[k] = 100 * (pres_levs[k] - pres_levs[k-1])
        except:
            print k
            raise

        # calculation of term1/term2 (level by level)
        # omits the final level (i = nk_mom)
        for k in xrange(nk_mom):
            print "Record Level: " + str(k+1)
            if k == 0:
                term1   = np.zeros((nk_mom, nj, ni))
                term2   = np.zeros((nk_mom, nj, ni))
                dt      = np.zeros((nk_mom, nj, ni))

            elif k != nk_mom - 1:
                # define the dthdy and dthdx based of thermo level
                dthdy[0,:] = (th[t,k,1,:] - th[t,k,0,:]) / dy[0,:]  # first dthdy level
                for j in xrange(1, nj-1):
                    dthdy[j,:] = ((th[t,k,j+1,:] - th[t,k,j-1,:]) / dy[j,:]) / 2.
                dthdy[-1,:] = (th[t,k,-1,:] - th[t,k,-2,:]) / dy[-1,:]  # last dthdy level
                dthdx[:,0] = (th[t,k,:,1] - th[t,k,:,0]) / dx[:,0]  # first dthdy level
                for i in xrange(1, ni-1):
                    dthdx[:,i] = ((th[t,k,:,i+1] - th[t,k,:,i-1])/ dx[:,i]) / 2.
                dthdx[:,-1] = (th[t,k,:,-1] - th[t,k,:,-2]) / dx[:,-1] # last dthdy level
                
            else:
                print "FINAL LEVEL REACHED"
                # final level (k = nk)
                dthdy[0,:]  = (th[t,k,1,:] - th[t,k,0,:]) / dy[0,:]
                dthdy[-1,:] = (th[t,k,-1,:] - th[t,k,-2,:]) / dy[-1,:]
                for j in xrange(1, nj-1):
                    dthdy[j,:] = ((th[t,k,j+1,:] - th[t,k,j-1,:]) / dy[j,:])/2

                dthdx[:,0]  = (th[t,k,:,1] - th[t,k,:,0]) / dx[:,0]
                dthdx[:,-1] = (th[t,k,:,-1] - th[t,k,:,-2]) / dx[:,-1]
                for i in xrange(1, ni-1):
                    dthdx[:,i] = ((th[t,k,:,i+1] - th[t,k,:,i-1]) / dx[:,i])/2

            '''
            from this point onwards, horizontal grid leapfrogging is complete.
            the next steps performed are going to be vertical integration of the terms
            onto pressure levels. these pressure levels will be integrated at "half" levels,
            not full levels. in the horizontal grid, leapfrogging was done with 
            level k-1 and level k+1 to give level k. this will be done on level k and k+1 to
            give level k + 0.5. this level doesnt really exist, it's just the inbetween of the
            two levels. 
            '''
            
            # PV values are on fields k-1/2 respective to the values obtained
            if k != 0:
                term1[k]    = g0 * ((vv[t,k] - vv[t,k-1])/dp[k]) * dthdx
                term2[k]    = g0 * ((uu[t,k] - uu[t,k-1])/dp[k]) * dthdy
                dt[k]       = th[t,k] - th[t,k-1]
        
        # these exist on half levels -0.5 of the original level
        PV      = ((-g0 * (dt/dp) * tempq) + term1 - term2) * 1e6

        PV  = np.transpose(PV, (0,2,1))
        PV  = shift_lon(PV)
        bin_coords(PV,t)
        PV  = vert_interp(pressures[t], PV)
        zonal_mean = get_zonal_mean(PV)  # change PV depending on which you choose to use

        pv_list.append(PV)
        zon_list.append(zonal_mean)
        
    # BOOKMARK: FSTFILE #
    
    ## this portion of the code simply opens a new file and ports PV onto it
    print "Total time taken: {0}".format(str(time.time() - start_time))
    #exit()
    l += 1
    for t, array in enumerate(pv_list):
        try:
            # copies the default record
            new_record = rmn.FST_RDE_META_DEFAULT.copy() 
            time_int = build_hhmmss(t)
            for rp1 in xrange(len(const_pressure)):  # writes a record for every level (as a different ip1)
                # converts rp1 into a ip1 with pressure kind
                ip1 = rmn.convertIp(rmn.CONVIP_ENCODE, const_pressure[rp1], rmn.KIND_PRESSURE)
                new_record.update(MACC_grid)
                new_record.update({  # Update with specific meta
                    'nomvar': 'PV',
                    'typvar': 'C', 
                    'etiket': 'MACCRean',
                    'ip2'   : t,
                    'ni'    : MACC_grid['ni'],
                    'nj'    : MACC_grid['nj'],
                    'ig1'   : tic_record['ip1'],
                    'ig2'   : tic_record['ip2'],
                    'ig3'   : tic_record['ip3'],
                    'ig4'   : tic_record['ig4'],
                    'dateo' : rmn.newdate(rmn.NEWDATE_PRINT2STAMP,20120100 + int(math.floor(t/4+1)), time_int),
                    'deet'  : int(86400/4),  # timestep in secs
                    'ip1'   : ip1
                    })
                
                tmp = array[rp1]
                tmp = np.asfortranarray(tmp)
                # data array is structured as tmp = monthly_mean[level] where monthly_mean is [level, lat, lon]
                # Updates with data array in the form (lon x lat)
                new_record.update({'d': tmp.astype(np.float32)}) 
                print "Defined a new record with dimensions ({0}, {1})".format(new_record['ni'], new_record['nj'])
                rmn.fstecr(file_id, new_record)  # write the dictionary record to the file as a new record
                
                zonal = zon_list[t]
                tmp = zonal[rp1]
                tmp = np.asfortranarray(tmp)
                zonal_record = new_record
                zonal_record.update({
                    'nomvar': 'ZO',
                    'd'     : tmp.astype(np.float32)
                    })
                print "Defined a zonal mean record with dimensions ({0}, {1})".format(new_record['ni'], new_record['nj'])       
                
                rmn.fstecr(file_id, zonal_record)  # write the dictionary record to the file as a new record
        except:
            rmn.fstfrm(file_id)
            rmn.fclos(file_id)
            raise

    rmn.fstfrm(file_id)
    rmn.fclos(file_id)


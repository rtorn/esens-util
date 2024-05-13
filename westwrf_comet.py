import os
import time
import sys
import cfgrib
import gzip
import urllib
import json
import datetime as dt
import numpy as np
import xarray as xr

def stage_grib_files(datea, config):
    '''
    This is a generic class for copying or linking grib file data from a specific location
    to a directory where calculations are performed.  No matter where these calculations are
    carried out, this routine must exist.
 
    This particular instance is for GEFS data on the NOAA jet machine.  In this 
    case, the individual ensemble member files are combined into one file per forecast
    time in the work directory, which is more efficient for reading.

    Attributes:
        datea (string):  The initialization time of the forecast (yyyymmddhh)
        config  (dict):  The dictionary with configuration information
    '''

    freq = config['model'].get('fcst_hour_int', 12)
    fmax = config['model'].get('fcst_hour_max', 120)

    init    = dt.datetime.strptime(datea, '%Y%m%d%H')
    init_s  = init.strftime("%y%j%H")

    #  Make the work directory if it does not exist
    if not os.path.isdir(config['locations']['work_dir']):
       try:
          os.makedirs(config['locations']['work_dir'])
       except OSError as e:
          raise e

def stage_atcf_files(datea, bbnnyyyy, config):
    '''
    This is a generic routine for copying or linking ATCF file data from a specific location
    to a directory where calculations are performed.  No matter where these calculations are
    carried out, this routine must exist.

    The result is a set of ATCF files in the work directory of the format atcf_NN.dat, 
    where NN is the ensemble member number.

    This particular instance is for GEFS data on the NOAA machine jet.  In this 
    case, all ATCF data for a particular storm is in one file, so the code waits for this
    initialization time to exist, then uses sed to get the lines attributed to each 
    ensemble member and places that data in a seperate file.

    Attributes:
        datea    (string):  The initialization time of the forecast (yyyymmddhh)
        bbnnyyyy (string):  TC Identification, where bb is the basin, nn is the number, yyyy is year
        config     (dict):  The dictionary with configuration information
    '''

    src  = '{0}/{1}.a{2}.dat'.format(config['locations']['atcf_dir'],datea,bbnnyyyy)

    #  Wait for the ensemble ATCF information to be placed in the file
    while ( len(os.popen('sed -ne /AP/p {0}/{1}.a{2}.dat'.format(config['locations']['atcf_dir'],datea,bbnnyyyy)).read()) == 0 ):
       time.sleep(20.7)

    #  Wait for the file to be finished being copied
    while ( (time.time() - os.path.getmtime(src)) < 60 ):
       time.sleep(10)

    for n in range(int(config['model']['num_ens']) + 1):

       if ( n > 0 ):
          modid = 'AP'
       else:
          modid = 'AC'

       nn = '%0.2i' % n
       file_name = '{0}/atcf_{1}.dat'.format(config['locations']['work_dir'],nn)

       #  If the specific member's ATCF file does not exist, copy from the source file with sed.
       if not os.path.isfile(file_name):

          fo = open(file_name,"w")
          fo.write(os.popen('sed -ne /{0}/p {1}/{0}.a{2}.dat | sed -ne /{3}{4}/p'.format(datea,config['locations']['atcf_dir'],bbnnyyyy,modid,nn)).read())
          fo.close()


def stage_best_file(bbnnyyyy, config):
    '''
    This is a generic routine for copying or linking ATCF best track file data from a specific location
    to a directory where calculations are performed.  No matter where these calculations are
    carried out, this routine must exist.

    The result is a single best track file in the work directory of the format bbnnyyyy.dat, 
    where bb is the basin, nn is the TC number.

    This particular instance copies data from the NHC data server.  

    Attributes:
        bbnnyyyy (string):  TC Identification, where bb is the basin, nn is the number, yyyy is year
        config     (dict):  The dictionary with configuration information
    '''

    try:    #  Look for the data in the real-time directory

      filei = urllib.request.urlopen('{0}/b{1}.dat'.format(config['locations'].get('best_dir','https://ftp.nhc.noaa.gov/atcf/btk'),bbnnyyyy))
      fileo = open('{0}/b{1}.dat'.format(config['locations']['work_dir'],bbnnyyyy), 'wb')
      fileo.write(filei.read())
      filei.close()
      fileo.close()

    except:    #  If the file is not present in the real-time directory, look in the archive.

      src  = '{0}/{1}/b{2}.dat.gz'.format(config['locations'].get('best_dir_alt','https://ftp.nhc.noaa.gov/atcf/archive'),bbnnyyyy[4:8],bbnnyyyy)

      #  Unzip the file from the NHC server, write the file to the work directory
      gzfile = gzip.GzipFile(fileobj=urllib.request.urlopen(src))
      uzfile = open('{0}/b{1}.dat'.format(config['locations']['work_dir'],bbnnyyyy), 'wb')
      uzfile.write(gzfile.read())
      gzfile.close()
      uzfile.close()


class ReadGribFiles:
    '''
    This is a generic class that is used to read specific forecast fields from a grib file at a specific 
    forecast hour.  This class is a generic class, but differs depending on the grib file format and fields.
    The init routine creates the class and constructs a field dictionary.

    This particular instance is for GEFS data on the NOAA jet machine.  Each GEFS file contains all
    members for a specific forecast time.  

    Attributes:
        datea (string):  The initialization time of the forecast (yyyymmddhh)
        fhr      (int):  Forecast hour
        config  (dict):  The dictionary with configuration information
    '''

    def __init__(self, datea, fhr, config):

       init    = dt.datetime.strptime(datea, '%Y%m%d%H')
       init_s  = init.strftime("%y%j%H")

       fdate   = init + dt.timedelta(hours=fhr)
       fdatestr = fdate.strftime("%Y-%m-%d_%H_%M_%S")

       self.datea = datea
       self.fdatestr = fdatestr
       self.config = config

       self.filetype = config['model'].get('file_format', 'westwrf')

       if 'ensemble_list' in self.config:

          elist = [e.strip() for e in config['model']['ensemble_list'].split(',')]         
          if len(elist) > 1:
             self.enslist = elist
          else:
             self.enslist = []
             print(elist[0])
             with open(elist[0].replace('{yyyymmddhh}',self.datea)) as f:
                for line in f:
                   self.enslist.append(line.rstrip())

       else:

          self.enslist = ['ecm000', 'ecm001', 'ecm002', 'ecm003', 'ecm004', 'ecm005', 'ecm006', 'ecm007', 'ecm008', 
                          'ecm009', 'ecm010', 'ecm011', 'ecm012', 'ecm013', 'ecm014', 'ecm015', 'ecm016', 'ecm017', 
                          'ecm018', 'ecm019', 'ecm020', 'ecm021', 'ecm022', 'ecm023', 'ecm024', 'ecm025', 
                          'gefs001', 'gefs002', 'gefs003', 'gefs004', 'gefs005', 'gefs006', 'gefs007', 'gefs008', 'gefs009', 
                          'gefs010', 'gefs011', 'gefs012', 'gefs013', 'gefs014', 'gefs015', 'gefs016', 'gefs017', 'gefs018', 
                          'gefs019', 'gefs020', 'gefs021', 'gefs022', 'gefs023', 'gefs024', 'gefs025']

       self.nens = len(self.enslist)

       file_name = self.member_name(0)

       try:
          self.domdict = xr.open_dataset(file_name, cache=False)
       except IOError as exc:
          raise RuntimeError('Failed to open {0}'.format(file_name)) from exc

       if self.filetype == 'minghua':
          self.horlist = {'lat_new': 'lat', 'lon_new': 'lon'}
          self.dimlist = {'lon': 'longitude', 'lat': 'latitude', 'pres_lev': 'isobaricInhPa', 'pres_new': 'isobaricInhPa'}
       elif self.filetype == 'westwrf':
          self.horlist = {'south_north': 'lat', 'west_east': 'lon'}
          self.dimlist = {'lon': 'longitude', 'lat': 'latitude', 'pressure': 'isobaricInhPa'}

       self.domdict = self.domdict.rename(self.dimlist)

       if np.max(self.domdict.coords['longitude']) > 180:
          self.domdict.coords['longitude'] = (self.domdict.coords['longitude'] + 180) % 360 - 180

       if config['model'].get('flip_lon','False') == 'True':
          self.domdict.coords['longitude'] = (self.domdict.coords['longitude'] + 360.) % 360.

       self.grid_type = config['model'].get('grid_type','Cassini')
       self.truelat1  = float(config['model'].get('grid_truelat1', '39.0'))
       self.truelat2  = float(config['model'].get('grid_truelat2', '50.0'))
       self.stdlon    = float(config['model'].get('grid_stdlon', '-125.0'))
       self.dx        = float(config['model'].get('grid_dx', '9000.0'))
       self.latinc    = float(config['model'].get('grid_latinc', '0.08'))
       self.loninc    = float(config['model'].get('grid_loninc', '0.08'))
       self.lat0      = float(config['model'].get('grid_lat0', '51.0'))
       self.lon0      = float(config['model'].get('grid_lon0', '180.0'))
       self.knowni    = float(config['model'].get('knowni', '0.0'))
       self.knownj    = float(config['model'].get('knownj', '0.0'))
       self.lat1      = self.domdict.latitude[0,0].values
       self.lon1      = self.domdict.longitude[0,0].values

       if self.grid_type == 'Cassini':
          self.set_cassini_domain()  
       elif self.grid_type == 'LambertConformal':
          self.set_lc_domain()

       #  This is a dictionary that maps from generic variable names to the name of variable in file
       if self.filetype == 'minghua':

          self.var_dict = {'zonal_wind': 'U',             \
                           'meridional_wind': 'V',        \
                           'zonal_wind_10m': 'U10',       \
                           'meridional_wind_10m': 'V10',  \
                           'geopotential_height': 'GEOPT',   \
                           'temperature': 'T',            \
                           'specific_humidity': 'Q',      \
                           'relative_humidity': 'RH',      \
                           'sea_level_pressure': 'PMSL', \
                           'precipitation': 'RAIN', \
                           'ivt': 'IVT' }

       elif self.filetype == 'westwrf':  

          self.var_dict = {'zonal_wind': 'u_tr_p',             \
                           'meridional_wind': 'v_tr_p',        \
                           'zonal_wind_10m': 'u_10m_tr',       \
                           'meridional_wind_10m': 'v_10m_tr',  \
                           'geopotential_height': 'Z_p',   \
                           'temperature': 'T_p',            \
                           'specific_humidity': 'q_p',      \
                           'relative_humidity': 'RH',      \
                           'sea_level_pressure': 'slp', \
                           'precipitation': 'precip', \
                           'ivt': 'IVT' }

       if self.var_dict['specific_humidity'] in self.domdict:
          self.has_specific_humidity = True
       else:
          self.has_specific_humidity = False

       self.has_total_precip = True


    def member_name(self, member):

       if self.filetype == 'minghua':
          file_name = "{0}/{1}/{2}/wrfout_d01_{3}_subset_variables.nc".format(self.config['locations']['model_dir'],self.datea,self.enslist[member],self.fdatestr)
       elif self.filetype == 'westwrf':
          file_name = "{0}/{2}/wrfcf_d01_{3}.nc".format(self.config['model_dir']['locations'].replace("{yyyymmddhh}",self.datea),self.datea,self.enslist[member],self.fdatestr)
#          file_name = "{0}/{1}/ecm{2}/wrfcf_d01_{3}.nc".format(self.config['model_dir'],self.datea,fmem,self.fdatestr)

       return file_name


    def set_lc_domain(self):

       if self.truelat1 < 0.0:
          self.hemi = -1.0
       else:
          self.hemi = 1.0
       self.rebydx = 6370000. / self.dx

       # Compute cone factor
       if abs(self.truelat1-self.truelat2) > 0.10:
          self.cone = np.log10(np.cos(np.radians(self.truelat1))) - np.log10(np.cos(np.radians(self.truelat2)))
          self.cone = self.cone / (np.log10(np.tan(np.radians(45.0 - abs(self.truelat1)/2.0))) - \
                      np.log10(np.tan(np.radians(45.0 - abs(self.truelat2)/2.0))))
       else:
          self.cone = np.sin(np.radians(abs(self.truelat1)))

       # Compute longitude differences and ensure we stay out of the forbidden "cut zone"
       deltalon1 = self.lon1 - self.stdlon
       if deltalon1 > 180.0:
          deltalon1 = deltalon1 - 360.0
       if deltalon1 < -180.0:
          deltalon1 = deltalon1 + 360.0

       # Convert truelat1 to radian and compute COS for later use
       ctl1r = np.cos(np.radians(self.truelat1))

       # Compute the radius to our known lower-left (SW) corner
       self.rsw = self.rebydx * ctl1r/ self.cone * \
              (np.tan(np.radians(90.0*self.hemi-self.lat1)/2.0) / \
               np.tan(np.radians(90.0*self.hemi-self.truelat1)/2.0))**self.cone

       # Find pole point
       arg = self.cone*np.radians(deltalon1)
       self.polei = 1.0 - self.hemi*self.knowni - self.hemi * self.rsw * np.sin(arg)
       self.polej = 1.0 + self.hemi*self.knownj + self.rsw * np.cos(arg)


    def set_cassini_domain(self):

       self.hemi = 1.0

       # Try to determine whether this domain has global coverage
       if np.abs(self.latinc/2.0 + 90.0) < 0.001 and \
          np.abs(np.mod(self.lon1 - self.loninc/2.0 - self.stdlon,360.0)) < 0.001:
          global_domain = True
       else:
          global_domain = False

       if np.abs(self.lat0) != 90.0 and (not global_domain):
          self.lat1, self.lon1 = self.rotate_coords(self.lat1, self.lon1, -1)


    def rotate_coords(self, ilat, ilon, direction):

       # Convert all angles to radians
       phi_np = np.radians(self.lat0)
       lam_np = np.radians(self.lon0)
       lam_0  = np.radians(self.stdlon)
       rlat   = np.radians(ilat)
       rlon   = np.radians(ilon)

       if direction < 0:
          dlam = np.pi - lam_0
       else:
          dlam = lam_np

       sinphi = np.cos(phi_np)*np.cos(rlat)*np.cos(rlon-dlam) + np.sin(phi_np)*np.sin(rlat)
       cosphi = np.sqrt(1.0-sinphi*sinphi)
       coslam = np.sin(phi_np)*np.cos(rlat)*np.cos(rlon-dlam) - np.cos(phi_np)*np.sin(rlat)
       sinlam = np.cos(rlat)*np.sin(rlon-dlam)

       if cosphi != 0.0:
          coslam = coslam/cosphi
          sinlam = sinlam/cosphi
       olat = np.degrees(np.arcsin(sinphi))
       olon = np.degrees(np.arctan2(sinlam,coslam)-dlam-lam_0+lam_np)
       if olon < -180.0:
          olon = olon + 360.0
       if olon > 180.0:
          olon = olon - 360.0

       return olat, olon


    def latlon_to_ij(self, lat, lon):

       if self.grid_type == 'Cassini':

          if np.abs(self.lat0) != 90.0:
             comp_lat, comp_lon = self.rotate_coords(lat, lon, -1)
          else:
             comp_lat = lat
             comp_lon = lon

          deltalat = comp_lat - self.lat1
          deltalon = comp_lon - self.lon1

          if deltalon <   0.0: 
             deltalon = deltalon + 360.0
          if deltalon > 360.0:
             deltalon = deltalon - 360.0

          i = deltalon / self.loninc
          j = deltalat / self.latinc

          if i <= 0.0:
             i = i + 360.0 / self.loninc
          if i > 360.0/self.loninc:
             i = i - 360.0 / self.loninc

          return int(i + self.knowni), int(j + self.knownj)

       elif self.grid_type == 'LambertConformal':

          deltalon = lon - self.stdlon
          if deltalon > 180.0:
             deltalon = deltalon - 360.0
          if deltalon < -180.0:
             deltalon = deltalon + 360.0

          ctl1r = np.cos(np.radians(self.truelat1))

          #  Radius to desired point
          rm = self.rebydx * ctl1r/self.cone * \
              (np.tan(np.radians(90.0*self.hemi-lat)/2.0) / \
               np.tan(np.radians(90.0*self.hemi-self.truelat1)/2.0))**self.cone

          arg = self.cone*np.radians(deltalon)
          i = int(self.hemi * (self.polei + self.hemi * rm * np.sin(arg)))
          j = int(self.hemi * (self.polej - rm * np.cos(arg)))

          return i, j


    def set_var_bounds(self, varname, vdict):
       '''
       This is a generic routine that is used to determine the appropriate starting and ending latitude,
       longitude, and if appropriate pressure levels to extract from the grib files.  This routine takes into
       account the order of the coordinate variables.  If the bounds are not specified, the routine uses all
       coordinate indicies.  The resulting bounds are added back into the dictionary object and returned for
       use in other routines within the class.

       Attributes:
           varname (string):  The name of the variable that will be extracted from file
           vdict     (dict):  The dictionary object with variable information
       '''

       if 'latitude' in vdict and 'longitude' in vdict:

          vdict['i_start'], vdict['j_start'] = self.latlon_to_ij(vdict['latitude'][0], vdict['longitude'][0])
          vdict['i_end'],   vdict['j_end']   = self.latlon_to_ij(vdict['latitude'][1], vdict['longitude'][1])
          vdict['i_start'] = vdict['i_start'] - 1
          vdict['j_start'] = vdict['j_start'] - 1

       else:

          vdict['i_start'] = 0
          vdict['j_start'] = 0
          vshape = self.domdict.get(self.var_dict[varname]).shape
          vdict['i_end'] = vshape[-1]
          vdict['j_end'] = vshape[-2]

       if 'isobaricInhPa' in vdict:

          #  See of pressure level values are reversed
          if self.domdict.isobaricInhPa[0] > self.domdict.isobaricInhPa[1]:
             vdict['pres_start'] = int(np.where(self.domdict.isobaricInhPa[:]==int(vdict['isobaricInhPa'][1]))[0])
             vdict['pres_end']   = int(np.where(self.domdict.isobaricInhPa[:]==int(vdict['isobaricInhPa'][0]))[0])+1
          else:
             vdict['pres_start'] = int(np.where(self.domdict.isobaricInhPa[:]==int(vdict['isobaricInhPa'][0]))[0])
             vdict['pres_end']   = int(np.where(self.domdict.isobaricInhPa[:]==int(vdict['isobaricInhPa'][1]))[0])+1

       return vdict


    def create_ens_array(self, varname, nens, vdict):
       '''
       This is a generic routine that is used to create an xarray object that contains all ensemble
       members for a particular field, with the proper units and coordinate variables.  The resulting
       array has dimensions (ensemble, latitude, longitude).  The routine
       takes into account the order of the lat/lon arrays

       Attributes:
           varname (string):  The name of the variable that will be extracted from file
           nens       (int):  number of ensemble members
           vdict     (dict):  The dictionary object with variable information
       '''

       #  Create attributes based on what is in the file
       attrlist = {}
       if 'description' in vdict:
         attrlist['description'] = vdict['description']
       if 'units' in vdict:
         attrlist['units'] = vdict['units']
       if '_FillValue' in vdict:
         attrlist['_FillValue'] = vdict['_FillValue']

       #  Create an empty xarray that can be used to copy data into
       lon = self.domdict.longitude[vdict['j_start']:vdict['j_end'], vdict['i_start']:vdict['i_end']].rename(self.horlist)
       lat = self.domdict.latitude[vdict['j_start']:vdict['j_end'], vdict['i_start']:vdict['i_end']].rename(self.horlist)


       ensarr = xr.DataArray(name='ensemble_data', data=np.zeros([nens, len(lat[:,0]), len(lat[0,:])]),
                             dims=['ensemble', 'lat', 'lon'], attrs=attrlist, 
                             coords={'ensemble': [i for i in range(nens)], 'latitude': lat, 'longitude': lon}) 

       return(ensarr)


    def read_pressure_levels(self, varname):

       nn = '%0.3i' % 1
       return self.domdict.isobaricInhPa[:]


    #  Function to read a single ensemble member's forecast field
    def read_grib_field(self, varname, member, vdict):
       '''
       This is a generic routine that is used to read a forecast field from a single ensemble member
       based on the information contained within the dictionary vdict.

       Attributes:
           varname (string):  Name of the variable that will be extracted from file
           member     (int):  ensemble member to read from file
           vdict     (dict):  The dictionary object with variable information
       '''

       ds = xr.open_dataset(self.member_name(member), cache=False).rename(self.dimlist)

       #  Read a single pressure level of data, if this is a variable that has pressure levels
       if 'isobaricInhPa' in vdict:

          vout = ds.get(self.var_dict[varname])[0,vdict['pres_start']:vdict['pres_end'], \
                                                vdict['j_start']:vdict['j_end'],vdict['i_start']:vdict['i_end']].squeeze().rename(self.horlist)

       #  Read the only level if it is a single level variable
       else:

          vout = ds.get(self.var_dict[varname])[0,vdict['j_start']:vdict['j_end'],vdict['i_start']:vdict['i_end']].squeeze().rename(self.horlist)

       ds.close()
       del ds

       return(vout)


    def close_files(self):

       self.domdict.close()
       del self.domdict


if __name__ == '__main__':

    src1 = "/Users/parthpatwari/RA_Atmospheric_Science/Old_Code/atcf_data"
    grib_src = "/Users/parthpatwari/RA_Atmospheric_Science/GRIB_files"
    dest1 = "/Users/parthpatwari/RA_Atmospheric_Science/New_Code/atcf_data"
    atcf_src = "/Users/parthpatwari/RA_Atmospheric_Science/Old_Code/atcf_data"
    # c1 = CopyFiles(src, dest)
    # if c1.checkandcreatedir():
    #     c1.copy_filestowork()
    g1 = ReadGribFiles(grib_src, '2019082900', 180)
    a1 = Readatcfdata(atcf_src)

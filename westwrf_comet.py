import os
import time
import sys
import cfgrib
import gzip
import urllib
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

    freq = config.get('fcst_hour_int', 12)
    fmax = config.get('fcst_hour_max', 120)

    init    = dt.datetime.strptime(datea, '%Y%m%d%H')
    init_s  = init.strftime("%y%j%H")

    #  Make the work directory if it does not exist
    if not os.path.isdir(config['work_dir']):
       try:
          os.makedirs(config['work_dir'])
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

    src  = '{0}/{1}.a{2}.dat'.format(config['atcf_dir'],datea,bbnnyyyy)

    #  Wait for the ensemble ATCF information to be placed in the file
    while ( len(os.popen('sed -ne /AP/p {0}/{1}.a{2}.dat'.format(config['atcf_dir'],datea,bbnnyyyy)).read()) == 0 ):
       time.sleep(20.7)

    #  Wait for the file to be finished being copied
    while ( (time.time() - os.path.getmtime(src)) < 60 ):
       time.sleep(10)

    for n in range(int(config['num_ens']) + 1):

       if ( n > 0 ):
          modid = 'AP'
       else:
          modid = 'AC'

       nn = '%0.2i' % n
       file_name = '{0}/atcf_{1}.dat'.format(config['work_dir'],nn)

       #  If the specific member's ATCF file does not exist, copy from the source file with sed.
       if not os.path.isfile(file_name):

          fo = open(file_name,"w")
          fo.write(os.popen('sed -ne /{0}/p {1}/{0}.a{2}.dat | sed -ne /{3}{4}/p'.format(datea,config['atcf_dir'],bbnnyyyy,modid,nn)).read())
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

      filei = urllib.request.urlopen('{0}/b{1}.dat'.format(config.get('best_dir','https://ftp.nhc.noaa.gov/atcf/btk'),bbnnyyyy))
      fileo = open('{0}/b{1}.dat'.format(config['work_dir'],bbnnyyyy), 'wb')
      fileo.write(filei.read())
      filei.close()
      fileo.close()

    except:    #  If the file is not present in the real-time directory, look in the archive.

      src  = '{0}/{1}/b{2}.dat.gz'.format(config.get('best_dir_alt','https://ftp.nhc.noaa.gov/atcf/archive'),bbnnyyyy[4:8],bbnnyyyy)

      #  Unzip the file from the NHC server, write the file to the work directory
      gzfile = gzip.GzipFile(fileobj=urllib.request.urlopen(src))
      uzfile = open('{0}/b{1}.dat'.format(config['work_dir'],bbnnyyyy), 'wb')
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

       self.filetype = config.get('file_format', 'minghua')
       self.nens = int(config['num_ens'])
       self.halfens = int((self.nens - 1) / 2)

       file_name = self.member_name(1)

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
          self.domdict.coords['longitude']  = (self.domdict.coords['longitude'] + 180) % 360 - 180


       self.truelat1 = float(config.get('truelat1', '39.0'))
       self.truelat2 = float(config.get('truelat2', '50.0'))
       self.knowni   = float(config.get('knowni', '0.0'))
       self.knownj   = float(config.get('knownj', '0.0'))
       self.stdlon   = float(config.get('stdlon', '-125.0'))
       self.dx       = float(config.get('dx', '9000.0'))
       self.lat1     = self.domdict.latitude[0,0]
       self.lon1     = self.domdict.longitude[0,0]

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


    def member_name(self, member):

       if member <= self.halfens:
          fmem = '%0.3i' % member
          if self.filetype == 'minghua':
             file_name = "{0}/{1}/ecm{2}/wrfout_d01_{3}_subset_variables.nc".format(self.config['model_dir'],self.datea,fmem,self.fdatestr)
          elif self.filetype == 'westwrf':
             file_name = "{0}/{1}/ecm{2}/wrfcf_d01_{3}.nc".format(self.config['model_dir'],self.datea,fmem,self.fdatestr)
       else:
          fmem = '%0.3i' % (member - self.halfens)
          if self.filetype == 'minghua':
             file_name = "{0}/{1}/gefs{2}/wrfout_d01_{3}_subset_variables.nc".format(self.config['model_dir'],self.datea,fmem,self.fdatestr)
          elif self.filetype == 'westwrf':
             file_name = "{0}/{1}/gefs{2}/wrfcf_d01_{3}.nc".format(self.config['model_dir'],self.datea,fmem,self.fdatestr)

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


    def latlon_to_ij(self, lat, lon):

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
             vdict['pres_start'] = np.where(self.domdict.isobaricInhPa[:]==int(vdict['isobaricInhPa'][0]))[0]
             vdict['pres_end']   = np.where(self.domdict.isobaricInhPa[:]==int(vdict['isobaricInhPa'][1]))[0]

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


    def read_pressure_levels(self, varname, vdict):

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
                                                vdict['j_start']:vdict['j_end'],vdict['i_start']:vdict['i_end']].squeeze()

       #  Read the only level if it is a single level variable
       else:

          vout = ds.get(self.var_dict[varname])[0,vdict['j_start']:vdict['j_end'],vdict['i_start']:vdict['i_end']].squeeze()

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

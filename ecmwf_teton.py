import os
import time
import shutil
import sys
import cfgrib
import datetime as dt
import glob
import gzip
import urllib
import pandas as pd
import numpy as np
import xarray as xr

def stage_grib_files(datea, config):
    '''
    This is a generic class for copying or linking grib file data from a specific location
    to a directory where calculations are performed.  No matter where these calculations are
    carried out, this routine must exist.
 
    This particular instance is for ECMWF data on the machine teton at UAlbany.  In this 
    case, the grib files are linked to files that are located in another directory on that
    machine.

    Attributes:
        datea (string):  The initialization time of the forecast (yyyymmddhh)
        config  (dict):  The dictionary with configuration information
    '''

    freq = config['model'].get('input_hour_int', 6)
    fmax = config['model'].get('fcst_hour_max', 120)

    #  Make the work directory if it does not exist
    if not os.path.isdir(config['locations']['work_dir']):
       try:
          os.makedirs(config['locations']['work_dir'])
       except OSError as e:
          raise e

    init   = dt.datetime.strptime(datea, '%Y%m%d%H')
    init_s = init.strftime("%m%d%H%M")

    #  Loop over all forecast times, link to the source file
    for fhr in range(0, int(fmax)+int(freq), int(freq)):

       datef   = init + dt.timedelta(hours=fhr)
       datef_s = datef.strftime("%m%d%H%M")

       grib_file = 'E1E{0}{1}1'.format(str(init_s), str(datef_s))
       infile    = '{0}/{1}'.format(config['locations']['model_dir'],grib_file)
       outfile   = '{0}/{1}'.format(config['locations']['work_dir'],grib_file)

       #  Only try to copy if the file is not there
       if ( not os.path.isfile(outfile) ):

          #  Wait for the source file to be present 
          while not os.path.exists(infile):
             time.sleep(20.1)

          #  Wait for the file to be finished being copied
          while ( (time.time() - os.path.getmtime(infile)) < 60 ):
             time.sleep(10)

          try:  #  Try to link from the source to the work directory
             os.symlink(infile, outfile)
          except Exception as err:
             raise err


def stage_atcf_files(datea, bbnnyyyy, config):
    '''
    This is a generic routine for copying or linking ATCF file data from a specific location
    to a directory where calculations are performed.  No matter where these calculations are
    carried out, this routine must exist.

    The result is a set of ATCF files in the work directory of the format atcf_NN.dat, 
    where NN is the ensemble member number.

    This particular instance is for ECMWF data on the machine teton at UAlbany.  In this 
    case, all ATCF data for a particular storm is in one file, so the code waits for this
    initialization time to exist, then uses sed to get the lines attributed to each 
    ensemble member and places that data in a seperate file.

    Attributes:
        datea    (string):  The initialization time of the forecast (yyyymmddhh)
        bbnnyyyy (string):  TC Identification, where bb is the basin, nn is the number, yyyy is year
        config     (dict):  The dictionary with configuration information
    '''

    src  = '{0}/a{1}.dat'.format(config['locations']['atcf_dir'],bbnnyyyy)
    nens = int(config['model']['num_ens'])

    if not os.path.isfile('{0}/atcf_{1}.dat'.format(config['locations']['work_dir'],'%0.2i' % nens)):

       #  Wait for the source file to be present 
       while not os.path.exists(src):
          time.sleep(20.5)

       #  Wait for the ensemble ATCF information to be placed in the file
       while ( len(os.popen('sed -ne /{0}/p {1} | sed -ne /EE/p'.format(datea,src)).read()) == 0 ):
          time.sleep(20.7)

       #  Wait for the file to be finished being copied
       while ( (time.time() - os.path.getmtime(src)) < 60 ):
          time.sleep(10)

       for n in range(nens + 1):

          nn = '%0.2i' % n
          file_name = '{0}/atcf_{1}.dat'.format(config['locations']['work_dir'],nn)

          fo = open(file_name,"w")
          fo.write(os.popen('sed -ne /{0}/p {1} | sed -ne /EE{2}/p'.format(datea,src,nn)).read())
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

    This particular instance is for ECMWF data on the UAlbany teton machine.  Each ECMWF file contains all
    members for a specific forecast time.  

    Attributes:
        datea (string):  The initialization time of the forecast (yyyymmddhh)
        fhr      (int):  Forecast hour
        config  (dict):  The dictionary with configuration information
    '''

    def __init__(self, datea, fhr, config):

        init    = dt.datetime.strptime(datea, '%Y%m%d%H')
        init_s  = init.strftime("%m%d%H%M")
        datef   = init + dt.timedelta(hours=fhr)
        datef_s = datef.strftime("%m%d%H%M")

        #  Construct the grib file dictionary for a particular forecast hour
        file_name = os.path.join(config['locations']['work_dir'], "E1E{0}{1}1".format(str(init_s), str(datef_s)))
        try:  
           ds = cfgrib.open_datasets(file_name)
           self.grib_dict = {}
           for d in ds:
              for tt in d:
                 if 'number' in d[tt].dims:
                    self.grib_dict.update({'{0}_pf'.format(tt): d[tt]})
                 else:
                    self.grib_dict.update({'{0}_cf'.format(tt): d[tt]})

        except IOError as exc:
           raise RuntimeError('Failed to open {0}'.format(file_name)) from exc

        #  This is a dictionary that maps from generic variable names to the name of variable in file
        self.var_dict = {'zonal_wind': 'u',           \
                         'meridional_wind': 'v',      \
                         'zonal_wind_10m': 'u10',      \
                         'meridional_wind_10m': 'v10', \
                         'geopotential_height': 'gh', \
                         'temperature': 't',          \
                         'relative_humidity': 'r',    \
                         'specific_humidity': 'q',    \
                         'sea_level_pressure': 'msl', \
                         'iwv': 'tcwv',               \
                         'precipitation': 'tp',       \
                         'precip_type': 'ptype'}

        for key in self.grib_dict:
           if np.max(self.grib_dict[key].coords['longitude']) > 180.:
              self.grib_dict[key] = self.grib_dict[key].assign_coords(longitude=(((self.grib_dict[key].longitude + 180.) % 360.) - 180.)).sortby('longitude')

        if config['model'].get('flip_lon','False') == 'True':
           for key in self.grib_dict:
              self.grib_dict[key] = self.grib_dict[key].assign_coords(longitude=((self.grib_dict[key].coords['longitude'] + 360.) % 360.)).sortby('longitude')

        if '{0}_cf'.format(self.var_dict['specific_humidity']) in self.grib_dict:
           self.has_specific_humidity = True
        else:
           self.has_specific_humidity = False

        self.has_frozen_fraction = False
        self.has_precip_category = False
        if '{0}_cf'.format(self.var_dict['precip_type']) in self.grib_dict:
           self.has_precip_category = True

        self.has_total_precip = True

        for var in self.var_dict.values():
           varname = '{0}_pf'.format(var)
           if varname in self.grib_dict:
              self.nens = int(self.grib_dict[varname].attrs['GRIB_totalNumber'])
              break


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

       vname = '{0}_cf'.format(self.var_dict[varname])

       if 'latitude' in vdict:

          #  See if the latitude values are ordered reversed
          if float(self.grib_dict[vname].attrs['GRIB_latitudeOfFirstGridPointInDegrees']) > \
             float(self.grib_dict[vname].attrs['GRIB_latitudeOfLastGridPointInDegrees']):
             vdict['lat_start'] = float(vdict['latitude'][1])
             vdict['lat_end']   = float(vdict['latitude'][0])
          else:
             vdict['lat_start'] = float(vdict['latitude'][0])
             vdict['lat_end']   = float(vdict['latitude'][1])

       else:

          #  Get all latitude values
          latvec = list(self.grib_dict[vname].latitude.data)
          vdict['lat_start'] = latvec[0]
          vdict['lat_end']   = latvec[-1]

       if 'longitude' in vdict:

          #  Take longitude from input values
          vdict['lon_start'] = float(vdict['longitude'][0])
          vdict['lon_end']   = float(vdict['longitude'][1])

       else:

          #  Use all longitude values
          lonvec = list(self.grib_dict[vname].longitude.data)
          vdict['lon_start'] = lonvec[0]
          vdict['lon_end']   = lonvec[-1]

       if 'isobaricInhPa' in vdict:

          #  See of pressure level values are reversed
          if self.grib_dict[vname].isobaricInhPa[0] > self.grib_dict[vname].isobaricInhPa[1]:
            vdict['pres_start'] = float(vdict['isobaricInhPa'][1])
            vdict['pres_end']   = float(vdict['isobaricInhPa'][0])
          else:
            vdict['pres_start'] = float(vdict['isobaricInhPa'][0])
            vdict['pres_end']   = float(vdict['isobaricInhPa'][1])

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

       vname = '{0}_cf'.format(self.var_dict[varname])

       #  Create attributes based on what is in the file
       attrlist = {}
       if 'description' in vdict:
         attrlist['description'] = vdict['description']
       if 'units' in vdict:
         attrlist['units'] = vdict['units']
       if '_FillValue' in vdict:
         attrlist['_FillValue'] = vdict['_FillValue']

       #  Create an empty xarray that can be used to copy data into
       lonvec = list(self.grib_dict[vname].sel(longitude=slice(vdict['lon_start'], vdict['lon_end'])).longitude.data)
       latvec = list(self.grib_dict[vname].sel(latitude=slice(vdict['lat_start'], vdict['lat_end'])).latitude.data)
       ensarr = xr.DataArray(name='ensemble_data', data=np.zeros([nens, len(latvec), len(lonvec)]),
                             dims=['ensemble', 'latitude', 'longitude'], attrs=attrlist, 
                             coords={'ensemble': [i for i in range(nens)], 'latitude': latvec, 'longitude': lonvec}) 

       return(ensarr)


    def read_pressure_levels(self, varname):

       vname = '{0}_cf'.format(self.var_dict[varname])
       return self.grib_dict[vname].isobaricInhPa[:]


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

       #  Read a single pressure level of data, if this is a variable that has pressure levels
       if 'isobaricInhPa' in vdict:

          if member == 0:  #  Control member
             vname = '{0}_cf'.format(self.var_dict[varname])
             vout  = self.grib_dict[vname].sel(latitude=slice(vdict['lat_start'], vdict['lat_end']),  \
                                               longitude=slice(vdict['lon_start'], vdict['lon_end']), \
                                               isobaricInhPa=slice(vdict['pres_start'], vdict['pres_end']))

          else:
             vname = '{0}_pf'.format(self.var_dict[varname])
             vout  = self.grib_dict[vname].sel(number=member,                                         \
                                               latitude=slice(vdict['lat_start'], vdict['lat_end']),  \
                                               longitude=slice(vdict['lon_start'], vdict['lon_end']), \
                                               isobaricInhPa=slice(vdict['pres_start'], vdict['pres_end']))

       #  Read the only level if it is a single level variable
       else:

          if member == 0:  #  Control member
             vname = '{0}_cf'.format(self.var_dict[varname])
             vout  = self.grib_dict[vname].sel(latitude=slice(vdict['lat_start'], vdict['lat_end']), \
                                               longitude=slice(vdict['lon_start'], vdict['lon_end']))

          else:

             vname = '{0}_pf'.format(self.var_dict[varname])
             vout  = self.grib_dict[vname].sel(number=member,                                        \
                                               latitude=slice(vdict['lat_start'], vdict['lat_end']), \
                                               longitude=slice(vdict['lon_start'], vdict['lon_end']))

       if 'valid_time' in list(vout.coords):
          return(vout.reset_coords('valid_time', drop=True))
       else:
          return(vout)


    def close_files(self):

       del self.grib_dict


    def read_static_field(self, static_file, varname, vdict):
       '''
       This is a generic routine that is used to read a static field from a file from 
       a user-provided file based on the information contained within the dictionary vdict.

       Attributes:
           static_file (string):  Name of file to extract variable from 
           varname     (string):  Name of the variable that will be extracted from file
           vdict         (dict):  The dictionary object with variable information
       '''

       static_dict = {'landmask': 'lsm'}

       ds = cfgrib.open_datasets(static_file)
       grib_dict = {}
       for d in ds:
          for tt in d:
             grib_dict.update({'{0}'.format(tt): d[tt]})

       if vdict.get('flip_lon','False') == 'True':
          for key in grib_dict:
             grib_dict[key] = grib_dict[key].assign_coords(longitude=((grib_dict[key].coords['longitude'] + 360.) % 360.)).sortby('longitude')

       vout = grib_dict[static_dict[varname]].sel(latitude=slice(vdict['lat_start'], vdict['lat_end']), \
                                                  longitude=slice(vdict['lon_start'], vdict['lon_end']))

       return(vout)


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

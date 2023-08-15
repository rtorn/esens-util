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

#    freq = config.get('fcst_hour_int', 12)
    freq = 6
    fmax = config.get('fcst_hour_max', 120)

    init    = dt.datetime.strptime(datea, '%Y%m%d%H')
    init_s  = init.strftime("%y%j%H")

    #  Make the work directory if it does not exist
    if not os.path.isdir(config['work_dir']):
       try:
          os.makedirs(config['work_dir'])
       except OSError as e:
          raise e

    for fhr in range(0, int(fmax)+int(freq), int(freq)):

       fbase = "{0}000{1}".format(init_s, '%0.3i' % fhr)
       fout  = '{0}/f{1}.grb2'.format(config['work_dir'],'%0.3i' % fhr)

       if not os.path.isfile(fout):

          #  Wait for the source file to be present 
          infile = '{0}/gec00/{1}'.format(config['model_dir'],fbase)
          while not os.path.exists(infile):
             time.sleep(20.1)

          #  Wait for the file to be finished being copied
          while ( (time.time() - os.path.getmtime(infile)) < 60 ):
             time.sleep(10)

          os.system('cat {0}/gec00/{1} {0}/gep*/{1} >& {2}'.format(config['model_dir'],fbase,fout))

          for n in range(int(config['num_ens'])+1):

             #  Construct the grib file dictionary for a particular forecast hour
             if n > 0:
                fdir = "gep{0}".format('%0.2i' % n)
             else:
                fdir = "gec00"

             #  read a few extra fields from the alternate files
             falt = '{0}/pgrb2b/{1}/{2}'.format(config['model_dir'],fdir,fbase)
             os.system('wgrib2 -s {0} | grep -e \"TMP:300 mb\" -e \"TMP:400 mb\" -e \"RH:300 mb\" -e \"RH:400 mb\" | wgrib2 -fix_ncep -i -append {0} -grib {1}'.format(falt,fout))

    #  Grab precipitation ooutput at more frequent intervals
    fmin = config.get('precip_hour_min', 6)
    freq = config.get('precip_hour_int', 12)
    fmax = config.get('precip_hour_max', -12)  

    for fhr in range(fmin, int(fmax)+int(freq), int(freq)):

       fbase = "{0}000{1}".format(init_s, '%0.3i' % fhr)
       fout  = '{0}/f{1}.grb2'.format(config['work_dir'],'%0.3i' % fhr)

       if not os.path.isfile(fout):

          for n in range(int(config['num_ens'])+1):

             #  Construct the grib file dictionary for a particular forecast hour
             if n > 0:
                fdir = "gep{0}".format('%0.2i' % n)
             else:
                fdir = "gec00"

             #  read a few extra fields from the alternate files
             falt = '{0}/{1}/{2}'.format(config['model_dir'],fdir,fbase)
             os.system('wgrib2 -s {0} | grep -e \"APCP\" | wgrib2 -fix_ncep -i -append {0} -grib {1}'.format(falt,fout))


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

#    src  = '{0}/a{1}.dat.gz'.format(config['atcf_dir'],bbnnyyyy)
#    nens = int(config['num_ens'])

    #  Unzip the file from the NHC server, write the file to the work directory
#    gzfile = gzip.GzipFile(fileobj=urllib.request.urlopen(src))
#    uzfile = open('{0}/a{1}.dat'.format(config['work_dir'],bbnnyyyy), 'wb')
#    uzfile.write(gzfile.read())
#    gzfile.close()
#    uzfile.close()

    init    = dt.datetime.strptime(datea, '%Y%m%d%H')
    datef   = init + dt.timedelta(hours=6)

    src  = '{0}/{1}.a{2}.dat'.format(config['atcf_dir'],datef.strftime('%Y%m%d%H'),bbnnyyyy)
    print('ATCF file',datef.strftime('%Y%m%d%H'),src)

#    #  Wait for the ensemble ATCF information to be placed in the file
#    while ( len(os.popen('sed -ne /AP/p {0}/{1}.a{2}.dat'.format(config['atcf_dir'],datea,bbnnyyyy)).read()) == 0 ):
#       time.sleep(20.7)

#    #  Wait for the file to be finished being copied
#    while ( (time.time() - os.path.getmtime(src)) < 60 ):
#       time.sleep(10)

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
#          fo.write(os.popen('sed -ne /{0}/p {1}/a{2}.dat | sed -ne /{3}{4}/p'.format(datea,config['work_dir'],bbnnyyyy,modid,nn)).read())
          fo.write(os.popen('sed -ne /{0}/p {1} | sed -ne /{2}{3}/p'.format(datea,src,modid,nn)).read())
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

        fbase = "{0}000{1}".format(init_s, '%0.3i' % fhr)
        self.grib_dict = {}

        #  Construct the grib file dictionary for a particular forecast hour
        file_name = "{0}/f{1}.grb2".format(config['work_dir'], '%0.3i' % fhr)
        try:  
           ds = cfgrib.open_datasets(file_name)

           for d in ds:
              for tt in d:
                 if d[tt].attrs['GRIB_dataType'] == 'pf':
                    self.grib_dict.update({'{0}_pf'.format(tt): d[tt]})
                 else:
                    self.grib_dict.update({'{0}_cf'.format(tt): d[tt]})

        except IOError as exc:
           raise RuntimeError('Failed to open {0}'.format(file_name)) from exc

        #  This is a dictionary that maps from generic variable names to the name of variable in file
        self.var_dict = {'zonal_wind': 'u',             \
                         'meridional_wind': 'v',        \
                         'zonal_wind_10m': 'u10',       \
                         'meridional_wind_10m': 'v10',  \
                         'geopotential_height': 'gh',   \
                         'temperature': 't',            \
                         'specific_humidity': 'q',      \
                         'relative_humidity': 'r',      \
                         'sea_level_pressure': 'prmsl', \
                         'precipitation': 'tp' }

        for key in self.grib_dict:
           if np.max(self.grib_dict[key].coords['longitude']) > 180:
              self.grib_dict[key].coords['longitude']  = (self.grib_dict[key].coords['longitude'] + 180) % 360 - 180
              self.grib_dict[key] = self.grib_dict[key].sortby('longitude')

        if config.get('flip_lon','False') == 'True':
           for key in self.grib_dict:
              self.grib_dict[key].coords['longitude'] = (self.grib_dict[key].coords['longitude'] + 360.) % 360.
              self.grib_dict[key] = self.grib_dict[key].sortby('longitude')

        if '{0}_cf'.format(self.var_dict['specific_humidity']) in self.grib_dict:
           self.has_specific_humidity = True
        else:
           self.has_specific_humidity = False

        self.nens = int(self.grib_dict['gh_pf'].attrs['GRIB_totalNumber']) + 1

        self.has_total_precip = False


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
       st_dict = {}
       for d in ds:
          for tt in d:
             st_dict.update({'{0}'.format(tt): d[tt]})

       for key in st_dict:
           if np.max(st_dict[key].coords['longitude']) > 180:
              st_dict[key] = st_dict[key].assign_coords(longitude=((st_dict[key].coords['longitude'] + 180) % 360 - 180)).sortby('longitude')

       if vdict.get('flip_lon','False') == 'True':
          for key in st_dict:
             st_dict[key] = st_dict[key].assign_coords(longitude=((st_dict[key].coords['longitude'] + 360.) % 360.)).sortby('longitude')

       vout = st_dict[static_dict[varname]].sel(latitude=slice(vdict['lat_start'], vdict['lat_end']), \
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

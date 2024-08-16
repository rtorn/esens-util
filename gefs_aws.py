import os
import time
import gzip
import sys
import warnings
import urllib
import datetime as dt
import numpy as np
import xarray as xr
import cfgrib
from herbie import Herbie

def stage_grib_files(datea, config):

    '''
    This is a generic class for copying or linking grib file data from a specific location
    to a directory where calculations are performed.  No matter where these calculations are
    carried out, this routine must exist.

    This particular instance employs GEFS data that is in the AWS cloud.  As a 
    consequence, this routine is a placeholder and does not need to do anything.

    Attributes:
        datea (string):  The initialization time of the forecast (yyyymmddhh)
        config  (dict):  The dictionary with configuration information
    '''

    #  Make the work directory if it does not exist
    if not os.path.isdir(config['locations']['work_dir']):
       try:
          os.makedirs(config['locations']['work_dir'])
       except OSError as e:
          raise e


def stage_atcf_files(datea, bbnnyyyy, config):
    '''
    This is a generic class for copying or linking ATCF file data from a specific location
    to a directory where calculations are performed.  No matter where these calculations are
    carried out, this routine must exist.

    The result is a set of ATCF files in the work directory of the format atcf_NN.dat,
    where NN is the ensemble member number.

    This particular instance is for GEFS data that is on the NHC FTP server. In this
    case, all ATCF data for a particular storm is in one file, so the code unzips and
    reads the storm file, then uses sed to get the lines attributed to each
    ensemble member and places that data in a seperate file.

    Attributes:
        datea    (string):  The initialization time of the forecast (yyyymmddhh)
        bbnnyyyy (string):  TC Identification, where bb is the basin, nn is the number, yyyy is year
        config     (dict):  The dictionary with configuration information
    '''

    src  = '{0}/a{1}.dat.gz'.format(config['locations']['atcf_dir'],bbnnyyyy)
    nens = int(config['model']['num_ens'])

    #  Unzip the file from the NHC server, write the file to the work directory
    gzfile = gzip.GzipFile(fileobj=urllib.request.urlopen(src))
    uzfile = open('{0}/a{1}.dat'.format(config['locations']['work_dir'],bbnnyyyy), 'wb')
    uzfile.write(gzfile.read())        
    gzfile.close()
    uzfile.close()

    #  Wait for the ensemble ATCF information to be placed in the file
#    while ( len(os.popen("sed -ne /" + self.init + "/p " + self.dest + "/a" + bbnnyyyy + ".dat | sed -ne /EE/p").read()) == 0 ):
#       time.sleep(20.7)

    #  Wait for the file to be finished being copied
#    while ( (time.time() - os.path.getmtime(self.src)) < 60 ):
#       time.sleep(10)

    for n in range(nens + 1):

       if ( n > 0 ):
          modid = 'AP'
       else:
          modid = 'AC'

       nn = '%0.2i' % n
       file_name = '{0}/atcf_{1}.dat'.format(config['locations']['work_dir'],nn)

       #  If the specific member's ATCF file does not exist, copy from the source file with sed.
       if not os.path.isfile(file_name):

          fo = open(file_name,"w")
          fo.write(os.popen('sed -ne /{0}/p {1}/a{2}.dat | sed -ne /{3}{4}/p'.format(datea,config['locations']['work_dir'],bbnnyyyy,modid,nn)).read())
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

    This particular instance is for GEFS data within Amazon Web Services.  This instance uses the Herbie 
    package to read subsets of individual member files and times.  

    Attributes:
        datea (string):  The initialization time of the forecast (yyyymmddhh)
        fhr      (int):  Forecast hour
        config  (dict):  The dictionary with configuration information
    '''

    def __init__(self, datea, fhr, config):

        self.init  = dt.datetime.strptime(datea, '%Y%m%d%H')
        self.datef = self.init + dt.timedelta(hours=fhr)
        self.config = config

        self.memlist = ['c00']
        for n in range(1,int(config['model']['num_ens'])+1):
          self.memlist.append('p{0}'.format('%0.2i' % n))

        self.grib_dict = {}
        self.grib_dictb = {}

        #  Create Herbie object for each primary and secondary ensemble file
        for n in range(len(self.memlist)):

           H = Herbie('{0}-{1}-{2} {3}:00'.format(datea[0:4],datea[4:6],datea[6:8],datea[8:10]), model='gefs', member=n, product='atmos.5', fxx=fhr)
           self.grib_dict.update({self.memlist[n]: H})

           H = Herbie('{0}-{1}-{2} {3}:00'.format(datea[0:4],datea[4:6],datea[6:8],datea[8:10]), model='gefs', member=n, product='atmos.5b', fxx=fhr)
           self.grib_dictb.update({self.memlist[n]: H})
    
        #  This is a dictionary that maps from generic variable names to the name of variable in file
        self.var_dict = {'zonal_wind':          {'hname': ':UGRD:\d+ mb',     'dname': 'u'    }, \
                         'meridional_wind':     {'hname': ':VGRD:\d+ mb',     'dname': 'v'    }, \
                         'zonal_wind_10m':      {'hname': ':UGRD:10 m above', 'dname': 'u10'  }, \
                         'meridional_wind_10m': {'hname': ':VGRD:10 m above', 'dname': 'v10'  }, \
                         'geopotential_height': {'hname': ':HGT:\d+ mb',      'dname': 'gh'   }, \
                         'temperature':         {'hname': 'TMP:\d+ mb',       'dname': 't'    }, \
                         'relative_humidity':   {'hname': 'RH:\d+ mb',        'dname': 'r'    }, \
                         'specific_humidity':   {'hname': 'Q:\d+ mb',         'dname': 'q'    }, \
                         'sea_level_pressure':  {'hname': 'PRMSL',            'dname': 'prmsl'}, \
                         'precipitation':       {'hname': 'APCP',             'dname': 'tp'   }, \
                         'precip_type':         {'hname': 'PTYPE',            'dname': 'ptype'}}


        self.has_specific_humidity = False

        self.has_total_precip = False

        self.nens = len(self.memlist)


    def set_longitude(self, datf):

       if np.max(datf.coords['longitude']) > 180.:
          datf.coords['longitude'] = (datf.coords['longitude'] + 180) % 360 - 180
          datf = datf.sortby(datf.longitude)

       if eval(self.config['model'].get('flip_lon','False')):
          datf.coords['longitude'] = (datf.coords['longitude'] + 360.) % 360.
          datf = datf.sortby(datf.longitude)

       return(datf)


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

       vname = self.var_dict[varname]
       df = self.set_longitude(self.grib_dict[self.memlist[0]].xarray(vname['hname']))

       df2 = self.grib_dict[self.memlist[0]].xarray(self.var_dict['zonal_wind']['hname'])

       if 'latitude' in vdict:

          if float(df.latitude.values[0]) > float(df.latitude.values[-1]):
             vdict['lat_start'] = int(vdict['latitude'][1])
             vdict['lat_end']   = int(vdict['latitude'][0])
          else:
             vdict['lat_start'] = int(vdict['latitude'][0])
             vdict['lat_end']   = int(vdict['latitude'][1])

       else:

          vdict['lat_start'] = df.latitude.values[0]
          vdict['lat_end']   = df.latitude.values[-1]

       vdict['num_lat'] = len(df.latitude.sel(latitude=slice(vdict['lat_start'], vdict['lat_end'])).values)


       if 'longitude' in vdict:

          vdict['lon_start'] = int(vdict['longitude'][0])
          vdict['lon_end']   = int(vdict['longitude'][1])

       else:

          vdict['lon_start'] = df.longitude.values[0]
          vdict['lon_end']   = df.longitude.values[-1]

       vdict['num_lon'] = len(df.longitude.sel(longitude=slice(vdict['lon_start'], vdict['lon_end'])).values)


       if 'isobaricInhPa' in vdict:

          if float(df2.isobaricInhPa.values[0]) > float(df2.isobaricInhPa.values[-1]):
            vdict['pres_start'] = int(vdict['isobaricInhPa'][1])
            vdict['pres_end']   = int(vdict['isobaricInhPa'][0])
          else:
            vdict['pres_start'] = int(vdict['isobaricInhPa'][0])
            vdict['pres_end']   = int(vdict['isobaricInhPa'][1])

          presall = list(df2.isobaricInhPa.sel(isobaricInhPa=slice(vdict['pres_start'], vdict['pres_end'])).values)
          vdict['pres_list'] = presall
          vdict['num_pres'] = len(vdict['pres_list'])

          presvar = list(df.isobaricInhPa.sel(isobaricInhPa=slice(vdict['pres_start'], vdict['pres_end'])).values)

          if presvar == presall:
  
             vdict['pres_match'] = True

          else:

             vdict['pres_match'] = False

             haspres = []
             for i in range(len(presall)):
               if presall[i] in presvar:
                 haspres.append(True)
               else:
                 haspres.append(False)

             vdict['pres_primary'] = haspres

       del df

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

       vname = self.var_dict[varname]
       df = self.set_longitude(self.grib_dict[self.memlist[0]].xarray(vname['hname']))

       #  Create an empty xarray that can be used to copy data into
       lonvec = df.longitude.sel(longitude=slice(vdict['lon_start'], vdict['lon_end'])).values
       latvec = df.latitude.sel(latitude=slice(vdict['lat_start'], vdict['lat_end'])).values
       ensarr = xr.DataArray(name='ensemble_data', data=np.zeros([nens, len(latvec), len(lonvec)]),
                             dims=['ensemble', 'latitude', 'longitude'], attrs=attrlist, 
                             coords={'ensemble': [i for i in range(nens)], 'latitude': latvec, 'longitude': lonvec}) 

       del df

       return(ensarr)


    def read_pressure_levels(self, varname):

       return self.grib_dict[self.memlist[0]].xarray(self.var_dict['zonal_wind']['hname']).isobaricInhPa.values


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

       vname = self.var_dict[varname]
       df = self.set_longitude(self.grib_dict[self.memlist[member]].xarray(vname['hname']))

       #  Read a single pressure level of data, if this is a variable that has pressure levels
       if 'isobaricInhPa' in vdict:

          #  Rename the lat, lon, and pressure arrays to common convention
          if vdict['pres_match']:

             vout  = df[vname['dname']].sel(latitude=slice(vdict['lat_start'], vdict['lat_end']),   \
                                            longitude=slice(vdict['lon_start'], vdict['lon_end']),   \
                                            isobaricInhPa=slice(vdict['pres_start'], vdict['pres_end'])).squeeze()

          else:

             dfb = self.grib_dictb[self.memlist[member]].xarray(vname['hname'])
             dfb = self.set_longitude(dfb)

             lonvec = df.longitude.sel(longitude=slice(vdict['lon_start'], vdict['lon_end'])).values
             latvec = df.latitude.sel(latitude=slice(vdict['lat_start'], vdict['lat_end'])).values
             vout = xr.DataArray(data=np.zeros([vdict['num_pres'], len(latvec), len(lonvec)]), \
                                 dims=['isobaricInhPa', 'latitude', 'longitude'], \
                                 coords={'isobaricInhPa': vdict['pres_list'], 'latitude': latvec, 'longitude': lonvec})

             for k in range(int(vdict['num_pres'])):

                if vdict['pres_primary'][k]:

                   vout[k,:,:] = df[vname['dname']].sel(latitude=slice(vdict['lat_start'], vdict['lat_end']),   \
                                                        longitude=slice(vdict['lon_start'], vdict['lon_end']),   \
                                                        isobaricInhPa=slice(vdict['pres_list'][k], vdict['pres_list'][k])).squeeze()

                else:

                   vout[k,:,:] = dfb[vname['dname']].sel(latitude=slice(vdict['lat_start'], vdict['lat_end']),   \
                                                         longitude=slice(vdict['lon_start'], vdict['lon_end']),   \
                                                         isobaricInhPa=slice(vdict['pres_list'][k], vdict['pres_list'][k])).squeeze()


       #  Read the only level if it is a single level variable
       else:

          #  Rename the lat and lon arrays to common convention
          vout  = df[vname['dname']].sel(latitude=slice(vdict['lat_start'], vdict['lat_end']), \
                                         longitude=slice(vdict['lon_start'], vdict['lon_end'])).squeeze()

       return(vout)


    def close_files(self):

       del self.grib_dict, self.grib_dictb


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


    @staticmethod
    def __check_file_exists(filename):
        isfile = False
        try:
            if os.path.isfile(filename):
                isfile = True
        except Exception as err:
            print(err)

        return isfile

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

import os
import time
import gzip
import sys
import warnings
import urllib
import datetime as dt
import numpy as np
import xarray as xr

def stage_grib_files(datea, config):

    '''
    This is a generic class for copying or linking grib file data from a specific location
    to a directory where calculations are performed.  No matter where these calculations are
    carried out, this routine must exist.

    This particular instance employs GEFS OpeNDAP data that is on NCEP servers.  As a 
    consequence, this routine is a placeholder and does not need to do anything.

    Attributes:
        datea (string):  The initialization time of the forecast (yyyymmddhh)
        config  (dict):  The dictionary with configuration information
    '''

    #  Make the work directory if it does not exist
    if not os.path.isdir(config['work_dir']):
       try:
          os.makedirs(config['work_dir'])
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

    src  = '{0}/a{1}.dat.gz'.format(config['atcf_dir'],bbnnyyyy)
    nens = int(config['num_ens'])

    #  Unzip the file from the NHC server, write the file to the work directory
    gzfile = gzip.GzipFile(fileobj=urllib.request.urlopen(src))
    uzfile = open('{0}/a{1}.dat'.format(config['work_dir'],bbnnyyyy), 'wb')
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
       file_name = '{0}/atcf_{1}.dat'.format(config['work_dir'],nn)

       #  If the specific member's ATCF file does not exist, copy from the source file with sed.
       if not os.path.isfile(file_name):

          fo = open(file_name,"w")
          fo.write(os.popen('sed -ne /{0}/p {1}/a{2}.dat | sed -ne /{3}{4}/p'.format(datea,config['work_dir'],bbnnyyyy,modid,nn)).read())
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

    This particular instance is for GEFS OpenNDAP data.  This instance reads the file that contains all 
    ensemble members and times.  The init routine reads the grib file dictionary.

    Attributes:
        datea (string):  The initialization time of the forecast (yyyymmddhh)
        fhr      (int):  Forecast hour
        config  (dict):  The dictionary with configuration information
    '''

    def __init__(self, datea, fhr, config):

        self.init  = dt.datetime.strptime(datea, '%Y%m%d%H')
        self.datef = self.init + dt.timedelta(hours=fhr)

        with warnings.catch_warnings():
            warnings.simplefilter(
                "ignore", category=xr.coding.variables.SerializationWarning)

            file_name = '{0}/gefs{1}/gefs_pgrb2ap5_all_{2}z'.format(config['model_dir'], datea[0:8], datea[8:10])

            try:
                self.ds_dict = xr.open_dataset(file_name)
            except (OSError, RuntimeError):
                print("Read Error; Sleeping before trying again")
                time.sleep(240)
                self.ds_dict = xr.open_dataset(file_name)
            except IOError as exc:
                raise RuntimeError('Failed to open {0}'.format(file_name)) from exc

            file_name = '{0}/gefs{1}/gefs_pgrb2bp5_all_{2}z'.format(config['model_dir'], datea[0:8], datea[8:10])

            try:
                self.ds_dictb = xr.open_dataset(file_name)
            except (OSError, RuntimeError):
                print("Read Error; Sleeping before trying again")
                time.sleep(240)
                self.ds_dictb = xr.open_dataset(file_name)
            except IOError as exc:
                raise RuntimeError('Failed to open {0}'.format(file_name)) from exc


        #  Put longitude into -180 to 180 format
        self.ds_dict.coords['lon']  = (self.ds_dict.coords['lon'] + 180.) % 360. - 180.
        self.ds_dict = self.ds_dict.sortby(self.ds_dict.lon)
#        self.ds_dictb.coords['lon'] = ((self.ds_dictb.coords['lon'] + 180.) % 360. - 180.).sortby('lon')
#        for key in self.ds_dict:
#           if np.max(self.ds_dict[key].coords['lon']) > 180:
#              self.ds_dict[key] = self.ds_dict[key].assign_coords(lon=(((self.ds_dict[key].lon + 180) % 360) - 180)).sortby('lon')

#        for key in self.ds_dictb:
#           if np.max(self.ds_dictb[key].coords['lon']) > 180:
#              self.ds_dictb[key] = self.ds_dictb[key].assign_coords(lon=(((self.ds_dictb[key].lon + 180) % 360) - 180)).sortby('lon')

        print(self.ds_dict)

        #  This is a dictionary that maps from generic variable names to the name of variable in file
        self.var_dict = {'zonal_wind': 'ugrdprs',          \
                         'meridional_wind': 'vgrdprs',     \
                         'geopotential_height': 'hgtprs',  \
                         'temperature': 'tmpprs',          \
                         'relative_humidity': 'rhprs',     \
                         'sea_level_pressure': 'prmslmsl', \
                         'precipitation': 'apcpsfc' }

        #print(self.ds_dict)

        self.has_specific_humidity = False

        self.has_total_precip = False

        self.nens = int(len(self.ds_dict['hgtprs'].coords['ens']))


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

       if 'latitude' in vdict:

          if float(self.ds_dict[vname].lat.data[0]) > float(self.ds_dict[vname].lat.data[-1]):
             vdict['lat_start'] = int(vdict['latitude'][1])
             vdict['lat_end']   = int(vdict['latitude'][0])
          else:
             vdict['lat_start'] = int(vdict['latitude'][0])
             vdict['lat_end']   = int(vdict['latitude'][1])

       else:

          try:
             latvec = list(self.ds_dict[vname].lat.data)
          except (OSError, RuntimeError):
             time.sleep(240)
             latvec = list(self.ds_dict[vname].lat.data) 
          vdict['lat_start'] = latvec[0]
          vdict['lat_end']   = latvec[-1]

       if 'longitude' in vdict:

          vdict['lon_start'] = int(vdict['longitude'][0])
          vdict['lon_end']   = int(vdict['longitude'][1])

       else:

          try:
             lonvec = list(self.ds_dict[vname].lon.data)
          except (OSError, RuntimeError):
             time.sleep(240)
             lonvec = list(self.ds_dict[vname].lon.data)
          vdict['lon_start'] = lonvec[0]
          vdict['lon_end']   = lonvec[-1]

       if 'isobaricInhPa' in vdict:

          if float(self.ds_dict[vname].lev.data[0]) > float(self.ds_dict[vname].lev.data[-1]):
            vdict['pres_start'] = int(vdict['isobaricInhPa'][1])
            vdict['pres_end']   = int(vdict['isobaricInhPa'][0])
          else:
            vdict['pres_start'] = int(vdict['isobaricInhPa'][0])
            vdict['pres_end']   = int(vdict['isobaricInhPa'][1])

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
       try:
          lonvec = list(self.ds_dict[self.var_dict[varname]].sel(lon=slice(vdict['lon_start'], vdict['lon_end'])).lon.data).copy()
       except (OSError, RuntimeError):
          time.sleep(240)
          lonvec = list(self.ds_dict[self.var_dict[varname]].sel(lon=slice(vdict['lon_start'], vdict['lon_end'])).lon.data).copy()
       try:
          latvec = list(self.ds_dict[self.var_dict[varname]].sel(lat=slice(vdict['lat_start'], vdict['lat_end'])).lat.data).copy()
       except (OSError, RuntimeError):
          time.sleep(240)
          latvec = list(self.ds_dict[self.var_dict[varname]].sel(lat=slice(vdict['lat_start'], vdict['lat_end'])).lat.data).copy()
       ensarr = xr.DataArray(name='ensemble_data', data=np.zeros([nens, len(latvec), len(lonvec)]),
                             dims=['ensemble', 'latitude', 'longitude'], attrs=attrlist, 
                             coords={'ensemble': [i for i in range(nens)], 'latitude': latvec, 'longitude': lonvec}) 

       return(ensarr)


    def read_pressure_levels(self, varname):

       return self.ds_dict[self.var_dict[varname]].lev[:]


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

       with warnings.catch_warnings():

          warnings.simplefilter("ignore", category=xr.coding.variables.SerializationWarning)

          vname = self.var_dict[varname]

          #  Read a single pressure level of data, if this is a variable that has pressure levels
          if 'isobaricInhPa' in vdict:

             #  Rename the lat, lon, and pressure arrays to common convention
             try:
                vout  = self.ds_dict[vname].sel(lat=slice(vdict['lat_start'], vdict['lat_end']),   \
                                                lon=slice(vdict['lon_start'], vdict['lon_end']),   \
                                                lev=slice(vdict['pres_start'], vdict['pres_end']), \
                                                time=slice(self.datef, self.datef),                \
                                                ens=slice(member+1, member+1)).squeeze().rename({'lon': 'longitude','lat': 'latitude','lev': 'isobaricInhPa'})
             except (OSError, RuntimeError):
                print("Read Field Error; Sleeping before trying again")
                time.sleep(240)      
                vout  = self.ds_dict[vname].sel(lat=slice(vdict['lat_start'], vdict['lat_end']),   \
                                                lon=slice(vdict['lon_start'], vdict['lon_end']),   \
                                                lev=slice(vdict['pres_start'], vdict['pres_end']), \
                                                time=slice(self.datef, self.datef),                \
                                                ens=slice(member+1, member+1)).squeeze().rename({'lon': 'longitude','lat': 'latitude','lev': 'isobaricInhPa'}).copy(deep=True)


             if vdict['pres_start'] != vdict['pres_end']:

                for k in range(len(vout[:,0,0])):

                   if np.isnan(vout[k,0,0]):

                      try:
                         vout[k,:,:] = self.ds_dictb[vname].sel(lat=slice(vdict['lat_start'], vdict['lat_end']),   \
                                                    lon=slice(vdict['lon_start'], vdict['lon_end']),   \
                                                    lev=slice(vout.isobaricInhPa.values[k], vout.isobaricInhPa.values[k]), \
                                                    time=slice(self.datef, self.datef),                \
                                                    ens=slice(member+1, member+1)).squeeze().rename({'lon': 'longitude','lat': 'latitude','lev': 'isobaricInhPa'})
                      except (OSError, RuntimeError):
                         print("Read Field Error; Sleeping before trying again")
                         time.sleep(240)
                         vout[k,:,:] = self.ds_dictb[vname].sel(lat=slice(vdict['lat_start'], vdict['lat_end']),   \
                                                    lon=slice(vdict['lon_start'], vdict['lon_end']),   \
                                                    lev=slice(vout.isobaricInhPa.values[k], vout.isobaricInhPa.values[k]), \
                                                    time=slice(self.datef, self.datef),                \
                                                    ens=slice(member+1, member+1)).squeeze().rename({'lon': 'longitude','lat': 'latitude','lev': 'isobaricInhPa'}).copy(deep=True)

             else:

                if np.isnan(vout[0,0]):

                   try:
                      vout[:,:] = self.ds_dictb[vname].sel(lat=slice(vdict['lat_start'], vdict['lat_end']),   \
                                               lon=slice(vdict['lon_start'], vdict['lon_end']),   \
                                               lev=slice(vdict['pres_start'], vdict['pres_end']), \
                                               time=slice(self.datef, self.datef),                \
                                               ens=slice(member+1, member+1)).squeeze().rename({'lon': 'longitude','lat': 'latitude','lev': 'isobaricInhPa'}).copy(deep=True)
                   except (OSError, RuntimeError):
                      print("Read Field Error; Sleeping before trying again")
                      time.sleep(240)
                      vout[:,:] = self.ds_dictb[vname].sel(lat=slice(vdict['lat_start'], vdict['lat_end']),   \
                                               lon=slice(vdict['lon_start'], vdict['lon_end']),   \
                                               lev=slice(vdict['pres_start'], vdict['pres_end']), \
                                               time=slice(self.datef, self.datef),                \
                                               ens=slice(member+1, member+1)).squeeze().rename({'lon': 'longitude','lat': 'latitude','lev': 'isobaricInhPa'}).copy(deep=True)


          #  Read the only level if it is a single level variable
          else:

             aout = self.ds_dict[vname].sel(ens=slice(member+1, member+1)).squeeze().rename({'lon': 'longitude','lat': 'latitude'})

             #  Rename the lat and lon arrays to common convention
             try:
                vout  = self.ds_dict[vname].sel(lat=slice(vdict['lat_start'], vdict['lat_end']), \
                                                lon=slice(vdict['lon_start'], vdict['lon_end']), \
                                                time=slice(self.datef, self.datef),              \
                                                ens=slice(member+1, member+1)).squeeze().rename({'lon': 'longitude','lat': 'latitude'}).copy(deep=True)
             except (OSError, RuntimeError):
                print("Read 2D Field Error; Sleeping before trying again") 
                time.sleep(240)
                vout  = self.ds_dict[vname].sel(lat=slice(vdict['lat_start'], vdict['lat_end']), \
                                                lon=slice(vdict['lon_start'], vdict['lon_end']), \
                                                time=slice(self.datef, self.datef),              \
                                                ens=slice(member+1, member+1)).squeeze().rename({'lon': 'longitude','lat': 'latitude'}).copy(deep=True)

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

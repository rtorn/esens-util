import os, sys
import argparse
import importlib
import configparser

sys.path.append('../esens-util')

def download_input(datea, paramfile):

  #  Read the configuration file
  conf = configparser.ConfigParser()
  conf.read(paramfile)

  #  Create new directory for input data
  conf['locations']['model_dir'] = '{0}/{1}'.format(conf['locations']['model_dir'],datea)
  conf['locations']['work_dir']  = conf['locations']['model_dir']
  os.makedirs(conf['locations']['model_dir'], exist_ok=True)

  #  Import the module that contains routines to Grib data specific to the model, download to input directory
  dpp = importlib.import_module(conf['model']['io_module'])

  os.chdir(conf['locations']['work_dir'])
  dpp.stage_grib_files(datea, conf)


if __name__ == '__main__':

  #  Read the initialization time and storm from the command line
  exp_parser = argparse.ArgumentParser()
  exp_parser.add_argument('--init',  action='store', type=str, required=True)
  exp_parser.add_argument('--param', action='store', type=str, default='example.parm')

  args = exp_parser.parse_args()

  download_input(args.init, args.param)

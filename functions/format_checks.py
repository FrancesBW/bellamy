#! /usr/bin/env python
from __future__ import print_function
import numpy as np
import os
#import textwrap
from scipy import interpolate
from scipy.interpolate import CloughTocher2DInterpolator
import astropy
from astropy import wcs
from astropy.io import fits
from astropy.io.votable import parse_single_table
from astropy.coordinates import SkyCoord, Angle, Latitude, Longitude, SkyOffsetFrame
from astropy.table import Table, hstack, join, vstack
import astropy.units as u
import sys
import glob
import argparse
#import psutil
import matplotlib.pyplot as plt
from matplotlib import gridspec
from copy import copy
import warnings
import logging
import logging.config

def check_column_exists(column,table,var,optional=False,length=None, variable=None,error=False):
        """
        Check that columns exist in a table and return that column. If interpolation is needed from the reference catalog return the interpolated column. If the column is an error column, return the appropriate interpolated error. Error will be raised if essential column cannot be found.
        
        :param column: the name of the column
        :param table: table being serached for the column name
        :param var: the variable pulled from the script needed for this function

        :param optional: If True and column cannot be found, a column of zeroes will be returned and error will not be raised, merely warning message will be logged (optional, default is False)
        :param length: Length of table, only required if optional is True (optional, default is None)
        :param variable: Name of data, logged in case of warning so that user is aware of what data is missing (optional, default is None)
        :param error: If True, the column is an error column for data that has been interpolated, the corresponding interpolated error needs to be returned

        :return: data from column if it exists, column of zeroes if data cannot be found or system exit is required column cannot be found.
        """
        #unpack var
        options, ref_min_freq, ref_max_freq, ref_p_or_s, num_freq, ref_freq,log=var
        log.info("Column name is: "+column)
        log.info("ref_p_or_s is: "+ref_p_or_s)
        print(ref_p_or_s=='suffix')
        
        #first try the column name as is
        try:
                test=table[column]
                return test
        #if it cannot be found, check if there is a frequency prefix, with or without the underscore
        except KeyError:
                log.info("Went to first except statement")
                if ref_p_or_s=='prefix':
                        log.info("Went down prefix path")
                        try:
                                if num_freq=='single':
                                        test=table[ref_freq+'_'+column]
                                        return test
                                else:
                                        #interpolate the data if it falls between frequencies
                                        test_1=table[ref_min_freq+'_'+column]
                                        test_2=table[ref_max_freq+'_'+column]
                                        if error==False:
                                                interpol=interpolate.interp1d([float(ref_min_freq),float(ref_max_freq)],[test_1,test_2],axis=0)
                                                return interpol(options.tar_freq)
                                        else:
                                                #if it is error, apply the appropriate error propagation given the interpolation
                                                return np.sqrt((test_1*((float(ref_max_freq)-options.tar_freq)/(float(ref_max_freq)-float(ref_min_freq))))**2+(test_2*((options.tar_freq-float(ref_min_freq))/(float(ref_max_freq)-float(ref_min_freq))))**2)
                        except KeyError:
                                try:
                                        if num_freq=='single':
                                                test=table[ref_freq+column]
                                                return test
                                        else:
                                                #interpolate the data if it falls between frequencies
                                                test_1=table[gref_min_freq+column]
                                                test_2=table[ref_max_freq+column]
                                                if error==False:
                                                        interpol=interpolate.interp1d([float(ref_min_freq),float(ref_max_freq)],[test_1,test_2],axis=0)
                                                        return interpol(options.tar_freq)
                                                else:
                                                        #if it is error, apply the appropriate error propagation given the interpolation
                                                        return np.sqrt((test_1*((float(ref_max_freq)-options.tar_freq)/(float(ref_max_freq)-float(ref_min_freq))))**2+(test_2*((options.tar_freq-float(ref_min_freq))/(float(ref_max_freq)-float(ref_min_freq))))**2)
                                except KeyError:
                                        #if no column matching either 3 names can be found, throw we can either return no data and continue or raise a fatal error
                                        if optional==True:
                                                #data is not required for algorithm to proceed, so just return zeroes and log a warning
                                                log.warning("No data found for the "+variable+" column. Margin of error will now be tighter and probability of matches will be lower.")
                                                return np.zeros(length)
                                        else:
                                                #data is required for algorithm to proceed, log error and exit the program
                                                log.error("Could not find column '{0}', '{1}_{0}' or '{1}{0}'. Please edit format file to reflect catalogue column names".format(column, ref_min_freq+'/'+ref_max_freq))
                                                sys.exit(1)
                #if it cannot be found, check if there is a frequency suffix, with or without the underscore
                elif ref_p_or_s=='suffix':
                        log.info("Went down suffix path")
                        try:
                                if num_freq=='single':
                                        log.info("Trying "+column+'_'+ref_freq)
                                        test=table[column+'_'+ref_freq]
                                        return test
                                else:
                                        #interpolate the data if it falls between frequencies
                                        log.info("Trying "+column+'_'+ref_min_freq)
                                        log.info("Trying "+column+'_'+ref_max_freq)
                                        test_1=table[column+'_'+ref_min_freq]
                                        test_2=table[column+'_'+ref_max_freq]
                                        if error==False:
                                                interpol=interpolate.interp1d([float(ref_min_freq),float(ref_max_freq)],[test_1,test_2],axis=0)
                                                return interpol(options.tar_freq)
                                        else:
                                                #if it is error, apply the appropriate error propagation given the interpolation
                                                return np.sqrt((test_1*((float(ref_max_freq)-options.tar_freq)/(float(ref_max_freq)-float(ref_min_freq))))**2+(test_2*((options.tar_freq-float(ref_min_freq))/(float(ref_max_freq)-float(ref_min_freq))))**2)
                        except KeyError:
                                try:
                                        if num_freq=='single':
                                                log.info("Trying "+column+ref_freq)
                                                test=table[column+ref_freq]
                                                return test
                                        else:
                                                #interpolate the data if it falls between frequencies
                                                log.info("Trying "+column+ref_min_freq)
                                                log.info("Trying "+column+ref_max_freq)
                                                test_1=table[column+ref_min_freq]
                                                test_2=table[column+ref_max_freq]
                                                if error==False:
                                                        interpol=interpolate.interp1d([float(ref_min_freq),float(ref_max_freq)],[test_1,test_2],axis=0)
                                                        return interpol(options.tar_freq)
                                                else:
                                                        #if it is error, apply the appropriate error propagation given the interpolation
                                                        return np.sqrt((test_1*((float(ref_max_freq)-options.tar_freq)/(float(ref_max_freq)-float(ref_min_freq))))**2+(test_2*((options.tar_freq-float(ref_min_freq))/(float(ref_max_freq)-float(ref_min_freq))))**2)
                                except KeyError:
                                        #if no column matching either 3 names can be found, we can either return no data and continue or raise a fatal error
                                        if optional==True:
                                                #data is not required for algorithm to proceed, so just return zeroes and log a warning
                                                log.warning("No data found for the "+variable+" column. Margin of error will now be tighter and probability of matches will be lower.")
                                                return np.zeros(length)
                                        else:
                                                #data is required for algorithm to proceed, log error and exit the program
                                                log.error("Could not find column '{0}', '{0}_{1}' or '{0}{1}'. Please edit format file to reflect catalogue column names".format(column, ref_min_freq+'/'+ref_max_freq))
                                                sys.exit(1)
                else:
                        #This path is executed if there is no prefix or suffix to the column name and the column name by itself cannot be found
                        if optional==True:
                                #data is not required for algorithm to proceed, so just return zeroes and log a warning
                                log.warning("No data found for the "+variable+" column. Margin of error will now be tighter and probability of matches will be lower.")
                                return np.zeros(length)
                        else:
                                #data is required for algorithm to proceed, log error and exit the program
                                log.error("Could not find column '{0}'. Please edit format file to reflect catalogue column names".format(column))
                                sys.exit(1)

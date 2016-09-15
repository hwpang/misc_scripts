#!/usr/bin/env python
# encoding: utf-8

"""
This script imports an individual RMG-Java kinetics library from a local directory and saves the output
kinetics library py file into a path of the user's choosing.  This library will be automatically
added to the 'libraryname' folder in the input/kinetics/libraries directory and can 
be used directly as an RMG-Py kinetics library.

usage:
importJavaKineticsLibrary.py [-h] INPUT LIBRARYNAME

positional arguments:
INPUT       the input path of the RMG-Java kinetics library directory
LIBRARYNAME      the libraryname for the RMG-Py format kinetics library

"""

import argparse
import os
from rmgpy.data.kinetics import KineticsFamily
from rmgpy import settings
          
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()
    parser.add_argument('inputPath', metavar='INPUT', type=str, nargs=1,
        help='the input path of the RMG-Java kinetics family directory')
    parser.add_argument('familyName', metavar='OUTPUT', type=str, nargs=1,
        help='the familyName for the RMG-Py format kinetics family')   
    
    args = parser.parse_args()
    inputPath = args.inputPath[0]
    familyName = args.familyName[0]
    
    family = KineticsFamily()
    family.loadOld(inputPath)
    
    try:
        os.makedirs(os.path.join(settings['database.directory'], 'kinetics', 'families', familyName))
    except:
        pass
    
    # Save in Py format
    family.save(os.path.join(settings['database.directory'], 'kinetics', 'families', familyName))

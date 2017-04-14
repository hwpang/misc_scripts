#!/usr/bin/env python
# encoding: utf-8

"""
Script to match a molecule to a node in a tree.
Takes SMILES string of molecule as input, along with desired kinetics family tree to descend.
"""

from rmgpy.data.kinetics import KineticsDatabase
from rmgpy.molecule import Molecule

def matchNode(speciesList,familyList,labelList=None):    
    
    if not isinstance(speciesList, list):
        speciesList = [speciesList]
    if not isinstance(familyList, list):
        familyList = [familyList]
    if labelList and not isinstance(labelList, list):
        labelList = [labelList]
    
    database = KineticsDatabase()
    database.load('/home/mjliu/Code/RMG-database/input/kinetics',families=familyList,libraries=[])
    
    matches = []
    
    for spec in speciesList:
        for family in familyList:
            if not isinstance(spec,Molecule):
                m = Molecule(SMILES=spec)
                if labelList:
                    for label in labelList:
                        print m.atoms
                        s = raw_input('Which atom to label as '+label+'? ')
                        try:
                            ind = int(s)
                        except ValueError:
                            s = raw_input('Invalid index, try again: ')
                            ind = int(s)
                        m.atoms[ind].label = label
            else:
                m = spec
            node = database.families[family].groups.descendTree(m,m.getLabeledAtoms(),strict=True)
            match = (spec,family,node)
            matches.append(match)
    
    return matches

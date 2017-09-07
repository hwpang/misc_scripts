#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os

import rmgpy
from rmgpy.data.rmg import RMGDatabase
from rmgpy.reaction import Reaction
from rmgpy.exceptions import UndeterminableKineticsError


source = ['Intra_R_Add_Endocyclic', 'Intra_R_Add_Exocyclic']
destination = ['Intra_R_Add_Polycyclic']
families = source + destination

databasePath = rmgpy.settings['database.directory']

database = RMGDatabase()
database.load(
    path=databasePath,
    thermoLibraries=[],
    reactionLibraries=[],
    seedMechanisms=[],
    kineticsFamilies=families,
    kineticsDepositories=['training']
    )


# Modified from addKineticsRulesFromTrainingSet

toMove = []


# Remove entries from source families
for label in source:
    family = database.kinetics.families[label]
    
    depository = family.getTrainingDepository()

    entries = depository.entries.values()
    entries.sort(key=lambda x: x.index)

    reverse_entries = []

    for entry in entries:
        for spc in entry.item.reactants + entry.item.products:
            if family.isMoleculeForbidden(spc.molecule[0]):
                print 'Found training reaction which is forbidden: {}'.format(entry)
                toMove.append(entry)
                depository.entries.pop(entry.index)

    for entry in entries:
        try:
            template = family.getReactionTemplate(entry.item)
        except UndeterminableKineticsError:
            reverse_entries.append(entry)
            continue

    for entry in reverse_entries:
        item = Reaction(reactants=entry.item.products, products=entry.item.reactants)

        try:
            template = family.getReactionTemplate(item)
        except UndeterminableKineticsError:
            print 'Found training reaction which does not fit: {}'.format(entry)
            toMove.append(entry)
            depository.entries.pop(entry.index)

# Add entries to destination families

print 'Moving {} training reactions total.'.format(len(toMove))

depository = database.kinetics.families[destination[0]].getTrainingDepository()

# start = 0
start = max(depository.entries.keys()) + 1

for i, entry in enumerate(toMove):
    valid = True
    for spc in entry.item.reactants + entry.item.products:
        if family.isMoleculeForbidden(spc.molecule[0]):
            print 'Training reaction is forbidden in destination: {}'.format(entry)
            valid = False

    try:
        template = family.getReactionTemplate(entry.item)
    except UndeterminableKineticsError:
        item = Reaction(reactants=entry.item.products, products=entry.item.reactants)

        try:
            template = family.getReactionTemplate(item)
        except UndeterminableKineticsError:
            print 'Training reaction does not fit in destination: {}'.format(entry)
            valid = False

    if valid:
        depository.entries[start + i] = entry

path = os.path.join(databasePath, 'kinetics/families')

for label in families:
    family = database.kinetics.families[label]
    for depository in family.depositories:
        family.saveDepository(depository, os.path.join(path, depository.label))


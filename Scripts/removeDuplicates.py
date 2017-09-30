#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This script removes duplicate reactions. Currently, keeps LibraryReactions
and removes RMG generated reactions.

To use, pass the paths of the Chemkin file on the command-line, e.g.

    $ python removeDuplicates.py /path/to/chem_annotated.inp

"""

import os
import argparse

import rmgpy.chemkin
from rmgpy.chemkin import loadChemkinFile, writeKineticsEntry
from rmgpy.rmg.model import CoreEdgeReactionModel
from rmgpy.data.kinetics.library import LibraryReaction

################################################################################

def main(chemkin):
    # Load Chemkin file
    reactionModel = CoreEdgeReactionModel()
    reactionModel.core.species, reactionModel.core.reactions = loadChemkinFile(chemkin)

    # Identify reactions to be removed
    reactionList = reactionModel.core.reactions
    duplicateReactionsToRemove = []
    for index1 in range(len(reactionList)):
        reaction1 = reactionList[index1]
        for index2 in range(index1 + 1, len(reactionList)):
            reaction2 = reactionList[index2]
            # Check if the two reactions are identical
            if (reaction1.reactants == reaction2.reactants and reaction1.products == reaction2.products) or (reaction1.reactants == reaction2.products and reaction1.products == reaction2.reactants):
                # Identify if exactly 1 of the reactions is a LibraryReaction
                if isinstance(reaction1, LibraryReaction) != isinstance(reaction2, LibraryReaction):
                    # Mark the non-LibraryReaction to be removed
                    if isinstance(reaction1, LibraryReaction):
                        duplicateReactionsToRemove.append(reaction2)
                    else:
                        duplicateReactionsToRemove.append(reaction1)

    # Remove the identified reactions
    for reaction in duplicateReactionsToRemove:
        reactionList.remove(reaction)

    # Write new Chemkin file
    path = 'chem_annotated_new.inp'
    speciesList = reactionModel.core.species + reactionModel.outputSpeciesList
    rxnList = reactionModel.core.reactions + reactionModel.outputReactionList
    with open(chemkin, 'rb') as old, open(path, 'wb') as new:
        # Copy species and reactions blocks to new Chemkin file
        line = old.readline()
        while 'REACTIONS' not in line.upper():
            new.write(line)
            line = old.readline()
        new.write('REACTIONS    KCAL/MOLE   MOLES\n\n')
        rmgpy.chemkin.__chemkin_reaction_count = 0
        for rxn in rxnList:
            new.write(writeKineticsEntry(rxn, speciesList=speciesList, verbose=True))
            new.write('\n')
        new.write('END\n\n')
    print "New Chemkin file contains {0} reactions.".format(rmgpy.chemkin.__chemkin_reaction_count)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('chemkin', metavar='CHEMKIN', type=str, nargs=1,
                        help='the Chemkin input file to visualize')

    args = parser.parse_args()

    chemkin = os.path.abspath(args.chemkin[0])

    main(chemkin)


#!/usr/bin/env python
# encoding: utf-8

"""
This script draws the entire group tree for the specified kinetics family,
including diagrams for each group definition.

Saves output as a (very large) pdf.
"""

import os
import argparse
import pydot
import rmgpy
from rmgpy.data.kinetics import KineticsDatabase
from rmgpy.molecule import Group

def convert(family):

    databasePath = rmgpy.settings['database.directory']
    familyPath = os.path.join(databasePath, 'kinetics', 'families')

    kinetics = KineticsDatabase()
    kinetics.loadFamilies(familyPath, families=[family])

    groups = []
    graph = pydot.Dot(graph_type='graph', dpi='52')
    for label, entry in kinetics.families[family].groups.entries.iteritems():
        if entry.index == -1:
            continue
        elif isinstance(entry.item, Group):
            imagePath = drawGroup(label, entry.item)
            nodeLabel = '<<TABLE border="0">\n<TR><TD>{1}</TD></TR>\n<TR><TD><IMG scale="True" SRC="{0}"/></TD></TR>\n</TABLE>>'.format(os.path.join(imagePath), label)
        else:
            nodeLabel = str(entry.item)
            nodeLabel = '"' + label + '\n' + nodeLabel + '"'

        groups.append(label)
        index1 = groups.index(label) + 1
        graph.add_node(pydot.Node(name=str(index1), label=nodeLabel, fontname="Helvetica", fontsize="16", shape="box"))
        if entry.parent is not None:
            index2 = groups.index(entry.parent.label) + 1
            graph.add_edge(pydot.Edge(src=str(index2), dst=str(index1), fontname="Helvetica", fontsize="16"))

    outputPath = family + '_groups.pdf'
    graph.write_pdf(outputPath, prog='dot')

def drawGroup(label, group):

    label = label.replace('/', '_')
    path = os.path.join(os.getcwd(), 'images', label + '.png')

    if os.path.exists(path):
        return path
    elif not os.path.exists('images'):
        os.makedirs('images')

    graph = pydot.Dot(graph_type='graph', dpi="52")
    for index, atom in enumerate(group.atoms):
        atomType = '{0!s} '.format(atom.label if atom.label != '' else '')
        atomType += ','.join([atomType.label for atomType in atom.atomType])
        atomType = '"' + atomType + '"'
        graph.add_node(pydot.Node(name=str(index + 1), label=atomType, fontname="Helvetica", fontsize="16"))
    for atom1 in group.atoms:
        for atom2, bond in atom1.bonds.iteritems():
            index1 = group.atoms.index(atom1)
            index2 = group.atoms.index(atom2)
            if index1 < index2:
                bondType = ','.join([order for order in bond.order])
                bondType = '"' + bondType + '"'
                graph.add_edge(pydot.Edge(src=str(index1 + 1), dst=str(index2 + 1), label=bondType, fontname="Helvetica", fontsize="16"))

    graph.write_png(path, prog='neato')

    return path


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('family', metavar='FAMILY', type=str, nargs=1,
        help='The name of the kinetics family')

    args = parser.parse_args()
    family = args.family[0]

    convert(family)

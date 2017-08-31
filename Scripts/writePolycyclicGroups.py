#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
from collections import OrderedDict

import rmgpy
import rmgpy.constants as constants
from rmgpy.data.base import Entry, LogicOr
from rmgpy.data.kinetics.family import KineticsFamily
from rmgpy.kinetics import Arrhenius, ArrheniusEP, ThirdBody, Lindemann, Troe, \
                           PDepArrhenius, MultiArrhenius, MultiPDepArrhenius, \
                           Chebyshev, KineticsData
from rmgpy.molecule.group import Group, GroupAtom, GroupBond

###############
# Load family #
###############

path = os.path.join(rmgpy.settings['database.directory'], 'kinetics', 'families')
label = 'Intra_R_Add_Polycyclic'

familyPath = os.path.join(path, label)

family = KineticsFamily(label=label)

local_context = {
        'KineticsData': KineticsData,
        'Arrhenius': Arrhenius,
        'ArrheniusEP': ArrheniusEP,
        'MultiArrhenius': MultiArrhenius,
        'MultiPDepArrhenius': MultiPDepArrhenius,
        'PDepArrhenius': PDepArrhenius,
        'Chebyshev': Chebyshev,
        'ThirdBody': ThirdBody,
        'Lindemann': Lindemann,
        'Troe': Troe,
        'R': constants.R,
    }

global_context = {}

family.load(path=familyPath, local_context=local_context, global_context=global_context, depositoryLabels=['!training'])

#######################
# Generate new groups #
#######################

# Parameters
s = [0, 1, 2, 3, 4, 5]  # Side chain
p = [1, 2, 3, 4]  # Attack position
r = [3, 4, 5, 6, 7, 8]  # Ring size
d = [0, 1]  # Direction - short or long

p_labels = {
    1: 'alpha',
    2: 'beta',
    3: 'gamma',
    4: 'delta',
}

d_labels = {
    0: 'short',
    1: 'long',
}

# Generate combinations
combos = []
logic_nodes = []
for s0 in s:
    logic_nodes.append('Rn{0}cx'.format(s0))
    for p0 in p:
        if ((s0 == 0 and p0 == 1) or
                (s0 == 1 and p0 > 1) or
                (s0 == 2 and p0 > 2) or
                (s0 == 3 and p0 > 3) or
                (s0 == 4 and p0 > 2) or
                (s0 == 5 and p0 > 1)):
            continue
        logic_nodes.append('Rn{0}cx_{1}'.format(s0, p_labels[p0]))
        for r0 in r:
            if p0 > r0/2.0 or (s0 == 0 and r0 == 3):
                continue
            elif p0 < r0 / 2.0:
                logic_nodes.append('Rn{0}c{1}_{2}'.format(s0, r0, p_labels[p0]))
                for d0 in d:
                    combos.append((s0, r0, p0, d0))
            else:  # p0 == r0/2.0
                combos.append((s0, r0, p0, None))

new_groups = {}

for s0, r0, p0, d0 in combos:
    # Generate label
    label = 'Rn{0}c{1}_{2}'.format(s0, r0, p_labels[p0])
    if d0 is not None:
        label += '_{0}'.format(d_labels[d0])
    print label

    # Create a new group, starting with the multiple bond
    new_group = Group().fromAdjacencyList("""
1 *3 R!H u0 r1 {2,[D,T,B]}
2 *2 R!H u0 r1 {1,[D,T,B]}
""")

    # Calculate number of atoms to add
    if d0:  # long
        a = s0 + r0 - p0 - 1  # Add to main chain
        b = p0 - 1  # To close other side of ring
    else:  # short
        a = s0 + p0 - 1  # Add to main chain
        b = r0 - p0 - 1  # To close other side of ring

    # Add on the side chain
    if a > 0:
        # Create atoms
        atoms = [GroupAtom(atomType=['R!H'], radicalElectrons=None, label='') for i in range(a)]
        # Set radical num to 1 for first atom
        atoms[0].radicalElectrons = [1]
        # Set ring properties for all atoms
        if s0 == 0:
            for atom in atoms:
                atom.props['inRing'] = True
        elif s0 > 1:
            for atom in atoms[0:s0 - 1]:
                atom.props['inRing'] = False
            for atom in atoms[s0:]:
                atom.props['inRing'] = True
        # Label atoms
        atoms[0].label = '*1'
        if a > 1:
            atoms[1].label = '*4'
        if a > 2:
            atoms[-1].label = '*5'
        if a > 3:
            for i, atom in enumerate(atoms[2:-1]):
                atom.label = '*{0}'.format(i + 6)
        # Add to group
        for atom in reversed(atoms):
            new_group.addAtom(atom)
        # Add bond to newly added atoms
        new_group.addBond(GroupBond(new_group.atoms[1], atoms[-1], order=[1, 2, 3, 1.5]))
        # Create bonds
        if a > 1:
            bonds = [GroupBond(atoms[i], atoms[i + 1], order=[1, 2, 3, 1.5]) for i in range(a - 1)]
            for bond in bonds:
                new_group.addBond(bond)

    # Add remaining atoms
    if b > 0:
        # Create atoms
        atoms = [GroupAtom(atomType=['R!H'], radicalElectrons=None, label='', props={'inRing': True}) for i in range(b)]
        # Add to group
        for atom in atoms:
            new_group.addAtom(atom)
        # Add bond to newly added atoms
        new_group.addBond(GroupBond(new_group.atoms[0], atoms[0], order=[1, 2, 3, 1.5]))
        # Close ring
        if d0:
            new_group.addBond(GroupBond(new_group.atoms[r0 - p0], atoms[-1], order=[1, 2, 3, 1.5]))
        else:
            new_group.addBond(GroupBond(new_group.atoms[p0], atoms[-1], order=[1, 2, 3, 1.5]))
        # Create bonds
        if b > 1:
            bonds = [GroupBond(atoms[i], atoms[i + 1], order=[1, 2, 3, 1.5]) for i in range(b - 1)]
            for bond in bonds:
                new_group.addBond(bond)
    else:  # b == 0
        # Close ring
        new_group.addBond(GroupBond(new_group.atoms[0], new_group.atoms[r0 - 1], order=[1, 2, 3, 1.5]))

    new_groups[label] = new_group

# Create tree structure
entry_dict = OrderedDict()

for s0 in s:
    label0 = 'Rn{0}cx'.format(s0)
    entry_dict[label0] = OrderedDict()

    for p0 in p:
        if ((s0 == 0 and p0 == 1) or
                (s0 == 1 and p0 > 1) or
                (s0 == 2 and p0 > 2) or
                (s0 == 3 and p0 > 3) or
                (s0 == 4 and p0 > 2) or
                (s0 == 5 and p0 > 1)):
            continue

        label1 = 'Rn{0}cx_{1}'.format(s0, p_labels[p0])
        entry_dict[label0][label1] = OrderedDict()

        for r0 in r:
            if p0 > r0/2.0 or (s0 == 0 and r0 == 3):
                continue

            elif p0 < r0 / 2.0:
                label2 = 'Rn{0}c{1}_{2}'.format(s0, r0, p_labels[p0])
                entry_dict[label0][label1][label2] = OrderedDict()

                for d0 in d:
                    label3 = 'Rn{0}c{1}_{2}_{3}'.format(s0, r0, p_labels[p0], d_labels[d0])
                    entry_dict[label0][label1][label2][label3] = new_groups[label3]

            else:  # p0 == r0/2.0
                label2 = 'Rn{0}c{1}_{2}'.format(s0, r0, p_labels[p0])
                entry_dict[label0][label1][label2] = new_groups[label2]


def get_children(top, item):

    entry0 = Entry(label=top)
    new_entries[top] = entry0

    children = []
    for key, value in item.iteritems():
        if isinstance(value, Group):
            entry = Entry(label=key, item=value)
            children.append(entry)
            new_entries[key] = entry
        elif isinstance(value, OrderedDict):
            children.append(get_children(key, value))

    entry0.item = LogicOr([child.label for child in children], invert=False)
    entry0.children = children

    for child in children:
        child.parent = entry0

    return entry0


new_entries = OrderedDict()
get_children('Rn_cyclic', entry_dict)

for key, value in new_entries.iteritems():
    family.groups.entries[key] = value

family.groups.entries = new_entries

family.groups.top[0] = new_entries['Rn_cyclic']

family.saveGroups('/home/mjliu/Code/groups.py')

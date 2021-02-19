#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This module defines the contains modified functions from rmgpy.chemkin in order to load chemkin files that contain fall off rates correctly when species dictionary is not provided.
"""

import os
from rmgpy.chemkin import remove_comment_from_line, read_species_block, ChemkinError, \
                        _read_kinetics_reaction, read_reaction_comments, read_thermo_block, \
                        _process_duplicate_reactions
from rmgpy.transport import TransportData
import rmgpy.kinetics as _kinetics
import logging
import re

def read_kinetics_entry(entry, species_dict, Aunits, Eunits):
    """
    Read a kinetics `entry` for a single reaction as loaded from a Chemkin
    file. The associated mapping of labels to species `species_dict` should also
    be provided. Returns a :class:`Reaction` object with the reaction and its
    associated kinetics.
    """
    Afactor = {
        'cm^3/(mol*s)': 1.0e6,
        'cm^3/(molecule*s)': 1.0e6,
        'm^3/(mol*s)': 1.0,
        'm^3/(molecule*s)': 1.0,
    }[Aunits[2]]

    lines = entry.strip().splitlines()

    # The first line contains the reaction equation and a set of
    # modified Arrhenius parameters
    reaction, third_body, kinetics, k_units, k_low_units = _read_kinetics_reaction(
        line=lines[0], species_dict=species_dict, Aunits=Aunits, Eunits=Eunits)

    if len(lines) == 1 and not third_body:
        # If there's only one line then we know to use the high-P limit kinetics as-is
        reaction.kinetics = kinetics['arrhenius high']
    else:
        # There's more kinetics information to be read
        kinetics.update({
            'chebyshev coefficients': [],
            'efficiencies': {},
        })

        # Note that the subsequent lines could be in any order
        for line in lines[1:]:
            kinetics = _read_kinetics_line(
                line=line, reaction=reaction, species_dict=species_dict, Eunits=Eunits,
                kunits=k_units, klow_units=k_low_units,
                kinetics=kinetics)

        # Decide which kinetics to keep and store them on the reaction object
        # Only one of these should be true at a time!
        if 'chebyshev' in kinetics:
            chebyshev = kinetics['chebyshev']
            if chebyshev.Tmin is None or chebyshev.Tmax is None:
                raise ChemkinError('Missing TCHEB line for reaction {0}'.format(reaction))
            if chebyshev.Pmin is None or chebyshev.Pmax is None:
                raise ChemkinError('Missing PCHEB line for reaction {0}'.format(reaction))
            if len(kinetics['chebyshev coefficients']) != (chebyshev.degreeT * chebyshev.degreeP):
                raise ChemkinError('Wrong number of Chebyshev coefficients '
                    'for reaction {0}'.format(reaction))
            index = 0
            for t in range(chebyshev.degreeT):
                for p in range(chebyshev.degreeP):
                    chebyshev.coeffs.value_si[t, p] = kinetics[
                        'chebyshev coefficients'][index]
                    index += 1
            # Don't forget to convert the Chebyshev coefficients to SI units!
            # This assumes that s^-1, cm^3/mol*s, etc. are compulsory
            chebyshev.coeffs.value_si[0, 0] -= (len(reaction.reactants) - 1) * math.log10(Afactor)
            reaction.kinetics = chebyshev
        elif 'pressure-dependent arrhenius' in kinetics:
            pdep_arrhenius = kinetics['pressure-dependent arrhenius']
            # Check for duplicates and combine them to MultiArrhenius objects
            duplicates_to_remove = []
            duplicates_to_add = []
            for index1 in range(len(pdep_arrhenius)):
                reaction1 = pdep_arrhenius[index1]
                p1, kinetics1 = reaction1
                if reaction1 in duplicates_to_remove:
                    continue
                for index2 in range(index1 + 1, len(pdep_arrhenius)):
                    reaction2 = pdep_arrhenius[index2]
                    p2, kinetics2 = reaction2
                    if p1 == p2:
                        if reaction1 not in duplicates_to_remove:
                            new_kinetics = _kinetics.MultiArrhenius()
                            duplicates_to_add.append([p1, new_kinetics])
                            new_kinetics.arrhenius = [kinetics1]
                            duplicates_to_remove.append(reaction1)
                        new_kinetics.arrhenius.append(kinetics2)
                        duplicates_to_remove.append(reaction2)
            for item in duplicates_to_remove:
                pdep_arrhenius.remove(item)
            pdep_arrhenius.extend(duplicates_to_add)

            pdep_arrhenius = sorted(pdep_arrhenius, key=lambda reaction: reaction[0])  # sort by ascending pressures

            reaction.kinetics = _kinetics.PDepArrhenius(
                pressures=([P for P, arrh in pdep_arrhenius], "atm"),
                arrhenius=[arrh for P, arrh in pdep_arrhenius],
            )
        elif 'troe' in kinetics:
            troe = kinetics['troe']
            troe.arrheniusHigh = kinetics['arrhenius high']
            troe.arrheniusLow = kinetics['arrhenius low']
            troe.efficiencies = kinetics['efficiencies']
            reaction.kinetics = troe
        elif third_body:
            reaction.kinetics = _kinetics.ThirdBody(
                arrheniusLow=kinetics['arrhenius low'])
            reaction.kinetics.efficiencies = kinetics['efficiencies']
        elif 'arrhenius low' in kinetics:
            reaction.kinetics = _kinetics.Lindemann(
                arrheniusHigh=kinetics['arrhenius high'],
                arrheniusLow=kinetics['arrhenius low'])
            reaction.kinetics.efficiencies = kinetics['efficiencies']
        elif 'explicit reverse' in kinetics or reaction.duplicate:
            # it's a normal high-P reaction - the extra lines were only either REV (explicit reverse) or DUP (duplicate)
            reaction.kinetics = kinetics['arrhenius high']
        elif 'sticking coefficient' in kinetics:
            reaction.kinetics = kinetics['sticking coefficient']
        else:
            raise ChemkinError(
                'Unable to understand all additional information lines for reaction {0}.'.format(entry))

        # These things may *also* be true
        if 'sri' in kinetics:
            reaction.kinetics.comment += "Warning: SRI parameters from chemkin file ignored on import. "

        if 'explicit reverse' in kinetics:
            reaction.kinetics.comment += \
                "Chemkin file stated explicit reverse rate: {0}".format(kinetics['explicit reverse'])

    return reaction

def read_reactions_block(f, species_dict, read_comments=True):
    """
    Read a reactions block from a Chemkin file stream.
    
    This function can also read the ``reactions.txt`` and ``pdepreactions.txt``
    files from RMG-Java kinetics libraries, which have a similar syntax.
    """
    energy_units = 'cal/mol'
    molecule_units = 'moles'
    volume_units = 'cm3'
    time_units = 's'

    line = f.readline()
    found = False
    while line != '' and not found:

        line = remove_comment_from_line(line)[0]
        line = line.strip()
        tokens = line.split()

        if len(tokens) > 0 and tokens[0].upper() == 'REACTIONS':
            # Regular Chemkin file
            found = True
            for token in tokens[1:]:
                unit = token.lower()
                if unit in ['molecules', 'moles', 'mole', 'mol', 'molecule']:
                    molecule_units = unit
                elif unit in ['kcal/mole', 'kcal/mol', 'cal/mole', 'cal/mol', 'kj/mole', 'kj/mol', 'j/mole', 'j/mol',
                              'kelvins']:
                    energy_units = unit
                else:
                    raise ChemkinError('Unknown unit type "{0}"'.format(unit))

        elif len(tokens) > 0 and tokens[0].lower() == 'unit:':
            # RMG-Java kinetics library file
            warnings.warn("The RMG-Java kinetic library files are"
                          " no longer supported and may be"
                          " removed in version 2.3.", DeprecationWarning)
            found = True
            while 'reactions:' not in line.lower():
                line = f.readline()
                line = remove_comment_from_line(line)[0]
                line = line.strip()

                if 'A:' in line or 'E:' in line:
                    units = line.split()[1]
                    if 'A:' in line:
                        molecule_units, volume_units, time_units = units.lower().split(
                            '/')  # Assume this is a 3-tuple: moles or molecules, volume, time
                    elif 'E:' in line:
                        energy_units = units.lower()
        else:
            line = f.readline()

    if not found:
        raise ChemkinError('Invalid reaction block.')

    # Check that the units are valid
    assert molecule_units in ['molecules', 'moles', 'mole', 'mol', 'molecule']
    assert volume_units in ['cm3', 'm3']
    assert time_units in ['s']
    assert energy_units in ['kcal/mole', 'kcal/mol', 'cal/mole', 'cal/mol', 'kj/mole', 'kj/mol', 'j/mole', 'j/mol',
                           'kelvins']

    # Homogenize units
    if molecule_units == 'molecules':
        molecule_units = 'molecule'
    elif molecule_units == 'moles' or molecule_units == 'mole':
        molecule_units = 'mol'
    volume_units = {'cm3': 'cm', 'm3': 'm'}[volume_units]
    if energy_units == 'kcal/mole':
        energy_units = 'kcal/mol'
    elif energy_units == 'cal/mole':
        energy_units = 'cal/mol'
    elif energy_units == 'kj/mole':
        energy_units = 'kj/mol'
    elif energy_units == 'j/mole':
        energy_units = 'j/mol'
    elif energy_units == 'kelvins':
        energy_units = 'K'
    energy_units = energy_units.replace('j/mol', 'J/mol')

    # Set up kinetics units
    Aunits = [
        '',  # Zeroth-order
        's^-1'.format(time_units),  # First-order
        '{0}^3/({1}*{2})'.format(volume_units, molecule_units, time_units),  # Second-order
        '{0}^6/({1}^2*{2})'.format(volume_units, molecule_units, time_units),  # Third-order
        '{0}^9/({1}^3*{2})'.format(volume_units, molecule_units, time_units),  # Fourth-order
    ]
    Eunits = energy_units

    kinetics_list = []
    comments_list = []
    kinetics = ''
    comments = ''

    line = f.readline()
    while line != '':

        line_starts_with_comment = line.lstrip().startswith('!') or line.lstrip().startswith('//')
        line, comment = remove_comment_from_line(line)
        line = line.strip()
        comment = comment.strip()

        if 'end' in line or 'END' in line:
            break

        # if 'rev' in line or 'REV' in line:
        # can no longer name reactants rev...
        #    line = f.readline()
        #    continue  # need to re-do the comment stripping!

        if '=' in line and not line_starts_with_comment:
            # Finish previous record
            kinetics_list.append(kinetics)
            comments_list.append(comments)
            kinetics = ''
            comments = ''

        if line: kinetics += line + '\n'
        if comment: comments += comment + '\n'

        line = f.readline()

    # Don't forget the last reaction!
    if kinetics.strip() != '':
        kinetics_list.append(kinetics)
        comments_list.append(comments)

    if len(kinetics_list) == 0 and len(comments_list) == 0:
        # No reactions found
        pass
    elif kinetics_list[0] == '' and comments_list[-1] == '':
        # True for Chemkin files generated from RMG-Py
        kinetics_list.pop(0)
        comments_list.pop(-1)
    elif kinetics_list[0] == '' and comments_list[0] == '':
        # True for Chemkin files generated from RMG-Java
        warnings.warn("RMG-Java loading is no longer supported and may be"
                      " removed in version 2.3.", DeprecationWarning)
        kinetics_list.pop(0)
        comments_list.pop(0)
    else:
        # In reality, comments can occur anywhere in the Chemkin
        # file (e.g. either or both of before and after the
        # reaction equation)
        # If we can't tell what semantics we are using, then just
        # throw the comments away
        # (This is better than failing to load the Chemkin file at
        # all, which would likely occur otherwise)
        if kinetics_list[0] == '':
            kinetics_list.pop(0)
        if len(kinetics_list) != len(comments_list):
            logging.warning("Discarding comments from Chemkin file because not sure which reaction they apply to")
            comments_list = ['' for kinetics in kinetics_list]

    reaction_list = []
    for kinetics, comments in zip(kinetics_list, comments_list):
        try:
            reaction = read_kinetics_entry(kinetics, species_dict, Aunits, Eunits)
            reaction = read_reaction_comments(reaction, comments, read=read_comments)
        except ChemkinError as e:
            if "Skip reaction!" in str(e):
                logging.warning("Skipping the reaction {0!r}".format(kinetics))
                continue
            else:
                raise
        reaction_list.append(reaction)

    return reaction_list

def load_chemkin_file(path, dictionary_path=None, transport_path=None, read_comments=True, thermo_path=None,
                      use_chemkin_names=False, check_duplicates=True):
    """
    Load a Chemkin input file located at `path` on disk to `path`, returning lists of the species
    and reactions in the Chemkin file. The 'thermo_path' point to a separate thermo file, or, if 'None' is
    specified, the function will look for the thermo database within the chemkin mechanism file
    """
    species_list = []
    species_dict = {}
    species_aliases = {}
    reaction_list = []

    # If the dictionary path is given, then read it and generate Molecule objects
    # You need to append an additional adjacency list for nonreactive species, such
    # as N2, or else the species objects will not store any structures for the final
    # HTML output.
    if dictionary_path:
        species_dict = load_species_dictionary(dictionary_path)

    with open(path, 'r') as f:
        previous_line = f.tell()
        line0 = f.readline()
        while line0 != '':
            line = remove_comment_from_line(line0)[0]
            line = line.strip()

            if 'SPECIES' in line.upper():
                # Unread the line (we'll re-read it in readReactionBlock())
                f.seek(previous_line)
                read_species_block(f, species_dict, species_aliases, species_list)

            elif 'THERM' in line.upper() and thermo_path is None:
                # Skip this if a thermo file is specified
                # Unread the line (we'll re-read it in read_thermo_block())
                f.seek(previous_line)
                formula_dict = read_thermo_block(f, species_dict)

            elif 'REACTIONS' in line.upper():
                # Reactions section
                # Unread the line (we'll re-read it in readReactionBlock())
                f.seek(previous_line)
                reaction_list = read_reactions_block(f, species_dict, read_comments=read_comments)

            previous_line = f.tell()
            line0 = f.readline()

    # Read in the thermo data from the thermo file        
    if thermo_path:
        with open(thermo_path, 'r') as f:
            line0 = f.readline()
            while line0 != '':
                line = remove_comment_from_line(line0)[0]
                line = line.strip()
                if 'THERM' in line.upper():
                    f.seek(f.tell() - len(line0), os.SEEK_SET)
                    formula_dict = read_thermo_block(f, species_dict)
                    break
                line0 = f.readline()
    # Index the reactions now to have identical numbering as in Chemkin
    index = 0
    for reaction in reaction_list:
        index += 1
        reaction.index = index

    # Process duplicate reactions
    if check_duplicates:
        _process_duplicate_reactions(reaction_list)

    # If the transport path is given, then read it to obtain the transport
    # properties
    
#     print(species_dict)
    if transport_path:
        load_transport_file(transport_path, species_dict)

    if not use_chemkin_names:
        # Apply species aliases if known
        for spec in species_list:
            try:
                spec.label = species_aliases[spec.label]
            except KeyError:
                pass

    # Attempt to extract index from species label
    indexPattern = re.compile(r'\(\d+\)$')
    for spec in species_list:
        if indexPattern.search(spec.label):
            label, sep, index = spec.label[:-1].rpartition('(')
            spec.label = label
            spec.index = int(index)

    reaction_list.sort(key=lambda reaction: reaction.index)
    return species_list, reaction_list, formula_dicts

def _read_kinetics_line(line, reaction, species_dict, Eunits, kunits, klow_units, kinetics):
    """
    Parse the subsequent lines of of a Chemkin reaction entry.
    """
    case_preserved_tokens = line.split('/')
    line = line.upper()
    tokens = line.split('/')

    if 'DUP' in line:
        # Duplicate reaction
        reaction.duplicate = True

    elif 'LOW' in line:
        # Low-pressure-limit Arrhenius parameters
        tokens = tokens[1].split()
        kinetics['arrhenius low'] = _kinetics.Arrhenius(
            A=(float(tokens[0].strip()), klow_units),
            n=float(tokens[1].strip()),
            Ea=(float(tokens[2].strip()), Eunits),
            T0=(1, "K"),
        )

    elif 'HIGH' in line:
        # What we thought was high, was in fact low-pressure
        kinetics['arrhenius low'] = kinetics['arrhenius high']
        kinetics['arrhenius low'].A = (
            kinetics['arrhenius low'].A.value, klow_units)
        # High-pressure-limit Arrhenius parameters
        tokens = tokens[1].split()
        kinetics['arrhenius high'] = _kinetics.Arrhenius(
            A=(float(tokens[0].strip()), kunits),
            n=float(tokens[1].strip()),
            Ea=(float(tokens[2].strip()), Eunits),
            T0=(1, "K"),
        )

    elif 'TROE' in line:
        # Troe falloff parameters
        tokens = tokens[1].split()
        alpha = float(tokens[0].strip())
        T3 = float(tokens[1].strip())
        T1 = float(tokens[2].strip())
        try:
            T2 = float(tokens[3].strip())
        except (IndexError, ValueError):
            T2 = None

        kinetics['troe'] = _kinetics.Troe(
            alpha=alpha,
            T3=(T3, "K"),
            T1=(T1, "K"),
            T2=(T2, "K") if T2 is not None else None,
        )
    elif line.strip().startswith('SRI'):
        kinetics['sri'] = True
        """To define an SRI pressure-dependent reaction, in addition to the LOW or HIGH parameters, the
        keyword SRI followed by three or five parameters must be included in the following order: a, b,
        c, d, and e [Eq. (74)]. The fourth and fifth parameters are options. If only the first three are
        stated, then by default d = 1 and e = 0.
        """
        # see eg. http://www.dipic.unipd.it/faculty/canu/files/Comb/Docs/chemkinCK.pdf

    elif 'CHEB' in line:
        # Chebyshev parameters
        chebyshev = kinetics.get(
            'chebyshev', _kinetics.Chebyshev(kunits=kunits))
        kinetics['chebyshev'] = chebyshev
        tokens = [t.strip() for t in tokens]
        if 'TCHEB' in line:
            index = tokens.index('TCHEB')
            tokens2 = tokens[index + 1].split()
            chebyshev.Tmin = Quantity(float(tokens2[0].strip()), "K")
            chebyshev.Tmax = Quantity(float(tokens2[1].strip()), "K")
        if 'PCHEB' in line:
            index = tokens.index('PCHEB')
            tokens2 = tokens[index + 1].split()
            chebyshev.Pmin = Quantity(float(tokens2[0].strip()), "atm")
            chebyshev.Pmax = Quantity(float(tokens2[1].strip()), "atm")
        if 'TCHEB' in line or 'PCHEB' in line:
            pass
        elif chebyshev.degreeT == 0 or chebyshev.degreeP == 0:
            tokens2 = tokens[1].split()
            chebyshev.degreeT = int(float(tokens2[0].strip()))
            chebyshev.degreeP = int(float(tokens2[1].strip()))
            chebyshev.coeffs = np.zeros((chebyshev.degreeT, chebyshev.degreeP), np.float64)
            # There may be some coefficients on this first line
            kinetics['chebyshev coefficients'].extend(
                [float(t.strip()) for t in tokens2[2:]])
        else:
            tokens2 = tokens[1].split()
            kinetics['chebyshev coefficients'].extend(
                [float(t.strip()) for t in tokens2])

    elif 'PLOG' in line:
        pdep_arrhenius = kinetics.get('pressure-dependent arrhenius', [])
        kinetics['pressure-dependent arrhenius'] = pdep_arrhenius
        tokens = tokens[1].split()
        pdep_arrhenius.append([float(tokens[0].strip()),
                               _kinetics.Arrhenius(
                                   A=(float(tokens[1].strip()), kunits),
                                   n=float(tokens[2].strip()),
                                   Ea=(float(tokens[3].strip()), Eunits),
                                   T0=(1, "K"),
                               )])
    elif tokens[0].startswith('REV'):
        reverse_A = float(tokens[1].split()[0])
        kinetics['explicit reverse'] = line.strip()
        if reverse_A == 0:
            logging.info("Reverse rate is 0 so making irreversible for reaction {0}".format(reaction))
            reaction.reversible = False
        else:
            logging.info("Ignoring explicit reverse rate for reaction {0}".format(reaction))
    elif line.strip() == 'STICK':
        # Convert what we thought was Arrhenius into StickingCoefficient
        k = kinetics['arrhenius high']
        kinetics['sticking coefficient'] = _kinetics.StickingCoefficient(
            A=k.A.value,
            n=k.n,
            Ea=k.Ea,
            T0=k.T0,
        )
        del kinetics['arrhenius high']
    else:
        # Assume a list of collider efficiencies
        try:
            for collider, efficiency in zip(case_preserved_tokens[0::2], case_preserved_tokens[1::2]):
                try:
                    efficiency = float(efficiency.strip())
                except ValueError:
                    raise ChemkinError("{0!r} doesn't look like a collision efficiency for species {1} in "
                                       "line {2!r}".format(efficiency, collider.strip(), line))
                if collider.strip() in species_dict:
                    kinetics['efficiencies'][species_dict[collider.strip()]] = efficiency
                else:  # try it with capital letters? Not sure whose malformed chemkin files this is needed for.
                    kinetics['efficiencies'][species_dict[collider.strip().upper()]] = efficiency
        except IndexError:
            error_msg = 'Could not read collider efficiencies for reaction: {0}.\n'.format(reaction)
            error_msg += 'The following line was parsed incorrectly:\n{0}'.format(line)
            error_msg += "\n(Case-preserved tokens: {0!r} )".format(case_preserved_tokens)
            raise ChemkinError(error_msg)
    return kinetics

def load_transport_file(path, species_dict):
    """
    Load a Chemkin transport properties file located at `path` and store the
    properties on the species in `species_dict`.
    """
    with open(path, 'r') as f:
        for line0 in f:
            line, comment = remove_comment_from_line(line0)
            line = line.strip()
            if line != '':
                # This line contains an entry, so parse it
                label = line[0:16].strip()
                data = line[16:].split()
                if label in species_dict:
                    species = species_dict[label]
                    species.transport_data = TransportData(
                        shapeIndex=int(data[0]),
                        sigma=(float(data[2]), 'angstrom'),
                        epsilon=(float(data[1]), 'K'),
                        dipoleMoment=(float(data[3]), 'De'),
                        polarizability=(float(data[4]), 'angstrom^3'),
                        rotrelaxcollnum=float(data[5]),
                        comment=comment.strip(),
                    )

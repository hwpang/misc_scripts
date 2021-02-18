#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This module defines the Model class, which contains a set of useful methods for
working with an RMG generated chemical mechanism.
"""

from __future__ import print_function
from collections import OrderedDict
from base64 import b64encode

from IPython.display import display, HTML
import matplotlib.pyplot as plt
import numpy as np

from rmgpy.chemkin import load_chemkin_file, load_species_dictionary, get_species_identifier
from rmgpy.molecule import Molecule
from rmgpy.species import Species


class Model(object):
    """
    Container class for an RMG model.

    Can be instantiated with a species dictionary only or a chemkin file
    plus a species dictionary.

    Provides a number of methods to search for species or reactions along
    with the ability to display them in an IPython notebook.
    """

    def __init__(self, chem_path=None, dict_path=None, update_labels=None):
        self.chem_path = chem_path
        self.dict_path = dict_path
        self.update_labels = update_labels
        self.chemkin_label_dict = dict()
        self.rmg_label_dict = dict()
        self.label_index_dict = dict()

        self.species = None
        self.reactions = None
        self.species_dict = None
        self.formula_dict = None

        if self.chem_path is not None:
            self.load_model()
        elif self.dict_path is not None:
            self.load_species()

    @staticmethod
    def display_species_html(species_list, data=False):
        """Make an html table of species labels and drawings"""
        html = ['<table>']
        html += ['<tr>'
                 '<th colspan="{0}" style="text-align:left;">Total {1} Species</th>'
                 '</tr>'.format(2, len(species_list))]
        for species in species_list:
            if species is not None:
                html += ['<tr><td colspan="{0}" style="text-align:center;">{1}</td>'.format(1, species.label)]
                html += ['<td colspan="{0}" style="text-align:center;">'
                         '<img style="display:inline;" src="data:image/png;base64,{1}">'
                         '</td></tr>'.format(1, b64encode(species._repr_png_()).decode('utf-8'))]
                if data:
                    html += ['<tr><td colspan="{0}" style="text-align:left;">{1}</td></tr>'.format(2, species.thermo)]
        html += ['</table>']

        display(HTML(''.join(html)))

    @staticmethod
    def display_reaction_html(reaction_list, data=False):
        """Make an html table of reaction labels and drawings"""
        html = ['<table>']
        html += ['<tr>'
                 '<th colspan="{0}" style="text-align:left;">Total {1} Reactions</th>'
                 '</tr>'.format(2, len(reaction_list))]
        for reaction in reaction_list:
            html += ['<tr><td colspan="{0}" style="text-align:center;">{1}</td>'.format(1, reaction.label)]
            html += ['<td colspan="{0}" style="text-align:center;">'
                     '<img style="display:inline;" src="data:image/png;base64,{1}">'
                     '</td></tr>'.format(1, b64encode(reaction._repr_png_()).decode('utf-8'))]
            if data:
                html += ['<tr><td colspan="{0}" style="text-align:left;">{1}</td></tr>'.format(2, reaction.kinetics)]
        html += ['</table>']

        display(HTML(''.join(html)))

    def find_species(self, identifier):
        """Find species based on provided identifier"""
        if identifier in self.species:
            # This is a Species object
            return identifier
        elif identifier in self.species_dict:
            # This is a species label
            return self.species_dict[identifier]
        else:
            try:
                spc = Species(molecule=[Molecule().from_smiles(identifier)])
            except ValueError:
                try:
                    spc = Species(molecule=[Molecule().from_inchi(identifier)])
                except ValueError:
                    raise ValueError('Unable to interpret provided identifier {0} as '
                                     'species label, SMILES, or InChI.'.format(identifier))

            spc.generate_resonance_structures()

            for species in self.species_dict.values():
                if spc.is_isomorphic(species):
                    return species
            else:
                raise ValueError('Unable to find any species matching the identifier {0}.'.format(identifier))

    def get_formulas(self):
        """Return sorted list of all available chemical formulas"""
        return sorted(self.formula_dict.keys())

    def get_labels(self, identifiers):
        """Return list of labels based on the provided list of identifiers"""
        return [species.label for species in self.get_species(identifiers)]

    def get_reactions_by_species(self, targets, require_all=False):
        """
        Find reactions containing the target species.

        If ``require_all=False``, return reactions containing any of the specified species.
        if ``require_all=True``, return reactions containing all of the specified species.
        """
        if not isinstance(targets, list):
            targets = [targets]
        species_list = self.get_species(targets)

        reaction_list = []
        for reaction in self.reactions:
            any_flag = False
            all_flag = True
            for species in species_list:
                if species in reaction.reactants or species in reaction.products:
                    any_flag = True
                else:
                    all_flag = False
            if not require_all and any_flag or require_all and all_flag:
                reaction_list.append(reaction)

        return reaction_list

    def get_smiles(self, identifiers):
        """Return list of SMILES based on the provided list of identifiers"""
        return [species.molecule[0].to_smiles() for species in self.get_species(identifiers)]

    def get_species(self, identifiers):
        """
        Return list of species based on the provided list of identifiers.

        If chemical formulas are provided in the list of identifiers, all species
        with that chemical formula will be included in the output species list.
        """
        if not isinstance(identifiers, list):
            identifiers = [identifiers]

        species_list = []
        for identifier in identifiers:
            if identifier in self.formula_dict:
                species_list.extend(self.formula_dict[identifier])
            else:
                try:
                    species_list.append(self.find_species(identifier))
                except ValueError:
                    raise ValueError('Unable to find any species matching the identifier {0}. Tried interpreting '
                                     'as chemical formula, species label, SMILES, and InChI.'.format(identifier))

        return species_list

    def load_model(self):
        """Load model from Chemkin file"""
        self.species, self.reactions = load_chemkin_file(self.chem_path, self.dict_path)
        if self.update_labels:
            self.update_all_labels()
        self.update_species_dict()
        self.update_formula_dict()
        self.update_label_index_dict()

    def update_label_index_dict(self):
        for i, spc in enumerate(self.species):
            self.label_index_dict[spc.label] = i

    def load_species(self):
        """
        Slightly modified version of chemkin.loadSpeciesDictionary

        Key changes: use OrderedDict, don't check for inerts, don't generate resonance structures
        """
        self.species_dict = OrderedDict()
        with open(self.dict_path, 'r') as f:
            adjlist = ''
            for line in f:
                if line.strip() == '' and adjlist.strip() != '':
                    # Finish this adjacency list
                    species = Species().from_adjacency_list(adjlist)
                    label = species.label
                    self.species_dict[label] = species
                    adjlist = ''
                else:
                    adjlist += line
            else:  # reach end of file
                if adjlist.strip() != '':
                    species = Species().from_adjacency_list(adjlist)
                    label = species.label
                    self.species_dict[label] = species

        self.species = self.species_dict.values()
        self.update_formula_dict()

    @staticmethod
    def plot_kinetics(reaction, T=None, P=None):
        """Plot reaction kinetics"""
        if T is None:
            tlist = np.linspace(298, 2000, 20)
        elif isinstance(T, (float, int)):
            tlist = np.array([T])
        else:
            # Assume its a list or array
            tlist = T
        if P is None:
            plist = [1e5]
        elif isinstance(P, (float, int)):
            plist = [P]
        else:
            # Assume its a list or array
            plist = P

        # plt.figure()

        x = 1000 / tlist
        for p in plist:
            k = []
            for t in tlist:
                k.append(reaction.kinetics.get_rate_coefficient(t, p))

            y = np.log10(k)
            plt.plot(x, y, label='P={0} Pa'.format(p))

        plt.xlabel('1000/T (K)')
        plt.ylabel('log(k)')
        plt.title(reaction.label)
        plt.legend()
        plt.show()

    def show_reactions(self, item=None, data=False):
        """
        Display table of reactions. Can optionally provide any of the following:

          - single Species object to search for
          - list of labels or identifiers to search for
          - chemical formula to search for

        If ``kinetics=True``, will also display reaction kinetics.
        """
        if item is None:
            species_list = self.species
        else:
            species_list = self.get_species(item)

        reaction_list = self.get_reactions_by_species(species_list)

        self.display_reaction_html(reaction_list, data=data)

    def show_species(self, item=None, data=False):
        """
        Display table of species. Can optionally provide any of the following:

          - single Species object to display
          - list of labels or identifiers to display
          - chemical formula to display

        If ``thermo=True``, will also display species thermo.
        """
        if item is None:
            species_list = self.species
        else:
            species_list = self.get_species(item)

        self.display_species_html(species_list, data=data)

    def update_all_labels(self):
        """Update species and reaction labels to match Chemkin labels"""
        for reaction in self.reactions:
            reaction.label = reaction.to_chemkin(kinetics=False)
        for species in self.species:
            chemkin_label = get_species_identifier(species)
            self.chemkin_label_dict[chemkin_label] = species.label
            self.rmg_label_dict[species.label] = chemkin_label
            species.label = chemkin_label
            

    def update_formula_dict(self):
        """Sort species into dictionary by chemical formula"""
        self.formula_dict = {}
        for species in self.species_dict.values():
            formula = species.molecule[0].get_formula()
            self.formula_dict.setdefault(formula, []).append(species)

    def update_species_dict(self):
        """Update self.species_dict based on self.species"""
        self.species_dict = OrderedDict([(species.label, species) for species in self.species])

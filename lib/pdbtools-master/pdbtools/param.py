#!/usr/bin/env python

# Copyright 2007, Michael J. Harms
# This program is distributed under General Public License v. 3.  See the file
# COPYING for a copy of the license.

__description__ = \
"""
pdb_param.py:

Calculates some molecular parameters given the sequence of the protein in a pdb
file.
    0) The molecular weight of the protein
    1) The pI of the protein assuming model compound pKa values
    2) The fraction titratable, acidic, basic
    3) The charge on the molecule assuming all ASP, GLU, HIS, LYS, and ARG are
       ionized
"""

__author__ = "Michael J. Harms"
__date__ = "080125"


import os, sys
from . import seq
from . import data
from .data.common import *

class PdbParamError(Exception):
    """
    General error class for this module.
    """
    pass

def calcChargeState(count_dict,pH):
    """
    Calculates the charge state of a protein at a particular pH.
    """

    charge_state = 0.
    for aa in list(count_dict.keys()):
        pKa = PKA_DICT[aa]
        charge = CHARGE_DICT[aa]
        charge_state += count_dict[aa]*charge/(1 + 10**(charge*(pH-pKa)))

    return charge_state


def pdbPi(aa_list,initial_pH_step=2.0,cutoff=0.001):
    """
    Calculate the PI of a sequence of amino acids.
    """

    # Count the number of titratable amino acids in a sequence
    count_dict = dict([(k,0) for k in list(PKA_DICT.keys())])
    for k in list(count_dict.keys()):
        count_dict[k] = len([aa for aa in aa_list if aa == k])
    count_dict["NTERM"] = 1
    count_dict["CTERM"] = 1

    # Determine the pI using a simple convergence algorithm
    #
    #  |<---------------no
    #  |                 |
    #  |             ------------
    # pH += step --> |Q(pH) < 0? | <-----------------|
    #                ------------                    |
    #                    |                           |
    #                   yes --> step *= 0.5 --> (pH -= step)

    step = initial_pH_step
    pH = 0
    while step > cutoff:
        if calcChargeState(count_dict,pH) > 0:
            pH += step
        else:
            step *= 0.5
            pH -= step

    return pH

def calcMW(aa_list):
    """
    Calculate the molecular weight of a sequence.
    """

    # Convert list to molecular weight
    try:
        mw = sum([MW_DICT[aa] for aa in aa_list])
    except KeyError:
        err = "Sequence contains non-standard amino acids!"
        raise PdbParamError(err)

    # Subtract the molecular for every bond in the protein
    mw = mw - (len(aa_list) - 1)*MW_H2O

    return mw


def pdbParam(pdb,chain="all",use_atoms=False):
    """
    Calculate some general electrostatic properties of a protein.
    """

    # Use pdb_seq to grab the sequence of the protein(s)
    chain_dict, seq_type = pdb_seq.pdbSeq(pdb,use_atoms)

    # Convert chain dictionary to list of amino acids
    if chain == "all":
        aa_list = []
        for c in list(chain_dict.keys()):
            aa_list.extend(chain_dict[c])
    else:
        if chain in list(chain_dict.keys()):
            aa_list = chain_dict[chain]
        else:
            err = "Chain \"%s\" is not in pdb file!" % chain
            raise PdbParamError(err)

    # Count number of each type of group
    count_dict = dict([(k,0) for k in list(MW_DICT.keys())])
    for k in list(count_dict.keys()):
        count_dict[k] = len([aa for aa in aa_list if aa == k])

    # Calculate molecular weight and pI
    mw = calcMW(aa_list)
    pI = pdbPi(aa_list)


    return count_dict, mw, pI, seq_type

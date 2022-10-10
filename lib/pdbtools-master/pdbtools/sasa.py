#!/usr/bin/env python

# Copyright 2007, Michael J. Harms
# This program is distributed under General Public License v. 3.  See the file
# COPYING for a copy of the license.

__description__ = \
"""
pdb_sasa.py

Calculate the solvent accessible surface area of every atom in a protein using
NACCESS.
"""
__author__ = "Michael J. Harms"
__date__ = "080129"

import os, sys, copy, shutil

global standards

class PdbSasaError(Exception):
    """
    General error class for this module.
    """

    pass

def runNaccess(pdb_file,probe=1.4,z_sample=0.05,vdw_file=None,keep_temp=False):
    """
    Wrapper to run naccess.
    """

    curr_dir = os.getcwd()
    tmp_dir = "%s_sasa-tmp" % pdb_file[:-4]

    try:
        os.mkdir(tmp_dir)
    except OSError as xxx_todo_changeme:
        (errno,errstr) = xxx_todo_changeme.args
        if errno == 17:
            pass
        else:
            err = "Problem creating temporary directory (%s)" % tmp_dir
            raise PdbSasaError(err)

    shutil.copy(pdb_file,tmp_dir)
    os.chdir(tmp_dir)

    # Actually run naccess
    args = ['naccess',pdb_file,"-p","%.2F" % probe,"-z","%.2F" % z_sample]
    if vdw_file != None:
        args.extend(["-r","%s" % vdw_file])
    status = os.spawnvp(os.P_WAIT,'naccess',args)

    if status != 0:
        err = "Naccess failed with error %i.\n" % status
        err += "Is naccess in your path?"
        raise PdbSasaError(err)

    # Read the .asa file
    root = pdb_file[:pdb_file.rfind(".")]
    f = open("%s.asa" % root,'r')
    lines = f.readlines()
    f.close()

    # Grab residue/atom names and calculated accessibility
    out = ["%s %s\n" % (l[13:26],l[54:62]) for l in lines if l[0:4] == "ATOM"]

    os.chdir("..")

    # Remove temporary files
    if not keep_temp:
        shutil.rmtree(tmp_dir)

    return out


def readStandards(standard_dir=None,probe_radius=1.4,z_sample=0.05,
                  vdw_file=None):
    """
    Calculate solvent accessibility of residues in blocked peptides.
    """

    # Determine where to look for the standard pdb files
    if standard_dir == None:
        script_dir = os.path.realpath(os.path.split(__file__)[0])
        standard_dir = os.path.join(script_dir,"pdb_data","peptides")

    # Move to the standards directory
    current_dir = os.getcwd()
    try:
        os.chdir(standard_dir)
    except OSError:
        err = "Standard dir \"%s\" does not exist!" % standard_dir
        raise PdbSasaError(err)

    # Create 2d dictionary, keying to residue type and atom within that residue
    standards = {}
    standard_list = [f for f in os.listdir('.') if f[-4:] == ".pdb"]
    for pdb in standard_list:

        residue_standard = runNaccess(pdb,probe_radius,z_sample,vdw_file)
        residue_standard = [l for l in residue_standard if l[10:13] == "  3"]

        residue_name = residue_standard[0][4:7]
        residue_sasa = dict([(l[0:3],float(l[13:])) for l in residue_standard])

        standards.update([(residue_name,copy.deepcopy(residue_sasa))])

    # Return to working directory
    os.chdir(current_dir)

    return standards


def pdbSASA(pdb_file,probe_radius=1.4,z_sample=0.05,vdw_file=None,
            keep_temp=False,standards=None):
    """
    Calculate the absolute and relative accessibility of every atom in a pdb
    file using NACCESS.
    """

    if standards == None:
        try:
            print("Reading standards")
            standards = readStandards(probe_radius=probe_radius,
                                      z_sample=z_sample,
                                      vdw_file=vdw_file)
        except PdbSasaError:
            print("Problem with standard calculation.  Only absolute")
            print("accessibility will be calculated.")
            standards = {}

    print(pdb_file)
    all_atoms = runNaccess(pdb_file,probe_radius,z_sample,vdw_file,keep_temp)

    residues = []
    for atom in all_atoms:
        if atom[4:13] not in residues:
            residues.append(atom[4:13])

    out = []
    for atom in all_atoms:

        absolute = float(atom[13:])
        try:
            absolute_peptide = standards[atom[4:7]][atom[:3]]
            fractional = "%10.3F" % (absolute/absolute_peptide)
        except (KeyError,ZeroDivisionError):
            fractional = "%10s" % "NA"

        atom_type = atom[:3]
        resid_type = atom[4:7]
        residue = "\"%s\"" % atom[7:13]
        out.append("%10s%10s%12s%10.3F%10s\n" % \
                   (atom_type,resid_type,residue,absolute,fractional))

    return out

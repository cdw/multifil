#!/usr/bin/env python
# encoding: utf-8
"""
hdf.py - allow the writing of sarcomere data to hdf5 formated files

Created by Dave Williams on 2011-10-27
"""

import numpy
import tables as hdf


## Attempt at HDF data saving
class Sarcomere(hdf.IsDescription):
    name            = hdf.StringCol(36)
    axialforce      = hdf.Float64Col()
    radialforce     = hdf.Float64Col(2)
    radialtension   = hdf.Float64Col()
    timestep        = hdf.Int32Col()
    frac_in_states  = hdf.Float64Col(3)
    hiding_line     = hdf.Float64Col()
    lattice_spacing = hdf.Float64Col()
    z_line          = hdf.Float64Col()
    class Thick(hdf.IsDescription):
        id            = hdf.Int32Col()
        axialforce    = hdf.Float64Col()
        radialforce   = hdf.Float64Col(2)
        radialtension = hdf.Float64Col()
        k             = hdf.Float64Col()
        rests         = hdf.Float64Col(60)
        bare_zone     = hdf.Float64Col()
        axial_locs    = hdf.Float64Col(60)
        axial_forces  = hdf.Float64Col(60)
        radial_forces = hdf.Float64Col((60, 2))
        class ThickFace(hdf.IsDescription):
            id           = hdf.Int32Col()
            thin_face    = hdf.Int32Col(2)
            radial_force = hdf.Int32Col()
            class XB(hdf.IsDescription):
                id          = hdf.Int32Col()
                state       = hdf.Int32Col()
                base        = hdf.Float64Col()
                rest        = hdf.Float64Col(2)
                bound_to    = hdf.Int32Col()
                energy      = hdf.Float64Col()
                axialforce  = hdf.Float64Col()
                radialforce = hdf.Float64Col()
    class Thin(hdf.IsDescription):
        id = hdf.Int32Col()
        axialforce    = hdf.Float64Col()
        radialforce   = hdf.Float64Col(2)
        k             = hdf.Float64Col()
        axial_locs    = hdf.Float64Col(90)
        axial_forces  = hdf.Float64Col(90)
        radial_forces = hdf.Float64Col((90, 2))
        class ThinFace(hdf.IsDescription):
            id = hdf.Int32Col()
            thick_face = hdf.Int32Col(2)
            radial_force = hdf.Int32Col()
            class BindingSite(hdf.IsDescription):
                id          = hdf.Int32Col()
                state       = hdf.Int32Col()
                base        = hdf.Float64Col()
                bound_to    = hdf.Int32Col()
                axialforce  = hdf.Float64Col()
                radialforce = hdf.Float64Col()

class HDF_file(object):
    def __init__(self, filename, filetitle='', overwrite=False):
        """ Open up the file to write or append to.

        Parameters:
            filename: the relative or absolute name of the hdf file
            filetitle: the title string to set on the root group
            overwrite: set true to not append to, but delete existing file
        """
        mode = 'a' if overwrite == False else 'w'
        self.file = hdf.openFile(filename, mode, filetitle)



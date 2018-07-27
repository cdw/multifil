#!/usr/bin/env python
# encoding: utf-8
"""
ti.py - A titin filament, with variable compliance

MORE ABOUT HOW TITIN WORKS

Created by Joe Powers and Dave Williams on 2017-02-17
"""


import numpy as np
import warnings
import numpy.random as random
random.seed() # Ensure proper seeding


class Titin:
    """This is all about the titin filament"""
    def __init__(self, parent_lattice, index, thick_face, thin_face,
        a=240, b=0.0045):
        """Initialize the titin filament.

        Parameters:
            parent_lattice: calling half-sarcomere instance
            index: which titin filament this is (0-23)
            thick_face: List of thick filament faces' numerical orientation (0-5)
            thin_face: List of thin filament faces' numerical orientation (0-2)

        Returns:
            None
        """
        # Name of the titin molecule
        self.index = index
        self.address = ('titin', index)
        # A connection to the parent lattice
        self.parent_lattice = parent_lattice
        # Location of the titin filament relative thick filament
        self.thick_face = thick_face
        # Link titin to the thick filament face
        self.thick_face.link_titin(self)
        # location of the titin filament relative to the thin filament
        self.thin_face = thin_face
        # Link titin to that face of the thin filament
        self.thin_face.link_titin(self)
        ## And now we declare titin properties that will be used to
        ## calculate forces
        self.rest = 120 #nm, no citation
        # Create the constants that determine force NOTE IMPROVE DOC
        self.a = a
        self.b = b

    def to_dict(self):
        """Create a JSON compatible representation of titin

        Usage example: json.dumps(titin.to_dict(), indent=1)

        Current output includes:
            a: force constant
            b: force constant
            address: largest to most local, indices for finding this
            rest: rest length
        """
        td = self.__dict__.copy()
        td.pop('index')
        td.pop('parent_lattice')
        td['thick_face'] = td['thick_face'].address
        td['thin_face'] = td['thin_face'].address
        return td

    def from_dict(self, td):
        """ Load values from a thin face dict. Values read in correspond to
        the current output documented in to_dict.
        """
        # Check for index mismatch
        read, current = tuple(td['address']), self.address
        assert read==current, "index mismatch at %s/%s"%(read, current)
        # Local keys
        self.a = td['a']
        self.b = td['b']
        self.rest = td['rest']
        self.thick_face = self.parent_lattice.resolve_address(
            td['thick_face'])
        self.thin_face = self.parent_lattice.resolve_address(
            td['thin_face'])

    def angle(self):
        """Caclulate the angle the titin makes relative to thick filament"""
        act_loc = self.thin_face.parent_thin.parent_lattice.z_line
        myo_loc = self.thick_face.get_axial_location(-1)
        ls = self.parent_lattice.lattice_spacing
        angle = np.arctan2(ls, act_loc-myo_loc)
        return angle

    def length(self):
        """Calculate the length of the titin filament"""
        act_loc = self.thin_face.parent_thin.parent_lattice.z_line
        myo_loc = self.thick_face.get_axial_location(-1)
        ls = self.parent_lattice.lattice_spacing
        length = np.sqrt( (act_loc-myo_loc)**2 + ls**2 )
        return length

    def stiffness(self):
        """Need instintanious stiffness at the current length to normalize
        force for settling at each timestep. We get this as the derivitive of
        force with respect to x. D[a*exp(b*x), x] = a*b*exp(b*x)
        """
        return self.force()*self.b

    def force(self):
        """Calculate the total force the titin filament exerts"""
        length = self.length()
        if length < self.rest:
            return 0 # titin buckles
        else:
            #Exponential Model
            x = length - self.rest
            return self.a*np.exp(self.b*x)
        #else:
            #Cubic model
            #x = length - self.rest
            #a = 3.25e-5
            #b = -5e-3
            #c = 9e-1
            #d = -7e-5
            #return a*(x**3) + b*(x**2) + c*x + d #Linke et al. PNAS 1998
        #else:
            # Sawtooth Model
            #extension = length - self.rest
            #a = 1.2e-6
            #b = 3
            #c = 0.5e-6#*random.rand()
            #d = 34
            #e = 3
            #return a*extension**b + c*np.fmod(extension,d)**e
        #else:
            #Hookean Model
            #return self.k * (self.length() - self.rest)

    def axialforce(self):
        """Return the total force the titin filament exerts on the
        thick filament's final node, (negate for the thin filament side)
        """
        return self.force() * np.cos(self.angle()) ## CHECK_JDP ##

    def radialforce(self):
        """Return the force in the radial direction (positive is compressive)
        TODO: The 'positive is compressive' part needs to be double checked
        """
        warnings.warn("Check radial force direction in titin")
        return self.force() * np.sin(self.angle())

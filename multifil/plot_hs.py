#!/usr/bin/env python
# encoding: utf-8
"""
plot_hs.py - Plot the half-sarcomere with mayavi

Created by Dave Williams on 2010-10-4
"""

import numpy as np
from enthought.mayavi import mlab

# Configure the graph
fig = mlab.figure(1, bgcolor=(0, 0, 0), size=(350, 350))
mlab.clf()
fig.scene.parallel_projection = True
mlab.view(-4.0, 84.5, 2423.2, (625.0, 21.4, -3.4))
fil_rad = 2
fil_seg = 12
fil_color = 'jet'
fil_lims = (-1.0, 1.0)

class plot_hs(object):
    
    def update_locs(self):
        # Get needed info from the half-sarcomere
        self.thick_xlocs = [t.axial for t in self.hs.thick]
        self.thin_xlocs = [t.axial for t in self.hs.thin]
        self.thick_s = [t.axialforce() for t in self.hs.thick]
        self.thin_s = [t.axialforce() for t in self.hs.thin]
        self.z_line = self.hs.z_line
        ls = self.hs.get_lattice_spacing()
        # Calculate y and z locations of fils
        ls_g = np.sqrt(3)/2 * ls
        ls_d = 0.5 * ls
        self.thick_yzlocs = [(0, 0), 
                             (0 + 2*ls_g, 0), 
                             (0 + ls_g,   0 - 3*ls_d), 
                             (0 + 3*ls_g, 0 - 3*ls_d)]
        act_a = lambda (y,z): (y - ls_g, z + ls_d)
        act_b = lambda (y,z): (y, z + 2*ls_d)
        act_c = lambda (y,z): (y + ls_g, z + ls_d)
        act_d = lambda (y,z): (y + ls_g, z - ls_d)
        act_e = lambda (y,z): (y, z - 2*ls_d)
        act_f = lambda (y,z): (y - ls_g, z - ls-d) 
        self.thin_yzlocs = [act_c(self.thick_yzlocs[1]),
                            act_b(self.thick_yzlocs[0]),
                            act_a(self.thick_yzlocs[1]),
                            act_b(self.thick_yzlocs[1]),
                            act_b(self.thick_yzlocs[3]),
                            act_c(self.thick_yzlocs[3]),
                            act_b(self.thick_yzlocs[2]),
                            act_c(self.thick_yzlocs[2])]
    
    def update_ends(self):
        """Update the effective forces at filament ends"""
        self.thick_end = [t.effective_axial_force() for t in self.hs.thick]
        self.thin_end = [t.effective_axial_force() for t in self.hs.thin]
    
    def update_bound(self):
        """Update which cross-bridges are bound and their states"""
        self.bound = []
        for thick in self.hs.thick:
            self.bound.append([])
            for face in thick.thick_faces:
                self.bound[-1].append([])
                for xb in face.xb:
                    if xb.get_numeric_state != 0:
                        self.bound[-1][-1].append(
                            (xb.face_index ,
                             xb.bound_to,
                             xb.get_numeric_state()))
    
    def __init__(self, hs):
        """Plot the half-sarcomere"""
        self.hs = hs
        # Trigger an update of location data
        self.update_locs()
        # Do initial plotting of the thick and thin fils
        self.thick_tubes = []
        for x, yz, s in zip(self.thick_xlocs, self.thick_yzlocs,
                            self.thick_s):
            y = np.repeat(yz[0], len(x))
            z = np.repeat(yz[1], len(x))
            self.thick_tubes.append(mlab.plot3d(x, y, z, s,
                tube_radius=fil_rad, tube_sides=fil_seg,
                colormap=fil_color, vmin=fil_lims[0], vmax=fil_lims[1]))
        self.thin_tubes = []
        for x, yz, s in zip(self.thin_xlocs, self.thin_yzlocs,
                            self.thin_s):
            y = np.repeat(yz[0], len(x))
            z = np.repeat(yz[1], len(x))
            self.thin_tubes.append(mlab.plot3d(x, y, z, s,
                tube_radius=0.6*fil_rad, tube_sides=fil_seg,
                colormap=fil_color, vmin=fil_lims[0], vmax=fil_lims[1]))
        # Plot the total force at the end of each filament
        self.update_ends()
        self.thick_end_cube = []
        for yz, s in zip(self.thick_yzlocs, self.thick_end):
            x = [0]
            y = [yz[0]]
            z = [yz[1]]
            s = [s]
            self.thick_end_cube.append(mlab.points3d(x,y,z,s,
                mode='cube', scale_mode='none', scale_factor=1.5*fil_rad,
                colormap=fil_color, vmin=-50, vmax=50))
        self.thin_end_cube = []
        for yz, s in zip(self.thin_yzlocs, self.thin_end):
            x = [self.z_line]
            y = [yz[0]]
            z = [yz[1]]
            s = [s]
            self.thin_end_cube.append(mlab.points3d(x,y,z,s,
                mode='cube', scale_mode='none', scale_factor=1.5*fil_rad,
                colormap=fil_color, vmin=-50, vmax=50))
        # Plot the cross-bridges
        self.update_bound()
        self.cross_bridges = []
        for fil, x, yz in zip(self.bound, self.thick_xlocs,
                              self.thick_yzlocs):
            y = [yz[0]]
            z = [yz[1]]
            for face in fil:
                # TODO: THIS IS WHERE I LEFT OFF ON IMPLEMENTING
                # CROSS-BRIDGE PLOTTING. THIS IS A BIT OF A STICKY FELLOW IN
                # THAT IT REQUIRES THE CROSS-BRIDGES BE EITHER ALL SHOWN
                # (WHICH IS VISUALLY CLUTTERED) OR DELETED AS NEEDED WITH 
                # EACH REDISPLAY (WHICH IS COMPLICATED). THE WAY TO GO IS
                # PROBABLY TO DELETE ALL WITH EACH REDISPLAY AND THEN PLOT 
                # ALL CROSS-BRIDGES ANEW EACH TIME.
                pass
    
    def update(self):
        """Update the visualization"""
        self.update_locs()
        self.update_ends()
        self.disable_rendering()
        for tube, x, yz, s in zip(self.thick_tubes, self.thick_xlocs,
                                  self.thick_yzlocs, self.thick_s):
            y = np.repeat(yz[0], len(x))
            z = np.repeat(yz[1], len(x))
            tube.mlab_source.set(x = x, y = y, z = z, scalars = s)
        for tube, x, yz, s in zip(self.thin_tubes, self.thin_xlocs,
                                  self.thin_yzlocs, self.thin_s):
            y = np.repeat(yz[0], len(x))
            z = np.repeat(yz[1], len(x))
            ts = tube.mlab_source
            ts.set(x = x, y = y, z = z, scalars = s)
        for cube, s in zip(self.thick_end_cube, self.thick_end):
            s = [s] 
            cube.mlab_source.set(scalars = s)
        for cube, s in zip(self.thin_end_cube, self.thin_end):
            s = [s] 
            cube.mlab_source.set(scalars = s)
        self.enable_rendering()
    
    def disable_rendering(self):
        """Kill rendering of the scene objects
        
        This makes things vastly faster if done when re-rendering
        things, as the whole scene will only be re-rendered once, 
        rather than as each """
        disable = lambda x: x.scene.set(disable_render=True)
        map(disable, self.thick_tubes)
        map(disable, self.thin_tubes)
        map(disable, self.thick_end_cube)
        map(disable, self.thin_end_cube)
    
    def enable_rendering(self):
        """Unkill rendering of the scene objects"""
        enable = lambda x: x.scene.set(disable_render=False)
        map(enable, self.thick_tubes)
        map(enable, self.thin_tubes)
        map(enable, self.thick_end_cube)
        map(enable, self.thin_end_cube)

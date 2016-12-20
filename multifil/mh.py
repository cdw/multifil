#!/usr/bin/env python
# encoding: utf-8
"""
mh.py - A single myosin head

Created by Dave Williams on 2010-01-04.
"""

import numpy.random as random
random.seed() # Ensure proper seeding
from numpy import pi, sqrt, log, radians
import math as m
import warnings

class Spring(object):
    """A generic spring, from which we make the myosin heads"""
    def __init__(self, config):
        ## Passed variables
        self.r_w = config['rest_weak']
        self.r_s = config['rest_strong']
        self.k_w = config['konstant_weak']
        self.k_s = config['konstant_strong']
        ## Diffusion governors
        # k_T = Boltzmann constant * temperature = (1.381E-23 J/K * 288 K)
        k_t = 1.381*10**-23 * 288 * 10**21 #10**21 converts J to pN*nM
        # Normalize: a factor used to normalize the PDF of the segment values
        self.normalize = sqrt(2*pi*k_t/self.k_w)
        self.stand_dev = sqrt(k_t/self.k_w) # of segment values
    
    def to_dict(self):
        """Create a JSON compatible representation of the spring """
        return self.__dict__.copy()
    
    def from_dict(self, sd):
        """ Load values from a spring dict. Values read in correspond 
        to the current output documented in to_dict.
        """
        self.r_w = sd['r_w']
        self.r_s = sd['r_s']
        self.k_w = sd['k_w']
        self.k_s = sd['k_s']
        self.normalize = sd['normalize']
        self.stand_dev = sd['stand_dev']
    
    def rest(self, state):
        """Return the rest value of the spring in state state 
        
        Takes:
            state: the state of the spring, ['free'|'loose'|'tight']
        Returns:    
            length/angle: rest length/angle of the spring in the given state
        """
        if state in ("free", "loose"):
            return self.r_w
        elif state == "tight":
            return self.r_s
        else:
            warnings.warn("Improper value for spring state")
    
    def constant(self, state):
        """Return the spring constant of the spring in state state 
        
        Takes:
            state: the state of the spring, ['free'|'loose'|'tight']
        Returns:    
            spring constant: for the spring in the given state
        """
        if state in ("free", "loose"):
            return self.k_w
        elif state == "tight":
            return self.k_s
        else:
            warnings.warn("Improper value for spring state")
    
    def energy(self, spring_val, state):
        """Given a current length/angle, return stored energy 
        
        Takes:
            spring_val: a spring length or angle
            state: a spring state, ['free'|'loose'|'tight']
        Returns:
            energy: the energy required to achieve the given value
        """
        if state in ("free", "loose"):
            return (0.5 * self.k_w * m.pow((spring_val-self.r_w), 2))
        elif state == "tight":
            return (0.5 * self.k_s * m.pow((spring_val-self.r_s), 2))
        else:
            warnings.warn("Improper value for spring state")
    
    def bop(self):
        """Bop for a new value, given an exponential energy dist
        
        A longer explanation is in singlexb/Crossbridge.py
        Takes:
            nothing: assumes the spring to be in the unbound state
        Returns:    
            spring_value: the length or angle of the spring after diffusion"""
        return (random.normal(self.r_w, self.stand_dev))
    

class SingleSpringHead(object):
    """A single-spring myosin head, as in days of yore"""
    def __init__(self):
        """Create the spring that makes up the head and set energy values"""
        self.state = "free"
        self.g = Spring({
            'rest_weak': 5,
            'rest_strong': 0,
            'konstant_weak': 5 / 3.976,          
            'konstant_strong': 5 / 3.976})
        # Free energy calculation helpers
        g_atp = 13 # In units of RT
        atp = 5  * 10**-3
        adp = 30 * 10**-6
        phos  = 3  * 10**-3
        self.deltaG = abs(-g_atp - log(atp / (adp * phos)))
        self.alpha = 0.28 
        self.eta = 0.68 
        # The time-step, master of all time
        self.timestep = 1 # ms
    
    def transition(self, bs):
        """Transition to a new state (or not)
        
        Takes:
            bs: relative Crown to Actin distance (x,y)
        Returns:
            boolean: transition that occurred (as string) or None
        """
        ## Transitions rates are checked against a random number
        check = random.rand()
        ## Check for transitions depending on the current state
        if self.state == "free":
            if self._r12(bs) > check:
                self.state = "loose"
                return '12'
        elif self.state == "loose":
            if self._r23(bs) > check:
                self.state = "tight"
                return '23'
            elif self._r21(bs) > check:
                self.state = "free"
                return '21'
        elif self.state == "tight":
            if self._r31(bs) > check:
                self.state = "free"
                return '31'
            elif self._r32(bs) > check:
                self.state = "loose"
                return '32'
        # Got this far? Than no transition occurred!
        return None
    
    def axialforce(self, tip_location):
        """Find the axial force a Head generates at a given location
        
        Takes:
            tip_location: relative Crown to Actin distance (x,y)
        Returns:
            f_x: the axial force generated by the Head
        """
        ## Get the Head length
        g_len = tip_location[0]
        ## Write all needed values to local variables
        g_s = self.g.rest(self.state)
        g_k = self.g.constant(self.state)
        ## Find and return force
        f_x = g_k * (g_len - g_s)
        return f_x
    
    def radialforce(self, tip_location):
        """Find the radial force a Head generates at a given location
        
        Takes:
            tip_location: relative Crown to Actin distance (x,y)
        Returns:
            f_y: the radial force generated by the Head
        """
        return 0.0
    
    def energy(self, tip_location, state=None):
        """Return the energy in the xb with the given parameters
        
        Takes:
            tip_location: relative Crown to Actin distance (x,y)
            state: kinetic state of the cross-bridge, ['free'|'loose'|'tight']
        Returns:
            xb_energy: the energy stored in the cross-bridge"""
        if state is None:
            state = self.state
        return self.g.energy(tip_location[0], state)
    
    def get_numeric_state(self):
        """Return the numeric state (0, 1, or 2) of the head"""
        lookup_state = {"free":0, "loose":1, "tight":2}
        return lookup_state[self.state]
    
    def _set_timestep(self, timestep):
        """Set the length of time step used to calculate transitions"""
        self.timestep = timestep
    
    def _r12(self, bs):
        """Binding rate, based on the distance from the Head tip to a Actin
        
        Takes:
            bs: relative Crown to Actin distance (x,y)
        Returns:
            probability: chance of binding occurring
        """
        ## Get needed values
        k_xb = self.g.constant("free")
        xb_0 = self.g.rest("free")
        A = 2000  # From Tanner, 2008 Pg 1209
        ## Calculate the binding probability
        rate = (A * sqrt(k_xb / (2 * pi)) * 
                m.exp(-.5 * k_xb * (bs[0] - xb_0)**2)) * self.timestep
        return float(rate)
    
    def _r21(self, bs):
        """The reverse transition, from loosely bound to unbound
        
        Takes:
            bs: relative Crown to Actin distance (x,y)
        Returns:
            rate: probability of transition occurring this timestep
        """
        ## The rate depends on the states' free energies 
        g_1 = self._free_energy(bs, "free")
        g_2 = self._free_energy(bs, "loose")
        ## Rate, as in pg 1209 of Tanner et al, 2007
        try:
            rate = self._r12(bs) / m.exp(g_1 - g_2)
        except ZeroDivisionError:
            rate = 1
        return float(rate)
    
    def _r23(self, bs):
        """Probability of becoming tightly bound if loosely bound
        
        Takes:
            bs: relative Crown to Actin distance (x,y)
        Returns:
            rate: probability of becoming tightly bound
        """
        ## Get other needed values
        k_xb = self.g.constant("loose")
        xb_0 = self.g.rest("loose")
        B = 100   # From Tanner, 2008 Pg 1209
        C = 1
        D = 1
        ## Rate taken from single cross-bridge work
        rate = (B / sqrt(k_xb) * (1 - m.tanh(C * sqrt(k_xb) *
                (bs[0] - xb_0))) + D) * self.timestep
        return float(rate)
    
    def _r32(self, bs):
        """The reverse transition, from tightly to loosely bound
        
        Takes:
            bs: relative Crown to Actin distance (x,y)
        Returns:
            rate: probability of becoming loosely bound
        """
        ## Governed as in self_r21
        g_2 = self._free_energy(bs, "loose")
        g_3 = self._free_energy(bs, "tight")
        try:
            rate = self._r23(bs) / m.exp(g_2 - g_3)
        except ZeroDivisionError:
            rate = 1
        return float(rate)
    
    def _r31(self, bs):
        """Probability of unbinding if tightly bound
        
        Takes:
            bs: relative Crown to Actin distance (x,y)
        Returns:
            rate: probability of detaching from the binding site
        """
        ## Get needed values
        k_xb = self.g.constant("tight")
        M = 3600 # From Tanner, 2008 Pg 1209
        N = 40
        P = 20
        ## Based on the energy in the tight state
        rate = (sqrt(k_xb) * (sqrt(M * (bs[0]-4.76)**2) - 
                N * (bs[0]-4.76)) + P) * self.timestep
        return float(rate)
    
    def _free_energy(self, tip_location, state):
        """Free energy of the Head
        
        Takes:
            tip_location: relative Crown to Actin distance (x,y)
            state: kinetic state of the cross-bridge, ['free'|'loose'|'tight']
        Returns:
            energy: free energy of the head in the given state
        """
        if state == "free":
            return 0
        elif state == "loose":
            k_xb = self.g.constant(state)
            xb_0 = self.g.rest(state)
            x = tip_location[0]
            return self.alpha * -self.deltaG + k_xb * (x - xb_0)**2
        elif state == "tight":
            k_xb = self.g.constant(state)
            x = tip_location[0]
            return self.eta * -self.deltaG + k_xb * x**2


class Head(object):
    """Head implements a single myosin head"""
    def __init__(self):
        """Create the springs that make up the head and set energy values
        Values are choosen for consistancy with single spring rest lengths
        and rest lattice spacings. More documentaion in the single spring
        code. All numerical values referenced are discussed in single
        crossbridge PLOS paper.
        """
        # Remember thine kinetic state
        self.state = "free"
        # Create the springs which make up the head
        self.c = Spring({   # the converter domain
            'rest_weak': radians(47.16),
            'rest_strong': radians(73.20),
            'konstant_weak': 40,
            'konstant_strong': 40})
        self.g = Spring({   # the globular domain
            'rest_weak': 19.93,
            'rest_strong': 16.47,
            'konstant_weak': 2,
            'konstant_strong': 2})
        # Free energy calculation helpers
        g_atp = 13 # In units of RT
        atp = 5  * 10**-3
        adp = 30 * 10**-6
        phos  = 3  * 10**-3
        deltaG = abs(-g_atp - log(atp / (adp * phos)))
        self.alphaDG = 0.28 * -deltaG
        self.etaDG = 0.68 * -deltaG
        # The time-step, master of all time
        self.timestep = 1 # ms
    
    def transition(self, bs, ap):
        """Transition to a new state (or not)
        
        Takes:
            bs: relative Crown to Actin distance (x,y)
            ap: Actin binding permissiveness, from 0 to 1
        Returns:
            boolean: transition that occurred (as string) or None
        """
        ## Transitions rates are checked against a random number
        check = random.rand()
        ## Check for transitions depending on the current state
        if self.state == "free":
            if self._prob(self._bind(bs))*ap > check:
                self.state = "loose"
                return '12'
        elif self.state == "loose":
            if self._prob(self._r23(bs)) > check:
                self.state = "tight"
                return '23'
            elif (1 - self._prob(self._r21(bs))) < check:
                self.state = "free"
                return '21'
        elif self.state == "tight":
            if self._prob(self._r31(bs)) > check:
                self.state = "free"
                return '31'
            elif (1 - self._prob(self._r32(bs))) < check:
                self.state = "loose"
                return '32'
        # Got this far? Than no transition occurred!
        return None
    
    def axialforce(self, tip_location):
        """Find the axial force a Head generates at a given location
        
        Takes:
            tip_location: relative Crown to Actin distance (x,y)
        Returns:
            f_x: the axial force generated by the Head
        """
        ## Get the Head length and angle
        (c_ang, g_len) = self._seg_values(tip_location)
        ## Write all needed values to local variables
        c_s = self.c.rest(self.state)
        g_s = self.g.rest(self.state)
        c_k = self.c.constant(self.state)
        g_k = self.g.constant(self.state)
        ## Find and return force
        f_x = (g_k * (g_len - g_s) * m.cos(c_ang) +
               1/g_len * c_k * (c_ang - c_s) * m.sin(c_ang))
        return f_x
     
    def radialforce(self, tip_location):
        """Find the radial force a Head generates at a given location
        
        Takes:
            tip_location: relative Crown to Actin distance (x,y)
        Returns:
            f_y: the radial force generated by the Head
        """
        ## Get the Head length and angle
        (c_ang, g_len) = self._seg_values(tip_location)
        ## Write all needed values to local variables
        c_s = self.c.rest(self.state)
        g_s = self.g.rest(self.state)
        c_k = self.c.constant(self.state)
        g_k = self.g.constant(self.state)
        ## Find and return force
        f_y = (g_k * (g_len - g_s) * m.sin(c_ang) +
               1/g_len * c_k * (c_ang - c_s) * m.cos(c_ang))
        return f_y
    
    def energy(self, tip_location, state=None):
        """Return the energy in the xb with the given parameters
        
        Takes:
            tip_location: relative Crown to Actin distance (x,y)
            state: kinetic state of the cross-bridge, ['free'|'loose'|'tight']
        Returns:
            xb_energy: the energy stored in the cross-bridge"""
        if state == None:
            state = self.state
        (ang, dist) = self._seg_values(tip_location)
        xb_energy = self.c.energy(ang, state) + self.g.energy(dist, state)
        return xb_energy
    
    def get_numeric_state(self):
        """Return the numeric state (0, 1, or 2) of the head"""
        lookup_state = {"free":0, "loose":1, "tight":2}
        return lookup_state[self.state]
    
    @property
    def timestep(self):
        return self._timestep
    
    @timestep.setter
    def timestep(self, timestep):
        """Set the length of time step used to calculate transitions"""
        self._timestep = timestep

    def _prob(self, rate):
        """Convert a rate to a probability, based on the current timestep
        length and the assumption that the rate is for a Poisson process.
        We are asking, what is the probability that at least one Poisson
        distributed value would occur during the timestep.

        Takes:
            rate: a per ms rate to convert to probability
        Returns:
            probability: the probability the event occurs during a timestep
                of length determined by self.timestep
        """
        return 1 - m.exp(-rate*self.timestep)
    
    def _bind(self, bs):
        """Bind (or don't) based on the distance from the Head tip to a Actin
        
        Takes:
            bs: relative Crown to Actin distance (x,y)
        Returns:
            probability: chance of binding occurring during a timestep
        """
        ## Flag indicates successful diffusion
        bop_right = False
        while bop_right is False:
            ## Bop the springs to get new values
            c_ang = self.c.bop()
            g_len = self.g.bop()
            ## Translate those values to an (x,y) position
            tip = (g_len * m.cos(c_ang), g_len * m.sin(c_ang))
            ## Only a bop that lands short of the thin fil is valid
            bop_right = bs[1] >= tip[1]
        ## Find the distance to the binding site
        distance = m.hypot(bs[0]-tip[0], bs[1]-tip[1])
    
        ## The binding rate is dependent on the exp of the dist
        # Rate = \tau * \exp^{-dist^2}
        rate = 72 * m.exp(-distance**2)
        ## Return the rate
        return rate
    
    def _r21(self, bs):
        """The reverse transition, from loosely bound to unbound
        
        This depends on the prob r12, the binding prob, which is given
        in a stochastic manner. Thus _p21 is returning not the prob of
        going from loosely bound to tightly bound, but the change that
        occurs in one particular timestep, the stochastic probability.
        Takes:
            bs: relative Crown to Actin distance (x,y)
        Returns:
            prob: probability of transition
        """
        ## The rate depends on the states' free energies
        unbound_free_energy = self._free_energy(bs, "free")
        loose_free_energy = self._free_energy(bs, "loose")
        ## Rate, as in pg 1209 of Tanner et al, 2007
        ## With added reduced-detachment factor, increases dwell time
        try:
            rate = self._bind(bs) / m.exp(
                unbound_free_energy - loose_free_energy)
        except ZeroDivisionError:
    
            rate = 1
        return float(rate)
    
    def _r23(self, bs):
        """Rate of becoming tightly bound if loosely bound
        
        Takes:
            bs: relative Crown to Actin distance (x,y)
        Returns:
            rate: per ms rate of becoming tightly bound
        """
        ## The transition rate depends on state energies
        loose_energy = self.energy(bs, "loose")
        tight_energy = self.energy(bs, "tight")
        ## Powerstroke rate, per ms
        rate = (0.6 * # reduce overall rate
                (1 +  # shift rate up to avoid negative rate
                m.tanh(6 + # move center of transition to right
                       0.2 * (loose_energy - tight_energy))))
        return float(rate)
    
    def _r32(self, bs):
        """The reverse transition, from tightly to loosely bound
        
        Takes:
            bs: relative Crown to Actin distance (x,y)
        Returns:
            rate: per ms rate of transition
        """
        ## Governed as in self_p21
        loose_free_energy = self._free_energy(bs, "loose")
        tight_free_energy = self._free_energy(bs, "tight")
        try:
            rate = self._r23(bs)/ m.exp(loose_free_energy - tight_free_energy)
        except ZeroDivisionError:
            rate = 1
        return float(rate)
    
    def _r31(self, bs):
        """Per ms rate of unbinding if tightly bound
        
        Takes:
            bs: relative Crown to Actin distance (x,y)
        Returns
            rate: per ms rate of detaching from the binding site
        """
        ## Based on the energy in the tight state
        loose_energy = self.energy(bs, "loose")
        tight_energy = self.energy(bs, "tight")
        rate = (2* # overall rate modifier
                m.tanh(0.02* # curve steepness parameter
                       abs(tight_energy-loose_energy-11))) # x offset shifter
        return float(rate)
    
    def _free_energy(self, tip_location, state):
        """Free energy of the Head
        
        Takes:
            tip_location: relative Crown to Actin distance (x,y)
            state: kinetic state of the cross-bridge, ['free'|'loose'|'tight']
        Returns:
            energy: free energy of the head in the given state
        """
        if state == "free":
            return 0
        elif state == "loose":
            return self.alphaDG + self.energy(tip_location, state)
        elif state == "tight":
            return self.etaDG + self.energy(tip_location, state)
    
    @staticmethod
    def _seg_values(tip_location):
        """Return the length and angle to the Head tip
        
        Takes:
            tip_location: relative Crown to Actin distance (x,y)
        Returns:
            (c_ang, g_len): the angle and length of the Head's springs
        """
        c_ang = m.atan2(tip_location[1], tip_location[0])
        g_len = m.hypot(tip_location[1], tip_location[0])
        return (c_ang, g_len)


class Crossbridge(Head):
    """A cross-bridge, including status of links to actin sites"""
    def __init__(self, index, parent_face, thin_face):
        """Set up the cross-bridge
    
        Parameters:
            index: the cross-bridge's index on the parent face
            parent_face: the associated thick filament face
            thin_face: the face instance opposite this cross-bridge
        """
        # Do that super() voodoo that instantiates the parent Head
        super(Crossbridge, self).__init__()
        # What is your name, where do you sit on the parent face?
        self.index = index
        # What log are you a bump upon?
        self.parent_face = parent_face
        # Remember who thou art squaring off against
        self.thin_face = thin_face
        # How can I ever find you?
        self.address = ('xb', self.parent_face.parent_filament.index,
                        self.parent_face.index, self.index)
        # Remember if thou art bound unto an actin
        self.bound_to = None # None if unbound, BindingSite object otherwise
    
    def __str__(self):
        """String representation of the cross-bridge"""
        out = '__XB_%02d__State_%s__Forces_%d_%d__'%(
            self.index, self.state,
            self.axialforce(), self.radialforce())
        return out
    
    def to_dict(self):
        """Create a JSON compatible representation of the crown
        
        Example usage: json.dumps(crown.to_dict(), indent=1)
        
        Current output includes:
            address: largest to most local, indices for finding this
            state: the free, loose, strong state of binding
            thin_face: the address of the opposing thin face
            bound_to: None or the address of the bound binding site
        """
        xbd = self.__dict__.copy()
        xbd.pop('_timestep')
        xbd.pop('index')
        xbd.pop('c')
        xbd.pop('g')
        xbd.pop('parent_face')
        if xbd['bound_to'] is not None:
            xbd['bound_to'] = xbd['bound_to'].address
        xbd['thin_face'] = xbd['thin_face'].address
        return xbd
       
    def from_dict(self, xbd):
        """ Load values from a crossbridge dict. Values read in correspond
        to the current output documented in to_dict.
        """
        # Check for index mismatch
        read, current = tuple(xbd['address']), self.address
        assert read==current, "index mismatch at %s/%s"%(read, current)
        # Local keys
        self.state = xbd['state']
        self.etaDG = xbd['etaDG']
        self.alphaDG = xbd['alphaDG']
        # Sub-structure and remote keys
        self.thin_face = self.parent_face.parent_filament.parent_lattice.\
                resolve_address(xbd['thin_face'])
        if xbd['bound_to'] is None:
            self.bound_to = None
        else:
            self.bound_to = self.parent_face.parent_filament.parent_lattice.\
                resolve_address(xbd['bound_to'])
    
    def transition(self):
        """Gather the needed information and try a transition
        
        Parameters:
            None
        Returns:
            transition: string of transition ('12', '32', etc.) or None
        """
        # When unbound, try to bind, otherwise just try a transition
        if self.bound_to is None:
            # Find the lattice spacing
            lattice_spacing = self._get_lattice_spacing()
            # Find this cross-bridge's axial location
            xb_axial_loc = self._get_axial_location()
            # Find the potential binding site
            actin_site = self.thin_face.nearest(xb_axial_loc)
            actin_axial_loc = actin_site.get_axial_location()
            actin_state = actin_site.permissiveness
            # Find the axial separation
            axial_sep = actin_axial_loc - xb_axial_loc
            # Combine the two distances
            distance_to_site = (axial_sep, lattice_spacing)
            # Allow the myosin head to take it from here
            trans = super(Crossbridge, self).transition(distance_to_site,
                                                        actin_state)
            # Process changes to bound state
            if trans == '12':
                self.bound_to = actin_site
                actin_site.bind_to(self)
            else:
                assert (trans is None), 'Bound state mismatch'
        else:
            # Get the distance to the actin site
            distance_to_site = self._dist_to_bound_actin()
            actin_state = self.bound_to.permissiveness
            # Allow the myosin head to take it from here
            trans = super(Crossbridge, self).transition(distance_to_site,
                                                        actin_state)
            # Process changes to the bound state
            if trans in set(('21', '31')):
                self.bound_to.bind_to(None)
                self.bound_to = None
            else:
                assert (trans in set(('23', '32', None))) , 'State mismatch'
        return trans
    
    def axialforce(self, base_axial_loc=None, tip_axial_loc = None):
        """Gather needed information and return the axial force
    
        Parameters:
            base_axial_location: location of the crown (optional)
            tip_axial_loc: location of an attached actin node (optional)
        Returns:
            f_x: the axial force generated by the cross-bridge
        """
        # Unbound? No force!
        if self.bound_to is None:
            return 0.0
        # Else, get the distance to the bound site and run with it
        distance = self._dist_to_bound_actin(base_axial_loc, tip_axial_loc)
        # Allow the myosin head to take it from here
        return super(Crossbridge, self).axialforce(distance)
    
    def radialforce(self):
        """Gather needed information and return the radial force
    
        Parameters:
            None
        Returns:
            f_y: the radial force generated by the cross-bridge
        """
        # Unbound? No force!
        if self.bound_to is None:
            return 0.0
        # Else, get the distance to the bound site and run with it
        distance_to_site = self._dist_to_bound_actin()
        # Allow the myosin head to take it from here
        return super(Crossbridge, self).radialforce(distance_to_site)
    
    def _get_axial_location(self):
        """Find the axial location of the thick filament attachment point
        
        Parameters:
            None
        Returns:
            axial: the axial location of the cross-bridge base
        """
        axial = self.parent_face.get_axial_location(self.index)
        return axial
    
    def _dist_to_bound_actin(self, xb_axial_loc=None, tip_axial_loc=None):
        
        """Find the (x,y) distance to the bound actin
        This is the distance format used by the myosin head.
        Parameters:
            xb_axial_loc: current axial location of the crown (optional)
            tip_axial_loc: location of an attached actin node (optional)
        Returns:
            (x,y): the axial distance between the cross-bridge base and
                   the actin site (x), and the lattice spacing (y)
        """
        # Are you really bound?
        assert (self.bound_to is not None) , "Lies, you're unbound!"
        # Find the lattice spacing
        lattice_spacing = self._get_lattice_spacing()
        # Find this cross-bridge's axial location if need be
        if xb_axial_loc is None:
            xb_axial_loc = self._get_axial_location()
        # Find the distance to the bound actin site if need be
        if tip_axial_loc is None:
            tip_axial_loc = self.bound_to.get_axial_location()
        # Combine the two distances
        return (tip_axial_loc - xb_axial_loc, lattice_spacing)
    
    def _get_lattice_spacing(self):
        """Ask our superiors for lattice spacing data"""
        return self.parent_face.get_lattice_spacing()


if __name__ == '__main__':
    print("mh.py is really meant to be called as a supporting module")

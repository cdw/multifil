#!/usr/bin/env python
# encoding: utf-8
"""
tm.py - A tropomyosin filament

Create and maintain a tropomyosin filament and the subgroups that comprise it.

Created by Dave Williams on 2017-12-03.
"""

import numpy as np

class TmSite:
    """A single tm site, located over a single actin binding site. 
    
    Individual tm sites keep track of how the tropomyosin chain is regulating
    the binding site: states, links to nearest neighbors, etc

    Kinetics
    --------

    The equilibrium state of a transition, K, is given by the ratio of forward
    to reverse transition rates, K = k_forward/k_reverse. Taking our kinetics
    from Tanner 2007 we define the three equilibrium states and their attendant
    forward transition rates. The equilibrium rates, K1, K2, and K3,
    correspond to the transitions thus:
    
        Tm+Tn+Ca <-K1-> Tm+Tn.Ca <-K2-> Tm.Tn.Ca <-K3-> Tm+Tn+Ca
    
    The forward transition rates are labeled as k_12, k_23, k_31. Reverse
    transitions are calculated from the above values. K1 and K2 are explicitly
    defined. In Tanner2007 K3 is stated to be obtained from the balancing 
    equation for equilibrium conditions, K1*K2*K3 = 1, but is instead defined
    directly. This implies non-equilibrium conditions, fair enough. 

    Only K1 is dependent on [Ca], although it would make sense for K3 to be so
    as well as Ca2+ detachment is surely dependent on [Ca]. 
   
    No rates include a temperature dependence. 
    """

    def __init__(self, parent_tm, binding_site, index):
        """ A single tropomyosin site, paired to a binding site

        Parameters
        ----------
        parent_tm: tropomyosin class instance
            the tropomyosin chain on which this site is located
        binding_site: actin binding site class instance
            the binding site this tm site regulates
        index: int
            where this tm site is along the tm chain
        """
        ## Who are you and where do you come from?
        self.parent_tm = parent_tm
        self.index = index
        self.address = ("tm_site", parent_tm.parent_thin.index, 
                        parent_tm.index, index)
        ## What are you regulating?
        self.binding_site = binding_site
        self.binding_site.tm_site = self
        self.state = 0
        ## Kinetics from Tanner 2007 and Tanner Thesis
        K1 = 1e5         # per mole Ca
        K2 = 10          # unitless
        K3 = 1e6         # unitless
        k_12 = 5e5       # per mole Ca per sec
        k_23 = 10        # per sec
        k_31 = 5         # per sec
        s_per_ms = 1e-3  # seconds per millisecond
        k_12 *= s_per_ms # convert rates from 
        k_23 *= s_per_ms # occurrences per second
        k_31 *= s_per_ms # to occurrences per ms
        self._K1, self._K2, self._K3 = K1, K2, K3
        self._k_12, self._k_23, self._k_31 = k_12, k_23, k_31

    def __str__(self):
        """Representation of the tmsite for printing"""
        out = "TmSite #%02i State:%i Loc:%04.0f"%(
            self.index, self.state, self.axial_location)
        return out

    def to_dict(self):
        """Create a JSON compatible representation of the tropomyosin site

        Usage example:json.dumps(tmsite.to_dict(), indent=1)

        Current output includes:
            address: largest to most local, indices for finding this
            binding_site: address of binding site
            state: kinetic state of site
            binding_influence: what the state tells you
            span: how far an activation spreads
            _k_12 - _k_31: transition rates
            _K1 - _K3: kinetic balances for reverse rates

        Returns
        -------
        tmsd: dict
            tropomyosin site dictionary
        """
        tmsd = self.__dict__.copy()
        tmsd.pop('parent_tm')
        tmsd.pop('index')
        tmsd['binding_site'] = tmsd['binding_site'].address
        return tmsd

    def from_dict(self, tmsd):
        """ Load values from a tropomyosin site dict. 
        
        Values read correspond to output documented in :to_dict:.
        """
        # Check for index mismatch
        read, current = tuple(tmsd['address']), self.address
        assert read==current, "index mismatch at %s/%s"%(read, current)
        # Local/remote keys
        self.binding_site = self.parent_tm.parent_thin.parent_lattice.\
                resolve_address(tmsd['binding_site']) 
        self.binding_site.tm_site = self
        # Local keys
        self.state = tmsd['_state']
        self._k_12 = tmsd['_k_12'] 
        self._k_23 = tmsd['_k_23']
        self._k_31 = tmsd['_k_31']
        self._K1 = tmsd['_K1']
        self._K2 = tmsd['_K2']
        self._K3 = tmsd['_K3']
        return

    @property
    def timestep(self):
        """Timestep size is stored at the half-sarcomere level"""
        return self.parent_tm.parent_thin.parent_lattice.timestep_len

    @property
    def pCa(self):
        """pCa stored at the half-sarcomere level"""
        pCa = self.parent_tm.parent_thin.parent_lattice.pCa 
        assert pCa>0, "pCa must be given in positive units by convention"
        return pCa

    @property
    def ca(self):
        """The calcium concentration stored at the half-sarcomere level"""
        Ca = 10.0**(-self.pCa)
        return Ca

    @property
    def axial_location(self):
        """Axial location of the bs we are piggybacking on"""
        return self.binding_site.axial_location

    @property
    def span(self):
        """What is the span of cooperative activation for this tm site?
        
        The span (state 2 coercion of adjacent sites to state 1 from 
        state 0) is based on the current tension at the binding site 
        co-located under this tropomyosin site. 

        Notes
        -----
        The functional form of the span is determined by:
            $$span = 0.5 * base (1 + tanh(steep*(force50 + f)))$$
        Where $span$ is the actual span, $base$ is the resting (no force) 
        span distance, $steep$ is how steep the decrease in span is, 
        $force50$ is the force at which the span has decreased by half, and 
        f is the current effective axial force of the thin filament, an 
        estimate of the tension along the thin filament. 

        These properties are stored at the tropomyosin chain level as they 
        are material properties of the entire chain.
        """
        base = self.parent_tm.span_base
        steep = self.parent_tm.span_steep
        f50 = self.parent_tm.span_force50
        f = self.binding_site.tension
        span = 0.5 * base * (1 - np.tanh(steep * (f50 + f)))
        return span
        
    @property
    def state(self):
        """Get the current state of the tm site"""
        return self._state

    @state.setter
    def state(self, new_state):
        """Set the state and thus the binding influence"""
        self._state = new_state
        self.binding_influence = {0:0, 1:0, 2:1}[new_state]

    def _r12(self):
        """Rate of Ca and TnC association, conditional on [Ca]"""
        forward = self._k_12 * self.ca
        return forward

    def _r21(self):
        """Rate of Ca detachment from TnC, conditional on [Ca]"""
        forward = self._r12()
        reverse = forward / self._K1
        return reverse

    def _r23(self):
        """Rate of TnI TnC association"""
        forward = self._k_23
        return forward

    def _r32(self):
        """Rate of TnI TnC detachment"""
        forward = self._r23()
        reverse = forward / self._K2
        return reverse

    def _r31(self):
        """Rate of Ca disassociation induced TnI TnC disassociation
        TODO: figure out calcium dependence
        """
        forward = self._k_31 
        return forward

    def _r13(self):
        """Rate of simultaneous Ca binding and TnI TnC association.
        Should be quite low. 
        TODO: figure out calcium dependence
        """
        forward = self._r31()
        reverse = forward / self._K3
        return reverse

    def _prob(self, rate):
        """ Convert a rate to a probability, based on the current timestep
        length and the assumption that the rate is for a Poisson process.
        We are asking, what is the probability that at least one Poisson
        distributed value would occur during the timestep.

        Parameters
        ----------
            rate: float
                a per millisecond rate to convert to probability
        Returns
        -------
            probability: float
                the probability the event occurs during a timestep
                of length determined by self.timestep
        """
        return 1 - np.exp(-rate*self.timestep)

    @staticmethod
    def _forward_backward(forward_p, backward_p, rand):
        """Transition forward or backward based on random variable, return
        "forward", "backward", or "none"
        """
        if rand<forward_p:
            return "forward"
        elif rand>(1-backward_p):
            return "backward"
        return "none"

    def transition(self):
        """Transition from one state to the next, or back, or don't """
        rand = np.random.random()
        if self.state==0:
            f, b = self._prob(self._r12()), self._prob(self._r13())
            trans_word = self._forward_backward(f, b, rand)
            self.state = {"forward":1, "backward":2, "none":0}[trans_word]
        elif self.state==1:
            f, b = self._prob(self._r23()), self._prob(self._r21())
            trans_word = self._forward_backward(f, b, rand)
            self.state = {"forward":2, "backward":0, "none":1}[trans_word]
        elif self.state==2:
            if self.binding_site.state == 1:
                return # can't transition if bound
            f, b = self._prob(self._r31()), self._prob(self._r32())
            trans_word = self._forward_backward(f, b, rand)
            self.state = {"forward":0, "backward":1, "none":2}[trans_word]
        assert self.state in (0,1,2), "Tropomyosin state has invalid value"
        return


class Tropomyosin:
    """Regulate the binding permissiveness of actin strands. 
    
    Tropomyosin stands in for both the Tm and the Tn. 
    
    Structure
    ---------

    The arrangement of the tropomyosin on the thin filament is represented as
    having an ability to be calcium regulated at each binding site with the
    option to spread that calcium regulation to adjacent binding sites. I am
    only grudgingly accepting of this as a representation of the TmTn
    interaction structure, but it is a reasonable first pass. 
    """

    def __init__(self, parent_thin, binding_sites, index):
        """A strand of tropomyosin chains
        
        Save the binding sites along a set of tropomyosin strands, 
        in preparation for altering availability of binding sites. 
        
        Parameters
        ----------
        parent_thin: thin filament instance
            thin filament on which the tropomyosin lives
        binding_sites: list of list of binding site instances
            binding sites on this tm string
        index: int
            which tm chain this is on the thin filament
        """
        ## Who are you and where do you come from?
        self.parent_thin = parent_thin
        self.index = index
        self.address = ("tropomyosin", parent_thin.index, index)
        ## What is your population?
        self.sites = [TmSite(self, bs, ind) for ind, bs in
                      enumerate(binding_sites)]
        ## How does activation spread? 
        # Material properties belong to tm chain, but actual span is 
        # calculated at the site level (where tension is experienced)
        self.span_base = 37  # influence span (Tanner 2007)
        self.span_base *= 1.2  # influence span increase 
        self.span_steep = 1  # how steep the influence curve is
        self.span_force50 = -20 # force at which span is decreased by half

    def to_dict(self):
        """Create a JSON compatible representation of the tropomyosin chain

        Usage example:json.dumps(tm.to_dict(), indent=1)

        Current output includes:
            address: largest to most local, indices for finding this
            sites: tm sites
        """
        tmd = self.__dict__.copy()
        tmd.pop('parent_thin')
        tmd.pop('index')
        tmd['sites'] = [site.to_dict() for site in tmd['sites']]
        return tmd

    def from_dict(self, tmd):
        """ Load values from a tropomyosin dict. Values read correspond to 
        the current output documented in to_dict.
        """
        # Check for index mismatch
        read, current = tuple(tmd['address']), self.address
        assert read==current, "index mismatch at %s/%s"%(read, current)
        # Sub-structures
        for data, site in zip(tmd['sites'], self.sites):
            site.from_dict(data)
        self.span = tmd['span']

    def resolve_address(self, address):
        """Give back a link to the object specified in the address
        We should only see addresses starting with 'tm_site'
        """
        if address[0] == 'tm_site':
            return self.sites[address[3]]
        import warnings
        warnings.warn("Unresolvable address: %s"%str(address))

    @property
    def states(self):
        """States of the contained TmSites (for monitoring)"""
        return [site.state for site in self.sites]

    def transition(self):
        """Chunk through all binding sites, transition if need be"""
        for site in self.sites:
            site.transition()
        # Spread activation
        self._spread_activation()

    def _spread_activation(self):
        """"Spread activation along the filament"""
        # If nothing is in state 2, we're done
        if max([site.state for site in self.sites])<2:
            return
        # Find all my axial locations
        locs = np.array([site.axial_location for site in self.sites])
        # Chunk through each site
        for site in self.sites:
            if site.state == 2:
                loc = site.axial_location
                span = site.span
                near_inds = np.nonzero(np.abs(locs-loc)<span)[0]
                near = [self.sites[index] for index in near_inds]
                for site in near:
                    if site.state != 2:
                        site.state = 1
        return



if __name__ == '__main__':
    print("tm.py is really meant to be called as a supporting module")

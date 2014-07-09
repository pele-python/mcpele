from __future__ import division
import numpy as np
from pele.utils.rotations import vector_random_uniform_hypersphere
from pele.potentials import Harmonic
from mcpele.monte_carlo import Metropolis_MCrunner
import unittest
import logging


class TestMetropolis(unittest.TestCase):
    
    def setUp(self):
        self.natoms = 4
        self.k=1
        self.Emax = 2 #this choice is fundamentally arbitrary, it's only used to generate the initial configuration
        #self.temperatures = [0.2,0.27,0.362,0.487,0.65,0.88,1.18,1.6]
        self.temperatures = [0.2,0.362,0.65,1.18]
        self.stepsize=1
        self.niter=2e6
        self.adjustf = 0.9
        self.adjust_niter = 0.3*self.niter
        self.radius = 1e10
            
    def test_heat_capacity_3D_com(self):
        self.bdim = 3
        self.ndim = self.natoms*self.bdim
        self.origin = np.zeros(self.ndim)
        potential = Harmonic(self.origin, self.k, bdim=self.bdim, com=True)
                        
        print "3D com"
        for T in self.temperatures:
                        
            mcrunner = Metropolis_MCrunner(potential, self.origin, T, self.stepsize, self.niter, hEmax = 100, 
                                           adjustf = self.adjustf, adjustf_niter = self.adjust_niter, radius=self.radius, 
                                           bdim=self.bdim)
            #MCMC 
            mcrunner.run()
            #collect the results
            binenergy, hist, mean, variance = mcrunner.get_histogram()
            
            cv =  variance/(T*T)
            cv_true = (self.natoms-1)*self.bdim/2.0
            
            print cv, cv_true
            
            self.assertLess(abs(cv-cv_true),1e-1,'failed for temperature {}, cv = {}'.format(T,cv))
    
    def test_heat_capacity_3D(self):
        self.bdim = 3
        self.ndim = self.natoms*self.bdim
        self.origin = np.zeros(self.ndim)
        potential = Harmonic(self.origin, self.k, bdim=self.bdim, com=False)
        
        print "3D"
        for T in self.temperatures:
                        
            mcrunner = Metropolis_MCrunner(potential, self.origin, T, self.stepsize, self.niter, hEmax = 100, 
                                           adjustf = self.adjustf, adjustf_niter = self.adjust_niter, radius=self.radius, 
                                           bdim=self.bdim)
            #MCMC 
            mcrunner.run()
            #collect the results
            binenergy, hist, mean, variance = mcrunner.get_histogram()
            
            cv =  variance/(T*T)
            cv_true = self.natoms*self.bdim/2.0
            
            print cv, cv_true
            
            self.assertLess(abs(cv-cv_true),1e-1,'failed for temperature {}, cv = {}'.format(T,cv))
    
    def test_heat_capacity_2D_com(self):
        self.bdim = 2
        self.ndim = self.natoms*self.bdim
        self.origin = np.zeros(self.ndim)
        potential = Harmonic(self.origin,self.k, bdim=self.bdim,com=True)
        #self.start_coords = vector_random_uniform_hypersphere(self.ndim) * np.sqrt(2*self.Emax) #coordinates sampled from Pow(ndim)
                        
        print "2D com"
        for T in self.temperatures:
                        
            mcrunner = Metropolis_MCrunner(potential, self.origin, T, self.stepsize, self.niter, hEmax = 100, 
                                           adjustf = self.adjustf, adjustf_niter = self.adjust_niter, radius=self.radius, 
                                           bdim=self.bdim)
            #MCMC 
            mcrunner.run()
            #collect the results
            binenergy, hist, mean, variance = mcrunner.get_histogram()
            
            cv =  variance/(T*T)
            cv_true = (self.natoms-1)*self.bdim/2.0
            
            print cv, cv_true
            
            self.assertLess(abs(cv-cv_true),1e-1,'failed for temperature {}, cv = {}'.format(T,cv))
    
    def test_heat_capacity_2D(self):
        self.bdim = 2
        self.ndim = self.natoms*self.bdim
        self.origin = np.zeros(self.ndim)
        self.potential = Harmonic(self.origin, self.k, bdim=self.bdim, com=False)
                        
        print "2D"
        for T in self.temperatures:
                        
            mcrunner = Metropolis_MCrunner(self.potential, self.origin, T, self.stepsize, self.niter, hEmax = 100, 
                                           adjustf = self.adjustf, adjustf_niter = self.adjust_niter, radius=self.radius, 
                                           bdim=self.bdim)
            #MCMC 
            mcrunner.run()
            #collect the results
            binenergy, hist, mean, variance = mcrunner.get_histogram()
            
            cv =  variance/(T*T)
            cv_true = self.natoms*self.bdim/2.0
            
            print cv, cv_true
            
            self.assertLess(abs(cv-cv_true),1e-1,'failed for temperature {}, cv = {}'.format(T,cv))

if __name__ == "__main__":
    logging.basicConfig(filename='Metropolis_mcrunner.log',level=logging.DEBUG)
    unittest.main()

from __future__ import division
import numpy as np
from pele.potentials import Harmonic
from mcpele.monte_carlo import Metropolis_MCrunner
import unittest
import logging


class TestMetropolisGlobal(unittest.TestCase):
    
    def setUp(self):
        self.natoms = 4
        self.k=1
        self.Emax = 2 #this choice is fundamentally arbitrary, it's only used to generate the initial configuration
        #self.temperatures = [0.2,0.27,0.362,0.487,0.65,0.88,1.18,1.6]
        self.temperatures = [0.2, 0.362, 0.65, 1.18]
        self.stepsize = 1
        self.niter = 2e6
        self.adjustf = 0.9
        self.adjust_niter = 0.3 * self.niter
        self.radius = 1e10
        self.tolerance = 2e-1
        self.bdim = 3
        self.ndim = self.natoms * self.bdim
        self.origin = np.zeros(self.ndim)
        self.seeds = dict(takestep=42, metropolis=44)
        self.potential = Harmonic(self.origin, self.k, bdim=self.bdim, com=False)
        self.mcrunner = Metropolis_MCrunner(self.potential, self.origin, 1, self.stepsize, self.niter, hEmax = 100, 
                                           adjustf = self.adjustf, adjustf_niter = self.adjust_niter, radius=self.radius, 
                                           bdim=self.bdim, single=False, seeds=self.seeds)
    
    def test_set_control(self):
        self.mcrunner.set_control(10)
        self.mcrunner.run()
        self.assertEqual(self.mcrunner.get_temperature(), self.mcrunner.temperature)
    
    def test_dump_histogram(self):
        self.mcrunner.run()
        mean, variance = self.mcrunner.dump_histogram("test_histogram.dat")
        data = np.genfromtxt("test_histogram.dat")
        #print data
        binenergy, hist, mean2, variance2 = self.mcrunner.get_histogram()
        self.assertListEqual(np.ndarray.tolist(binenergy), np.ndarray.tolist(data[:,0]))
        self.assertListEqual(np.ndarray.tolist(hist), np.ndarray.tolist(data[:,1]))
        self.assertEqual(mean, mean2)
        self.assertEqual(variance, variance2)
        count = self.mcrunner.histogram.get_count()
        self.assertEqual(count, 1400000)
    
    def test_heat_capacity_3D_com(self):
        self.bdim = 3
        self.ndim = self.natoms*self.bdim
        self.origin = np.zeros(self.ndim)
        potential = Harmonic(self.origin, self.k, bdim=self.bdim, com=True)
        mcrunner = Metropolis_MCrunner(potential, self.origin, 1, self.stepsize, self.niter, hEmax = 100, 
                                           adjustf = self.adjustf, adjustf_niter = self.adjust_niter, radius=self.radius, 
                                           bdim=self.bdim, seeds=self.seeds)
        mcrunner.run()
                        
        print "3D com"
        for T in self.temperatures:
                        
            mcrunner = Metropolis_MCrunner(potential, self.origin, T, self.stepsize, self.niter, hEmax = 100, 
                                           adjustf = self.adjustf, adjustf_niter = self.adjust_niter, radius=self.radius, 
                                           bdim=self.bdim, seeds=self.seeds)
            #MCMC 
            mcrunner.run()
            #collect the results
            binenergy, hist, mean, variance = mcrunner.get_histogram()
            
            cv =  variance/(T*T)
            cv_true = (self.natoms-1)*self.bdim/2.0
            
            print cv, cv_true
            
            self.assertLess(abs(cv-cv_true),self.tolerance,'failed for temperature {}, cv = {}'.format(T,cv))
    
    def test_heat_capacity_3D(self):
        self.bdim = 3
        self.ndim = self.natoms*self.bdim
        self.origin = np.zeros(self.ndim)
        potential = Harmonic(self.origin, self.k, bdim=self.bdim, com=False)
        
        print "3D"
        for T in self.temperatures:
                        
            mcrunner = Metropolis_MCrunner(potential, self.origin, T, self.stepsize, self.niter, hEmax = 100, 
                                           adjustf = self.adjustf, adjustf_niter = self.adjust_niter, radius=self.radius, 
                                           bdim=self.bdim, seeds=self.seeds)
            #MCMC 
            mcrunner.run()
            #collect the results
            binenergy, hist, mean, variance = mcrunner.get_histogram()
            
            cv =  variance/(T*T)
            cv_true = self.natoms*self.bdim/2.0
            
            print cv, cv_true
            
            self.assertLess(abs(cv-cv_true),self.tolerance,'failed for temperature {}, cv = {}'.format(T,cv))
    
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
                                           bdim=self.bdim, seeds=self.seeds)
            #MCMC 
            mcrunner.run()
            #collect the results
            binenergy, hist, mean, variance = mcrunner.get_histogram()
            
            cv =  variance/(T*T)
            cv_true = (self.natoms-1)*self.bdim/2.0
            
            print cv, cv_true
            
            self.assertLess(abs(cv-cv_true),self.tolerance,'failed for temperature {}, cv = {}'.format(T,cv))
    
    def test_heat_capacity_2D(self):
        self.bdim = 2
        self.ndim = self.natoms*self.bdim
        self.origin = np.zeros(self.ndim)
        self.potential = Harmonic(self.origin, self.k, bdim=self.bdim, com=False)
                        
        print "2D"
        for T in self.temperatures:
                        
            mcrunner = Metropolis_MCrunner(self.potential, self.origin, T, self.stepsize, self.niter, hEmax = 100, 
                                           adjustf = self.adjustf, adjustf_niter = self.adjust_niter, radius=self.radius, 
                                           bdim=self.bdim, seeds=self.seeds)
            #MCMC 
            mcrunner.run()
            #collect the results
            binenergy, hist, mean, variance = mcrunner.get_histogram()
            
            cv =  variance/(T*T)
            cv_true = self.natoms*self.bdim/2.0
            
            print cv, cv_true
            
            self.assertLess(abs(cv-cv_true),self.tolerance,'failed for temperature {}, cv = {}'.format(T,cv))
            
    def test_composite_moves(self):
        self.bdim = 2
        self.ndim = self.natoms*self.bdim
        self.origin = np.zeros(self.ndim)
        self.potential = Harmonic(self.origin, self.k, bdim=self.bdim, com=False)

class TestMetropolisSingle(unittest.TestCase):
    
    def setUp(self):
        self.natoms = 4
        self.k=1
        self.Emax = 2 #this choice is fundamentally arbitrary, it's only used to generate the initial configuration
        #self.temperatures = [0.2,0.27,0.362,0.487,0.65,0.88,1.18,1.6]
        self.temperatures = [0.2, 0.362, 0.65, 1.18]
        self.stepsize = 1
        self.niter = 2e6
        self.adjustf = 0.9
        self.adjust_niter = 0.3 * self.niter
        self.radius = 1e10
        self.tolerance = 2e-1
        self.bdim = 3
        self.ndim = self.natoms * self.bdim
        self.origin = np.zeros(self.ndim)
        self.potential = Harmonic(self.origin, self.k, bdim=self.bdim, com=False)
        self.mcrunner = Metropolis_MCrunner(self.potential, self.origin, 1, self.stepsize, self.niter, hEmax = 100, 
                                           adjustf = self.adjustf, adjustf_niter = self.adjust_niter, radius=self.radius, 
                                           bdim=self.bdim, single=True)
    
    def test_set_control(self):
        self.mcrunner.set_control(10)
        self.mcrunner.run()
        self.assertEqual(self.mcrunner.get_temperature(), self.mcrunner.temperature)
    
    def test_dump_histogram(self):
        self.mcrunner.run()
        mean, variance = self.mcrunner.dump_histogram("test_histogram.dat")
        data = np.genfromtxt("test_histogram.dat")
        #print data
        binenergy, hist, mean2, variance2 = self.mcrunner.get_histogram()
        self.assertListEqual(np.ndarray.tolist(binenergy), np.ndarray.tolist(data[:,0]))
        self.assertListEqual(np.ndarray.tolist(hist), np.ndarray.tolist(data[:,1]))
        self.assertEqual(mean, mean2)
        self.assertEqual(variance, variance2)
    
    def test_heat_capacity_3D_com(self):
        self.bdim = 3
        self.ndim = self.natoms*self.bdim
        self.origin = np.zeros(self.ndim)
        potential = Harmonic(self.origin, self.k, bdim=self.bdim, com=True)
        mcrunner = Metropolis_MCrunner(potential, self.origin, 1, self.stepsize, self.niter, hEmax = 100, 
                                           adjustf = self.adjustf, adjustf_niter = self.adjust_niter, radius=self.radius, 
                                           bdim=self.bdim)
        mcrunner.run()
                        
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
            
            self.assertLess(abs(cv-cv_true),self.tolerance,'failed for temperature {}, cv = {}'.format(T,cv))
    
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
            
            self.assertLess(abs(cv-cv_true),self.tolerance,'failed for temperature {}, cv = {}'.format(T,cv))
    
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
            
            self.assertLess(abs(cv-cv_true),self.tolerance,'failed for temperature {}, cv = {}'.format(T,cv))
    
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
            
            self.assertLess(abs(cv-cv_true),self.tolerance,'failed for temperature {}, cv = {}'.format(T,cv))
            
    def test_composite_moves(self):
        self.bdim = 2
        self.ndim = self.natoms*self.bdim
        self.origin = np.zeros(self.ndim)
        self.potential = Harmonic(self.origin, self.k, bdim=self.bdim, com=False)

if __name__ == "__main__":
    logging.basicConfig(filename='Metropolis_mcrunner.log',level=logging.DEBUG)
    unittest.main()

"""
Useful functions to perform quick calculations.

"""
import numpy as np 

# Default parameters
c = 343                    # Speed ​​of sound in air.

def cuttoff_frec(a):
    """ Cutoff frequency calculation for validity of the plane wave hypothesis.
    Verify that you are working below it (above this there are modes in other 
    dimensions).
    This frequency is defined by the radius "a" of the stage."""
    return (1.84*c)/(2*np.pi*a)

def expansion_chamber(f,n=0):
    """ It is used to calculate the length that an expansion chamber must have,
    for a certain frequency value. """
    return ((2*n + 1)*c)/(4*f)
    
def side_branch(f,n=0):
    """ It is used to calculate the length that an side branch must have,
    for a certain frequency value. """
    return ((2*n + 1)*c)/(4*f)
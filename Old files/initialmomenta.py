#Initial momenta module
import numpy as np
import random as rn

def init_mom(N, temp):
  momenta = np.random.normal(0.0, np.sqrt(temp), (N,3))
  avmom = np.zeros((1,3),dtype=float)
  avmom = np.sum(momenta, axis=0)/N
  momenta = momenta - avmom
  return momenta


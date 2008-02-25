import numpy as npy

from ase.optimize import Dynamics
from ase.data import atomic_masses


class MolecularDynamics(Dynamics):
    def __init__(self, atoms):
        Dynamics.__init__(self, atoms, logfile=None)

    def run(self, dt, steps=50):
        """Integrate equation of motion.

        Parameters
        ----------
        dt: float
            Time step.
        steps: int
            Number of steps (defaults to 50).
        """
        
        m = self.atoms.get_masses()
        if m is None:
            m = atomic_masses[self.atoms.get_atomic_numbers()]
        m = m.reshape((-1, 1))

        f = self.atoms.get_forces()

        if self.atoms.get_momenta() is None:
            self.atoms.set_momenta(npy.zeros_like(f))

        for step in xrange(steps):
            f = self.step(f, dt, m)
            self.call_observers(step)

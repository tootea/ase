# -*- coding: utf-8 -*-

"""Vibrational modes."""

import pickle
from math import sin, pi, sqrt
from os import remove
from os.path import isfile

import numpy as npy

import ase.units as units
from ase.data import atomic_masses
from ase.io.trajectory import PickleTrajectory
from ase.parallel import rank, barrier


class Vibrations:
    """Class for calculating vibrational modes using finite difference."""
    def __init__(self, atoms, indices=None, name='vib', delta=0.01):
        """Create Vibrations object.

        Parameters
        ==========
        atoms: Atoms object
            The atoms to work on.
        indices: list of int
            List of indices of atoms to vibrate.  Default behavior is
            to vibrate all atoms.
        name: str
            Name to use for files.
        delta: float
            Magnitude of displacements.

        Example
        =======

        >>> from ase import *
        >>> from ase.vibrations import Vibrations
        >>> n2 = Atoms('N2', [(0, 0, 0), (0, 0, 1.1)],
        ...            calculator=EMT())
        >>> QuasiNewton(n2).run(fmax=0.01)
        QuasiNewton:   0        0.042171       2.9357
        QuasiNewton:   1        0.016313       1.6546
        QuasiNewton:   2        0.000131       0.1534
        QuasiNewton:   3        0.000000       0.0093
        >>> vib = Vibrations(n2)
        >>> vib.run()
        >>> vib.summary()
        ---------------------
          #    meV     cm^-1
        ---------------------
          0    1.7i     13.5i
          1    1.7i     13.5i
          2    0.0i      0.0i
          3    0.0       0.0 
          4    0.0       0.0 
          5  232.8    1877.9 
        ---------------------
        Zero-point energy: 0.116 eV
        
        >>> vib.write_mode(-1)  # write last mode to trajectory file

        """
        
	self.atoms = atoms
        if indices is None:
            indices = range(len(atoms))
        self.indices = npy.asarray(indices)
        self.name = name
        self.delta = delta
        self.H = None

    def run(self):
        """Run the vibration calculations.

        This will calculate the forces for 6 dislpacements per atom
        ±x, ±y, ±z.  Only those calculations that are not already done
        will be started."""
        
        p = self.atoms.positions.copy()
        for a in self.indices:
            for i in range(3):
                for sign in [-1, 1]:
                    filename = '%s.%d%s%s.pckl' % (self.name, a,
                                                   'xyz'[i], ' +-'[sign])
                    if isfile(filename):
                        continue
                    barrier()
                    if rank == 0:
                        fd = open(filename, 'w')
                    self.atoms.positions[a, i] = p[a, i] + sign * self.delta
                    forces = self.atoms.get_forces()
                    if rank == 0:
                        pickle.dump(forces, fd)
                        fd.close()
                    self.atoms.positions[a, i] = p[a, i]
        self.atoms.set_positions(p)

    def clean(self):
        for a in self.indices:
            for i in 'xyz':
                for sign in '-+':
                    name = '%s.%d%s%s.pckl' % (self.name, a, i, sign)
                    if isfile(name):
                        remove(name)
        
    def read(self):
        n = 3 * len(self.indices)
        H = npy.empty((n, n))
        r = 0
        for a in self.indices:
            for i in 'xyz':
                name = '%s.%d%s' % (self.name, a, i)
                fminus = pickle.load(open(name + '-.pckl'))[self.indices]
                fplus = pickle.load(open(name + '+.pckl'))[self.indices]
                H[r] = (fminus - fplus).ravel() / (4 * self.delta)
                r += 1
        H += H.copy().T
        self.H = H
        try:
            m = self.atoms.get_masses()
        except KeyError:
            m = atomic_masses[self.atoms.get_atomic_numbers()]
        self.im = npy.repeat(m[self.indices]**-0.5, 3)
        Q = npy.diag(self.im)
        omega2, modes = npy.linalg.eigh(npy.dot(Q, npy.dot(H, Q)))
        self.modes = modes.T.copy()

        # Conversion factor:
        s = units._hbar * 1e10 / sqrt(units._e * units._amu)
        self.hnu = s * omega2.astype(complex)**0.5

    def get_energies(self):
        """Get vibration energies in eV."""
        if self.H is None:
            self.read()
        return self.hnu

    def get_frequencies(self):
        """Get vibration frequencies in cm^-1."""
        s = 0.01 * units._e / units._c / units._hplanck
        return s * self.get_energies()

    def summary(self):
        hnu = self.get_energies()
        s = 0.01 * units._e / units._c / units._hplanck
        print '---------------------'
        print '  #    meV     cm^-1'
        print '---------------------'
        for n, e in enumerate(hnu):
            if e.imag != 0:
                c = 'i'
                e = e.imag
            else:
                c = ' '
            print '%3d %6.1f%s  %7.1f%s' % (n, 1000 * e, c, s * e, c)
        print '---------------------'
        print 'Zero-point energy: %.3f eV' % self.get_zero_point_energy()
        print

    def get_zero_point_energy(self):
        return 0.5 * self.hnu.real.sum()

    def get_mode(self, n):
        mode = npy.zeros((len(self.atoms), 3))
        mode[self.indices] = (self.modes[n] * self.im).reshape((-1, 3))
        return mode

    def write_mode(self, n, kT=units.kB * 300, nimages=30):
        """Write mode to trajectory file."""
        mode = self.get_mode(n) * sqrt(kT / self.hnu[n])
        p = self.atoms.positions.copy()
        n %= 3 * len(self.indices)
        traj = PickleTrajectory('%s.%d.traj' % (self.name, n), 'w')
        calc = self.atoms.get_calculator()
        self.atoms.set_calculator()
        for x in npy.linspace(0, 2 * pi, nimages, endpoint=False):
            self.atoms.set_positions(p + sin(x) * mode)
            traj.write(self.atoms)
        self.atoms.set_positions(p)
        self.atoms.set_calculator(calc)
        traj.close()

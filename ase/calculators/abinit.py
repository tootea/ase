b"""This module defines an ASE interface to ABINIT.

http://www.abinit.org/
"""

import os
from glob import glob
from os.path import join, isfile, islink

import numpy as np

from ase.data import chemical_symbols
from ase.data import atomic_numbers
from ase.units import Bohr, Hartree
from ase.calculators.calculator import FileIOCalculator, Parameters
from ase.io.abinit import read_abinit


class Abinit(FileIOCalculator):
    """Class for doing ABINIT calculations.

    The default parameters are very close to those that the ABINIT
    Fortran code would use.  These are the exceptions::

      calc = Abinit(label='abinit', xc='LDA', diemix=0.1)
    """
    def __init__(self, label='abinit', atoms=None, 
                 xc='LDA',
                 width=0.04 * Hartree, smearing='fermi-dirac',
                 ecut=None,
                 charge=0.0,
                 npulayit=7, diemix=0.1, diemac=1.e6,
                 pps='fhi',
                 toldfe=1.0e-6,
                 scratch=None,
                 **kwargs):
        """Construct ABINIT-calculator object.

        Parameters
        ==========
        label: str
            Prefix to use for filenames (label.in, label.txt, ...).
            Default is 'abinit'.
        xc: str
            Exchange-correlation functional.  Must be one of LDA, PBE,
            revPBE, RPBE.
        kpts: list of three int
            Monkhost-Pack sampling.
        nbands: int
            Number of bands.
            For the values of occopt not equal to 0 or 2, nbands can be omitted.
        nstep: int
            Number of self-consistent field STEPS.
        width: float
            Fermi-distribution width in eV.
            Default is 0.04 Hartree.
        ecut: float
            Planewave cutoff energy in eV.
            No default.
        charge: float
            Total charge of the system.
            Default is 0.
        npulayit: int
            Number of old densities to use for Pulay mixing.
        diemix: float
            Mixing parameter between zero and one for density mixing.
        diemac: float
            Model DIElectric MACroscopic constant.
            The value of diemac should usually be bigger than 1.0d0,
            on physical grounds. If you let diemac to its default value,
            you might even never obtain the self-consistent convergence!

        Examples
        ========
        Use default values:

        >>> h = Atoms('H', calculator=Abinit())
        >>> h.center(vacuum=3.0)
        >>> e = h.get_potential_energy()

        """

        if ecut is None:
            raise ValueError('Planewave cutoff energy in eV (ecut) not set')

        if pps not in ['fhi', 'hgh', 'hgh.sc', 'hgh.k', 'tm', 'paw']:
            raise ValueError('Unexpected PP identifier %s' % pps)

        self.scratch = scratch
        
        self.species = None
        self.ppp_list = None

        self.n_entries_int = 20  # integer entries per line
        self.n_entries_float = 8  # float entries per line

        FileIOCalculator.__init__(self, label, atoms,
                                  xc=xc,
                                  width=width, smearing=smearing,
                                  ecut=ecut,
                                  charge=charge,
                                  npulayit=npulayit,
                                  diemix=diemix, diemac=diemac,
                                  pps=pps,
                                  toldfe=toldfe,
                                  **kwargs)

    def write_input(self, atoms, properties, system_changes):
        """Write input parameters to files-file."""

        FileIOCalculator.write_input(self, atoms, properties, system_changes)

        if 'numbers' in system_changes or 'magmoms' in system_changes:
            self.initialize(atoms)

        dir, prefix = self.split_label()

        fh = open(self.label + '.files', 'w')

        fh.write('%s\n' % (self.label+'.in')) # input
        fh.write('%s\n' % (self.label+'.txt')) # output
        fh.write('%s\n' % (self.label+'i')) # input
        fh.write('%s\n' % (self.label+'o')) # output
        
        # scratch files
        scratch = self.scratch
        if scratch is None:
            scratch = dir
        if not os.path.exists(scratch):
            os.makedirs(scratch)
        fh.write('%s\n' % (os.path.join(scratch, prefix + '.abinit')))
        # Provide the psp files
        for ppp in self.ppp_list:
            fh.write('%s\n' % (ppp)) # psp file path

        fh.close()

        param = self.parameters
        param.write(self.label + '.param')

        fh = open(self.label + '.in', 'w')
        inp = {}
        inp.update(param)
        del inp['width']
        del inp['smearing']
        del inp['xc']
        del inp['kpts']
        del inp['toldfe']
        del inp['pps']

        inp['tsmear'] = param.width
        inp['occopt'] = 3
        assert param.smearing == 'fermi-dirac'
        inp['natom'] = len(atoms)

        if 'nbands' in param:
            inp['nband'] = param.nbands
            del inp['nbands']

        inp['ixc'] = {'LDA': 7,
                      'PBE': 11,
                      'revPBE': 14,
                      'RPBE': 15,
                      'WC': 23}[param.xc]

        magmoms = atoms.get_initial_magnetic_moments()
        if magmoms.any():
            inp['nsppol'] = 2
            fh.write('spinat\n')
            for n, M in enumerate(magmoms):
                fh.write('%.14f %.14f %.14f\n' % (0, 0, M))
        else:
            inp['nsppol'] = 1

        for key in sorted(inp.keys()):
            value = inp[key]
            print key,value
            if isinstance(value, list):
                fh.write('%s %s\n' % ('%block', key))  # don't confuse '%b'
                for line in value:
                    fh.write(' '.join(['%s' % x for x in line]) + '\n')
                fh.write('%endblock %s\n' % key)

            unit = keys_with_units.get(key)
            if unit is None:
                fh.write('%s %s\n' % (key, value))
            else:
                if 'fs**2' in unit:
                    value /= fs**2
                elif 'fs' in unit:
                    value /= fs
                fh.write('%s %f %s\n' % (key, value, unit))

        fh.write('#Definition of the unit cell\n')
        fh.write('acell\n')
        fh.write('%.14f %.14f %.14f Angstrom\n' %  (1.0, 1.0, 1.0))
        fh.write('rprim\n')
        for v in atoms.cell:
            fh.write('%.14f %.14f %.14f\n' %  tuple(v))

        fh.write('chkprim 0 # Allow non-primitive cells\n')

        fh.write('#Definition of the atom types\n')
        fh.write('ntypat %d\n' % (len(self.species)))
        fh.write('znucl')
        for n, Z in enumerate(self.species):
            fh.write(' %d' % (Z))
        fh.write('\n')
        fh.write('#Enumerate different atomic species\n')
        fh.write('typat')
        fh.write('\n')
        self.types = []
        for Z in atoms.numbers:
            for n, Zs in enumerate(self.species):
                if Z == Zs:
                    self.types.append(n+1)
        for n, type in enumerate(self.types):
            fh.write(' %d' % (type))
            if n > 1 and ((n % self.n_entries_int) == 1):
                fh.write('\n')
        fh.write('\n')

        fh.write('#Definition of the atoms\n')
        fh.write('xangst\n')
        a = 0
        for pos, Z in zip(atoms.positions, atoms.numbers):
            a += 1
            fh.write('%.14f %.14f %.14f\n' %  tuple(pos))

        if param.kpts is not None:
            fh.write('kptopt 1\n')
            fh.write('ngkpt %d %d %d\n' % tuple(param.kpts))
            fh.write('nshiftk 1\n')
            fh.write('shiftk\n')
            fh.write('%.1f %.1f %.1f\n' %
                     tuple((np.array(param.kpts) + 1) % 2 * 0.5))

        fh.write('#Definition of the SCF procedure\n')
        fh.write('toldfe %.1g\n' % param.toldfe)
        fh.write('chkexit 1 # abinit.exit file in the running directory terminates after the current SCF\n')

        fh.close()

    def read(self):
        """Read results from ABINIT's text-output file."""
        filename = self.label + '.txt'
        if not os.path.isfile(filename):
            return

        self.state = read_abinit(self.label + '.in')
        self.parameters = Parameters.read(self.label + '.param')

        self.initialize(self.state)
        self.read_results()

    def read_results(self):
        filename = self.label + '.txt'
        text = open(filename).read().lower()
        assert 'error' not in text
        assert 'was not enough scf cycles to converge' not in text

        for line in iter(text.split('\n')):
            if line.rfind('natom  ') > -1:
                natoms = int(line.split()[-1])

        lines = iter(text.split('\n'))
        # Stress:
        # Printed in the output in the following format [Hartree/Bohr^3]:
        # sigma(1 1)=  4.02063464E-04  sigma(3 2)=  0.00000000E+00
        # sigma(2 2)=  4.02063464E-04  sigma(3 1)=  0.00000000E+00
        # sigma(3 3)=  4.02063464E-04  sigma(2 1)=  0.00000000E+00
        for line in lines:
            if line.rfind(
                'cartesian components of stress tensor (hartree/bohr^3)') > -1:
                stress = np.empty(6)
                for i in range(3):
                    entries = lines.next().split()
                    stress[i] = float(entries[2])
                    stress[i + 3] = float(entries[5])
                self.results['stress'] = stress * Hartree / Bohr**3
                break
        else:
            raise RuntimeError

        # Energy [Hartree]:
        # Warning: Etotal could mean both electronic energy and free energy!
        for line in iter(text.split('\n')):
            if line.rfind('>>>>> internal e=') > -1:
                etotal = float(line.split('=')[-1])*Hartree
                for line1 in iter(text.split('\n')):
                    if line1.rfind('>>>>>>>>> etotal=') > -1:
                        efree = float(line1.split('=')[-1])*Hartree
                        break
                else:
                    raise RuntimeError
                break
        else:
            for line2 in iter(text.split('\n')):
                if line2.rfind('>>>>>>>>> etotal=') > -1:
                    etotal = float(line2.split('=')[-1])*Hartree
                    efree = self.etotal
                    break
            else:
                raise RuntimeError

        self.results['energy'] = etotal
        # Energy extrapolated to zero Kelvin:
        self.results['free energy'] = (etotal + efree) / 2

        # Forces:
        for line in lines:
            if line.rfind('cartesian forces (ev/angstrom) at end:') > -1:
                forces = []
                for i in range(natoms):
                    forces.append(np.array(
                            [float(f) for f in lines.next().split()[1:]]))
                self.results['forces'] = np.array(forces)
                break
        else:
            raise RuntimeError
        #
        self.nband = self.read_number_of_bands()
        self.niter = self.read_number_of_iterations()
        self.nelect = self.read_number_of_electrons()
        self.results['magmom'] = self.read_magnetic_moment()

    def initialize(self, atoms):
        numbers = atoms.get_atomic_numbers().copy()
        self.species = []
        for a, Z in enumerate(numbers):
            if Z not in self.species:
                self.species.append(Z)

        self.spinpol = atoms.get_initial_magnetic_moments().any()

        if 'ABINIT_PP_PATH' in os.environ:
            pppaths = os.environ['ABINIT_PP_PATH'].split(':')
        else:
            pppaths = []

        self.ppp_list = []
        if self.parameters.xc != 'LDA':
            xcname = 'GGA'
        else:
            xcname = 'LDA'

        for Z in self.species:
            symbol = chemical_symbols[abs(Z)]
            number = atomic_numbers[symbol]

            pps = self.parameters.pps
            if pps == 'fhi':
                name = '%02d-%s.%s.fhi' % (number, symbol, xcname)
            elif pps in ['paw']:
                hghtemplate = '%s-%s-%s.paw'  # E.g. "H-GGA-hard-uspp.paw"
                name = hghtemplate % (symbol, xcname, '*')
            elif pps in ['hgh.k']:
                hghtemplate = '%s-q%s.hgh.k'  # E.g. "Co-q17.hgh.k"
                name = hghtemplate % (symbol, '*')
            elif pps in ['tm']:
                hghtemplate = '%d%s%s.pspnc'  # E.g. "44ru.pspnc"
                name = hghtemplate % (number, symbol.lower(), '*')
            elif pps in ['hgh', 'hgh.sc']:
                hghtemplate = '%d%s.%s.hgh'  # E.g. "42mo.6.hgh"
                # There might be multiple files with different valence
                # electron counts, so we must choose between
                # the ordinary and the semicore versions for some elements.
                #
                # Therefore we first use glob to get all relevant files,
                # then pick the correct one afterwards.
                name = hghtemplate % (number, symbol.lower(), '*')

            found = False
            for path in pppaths:
                if (pps.startswith('paw') or
                    pps.startswith('hgh') or
                    pps.startswith('tm')):
                    filenames = glob(join(path, name))
                    if not filenames:
                        continue
                    assert len(filenames) in [0, 1, 2]
                    if pps == 'paw':
                        selector = max  # Semicore or hard
                        # warning: see download.sh in
                        # abinit-pseudopotentials*tar.gz for additional
                        # information!
                        S = selector(
                            [str(os.path.split(name)[1].split('-')[2][:-4])
                             for name in filenames])
                        name = hghtemplate % (symbol, xcname, S)
                    elif pps == 'hgh':
                        selector = min  # Lowest valence electron count
                        Z = selector([int(os.path.split(name)[1].split('.')[1])
                                      for name in filenames])
                        name = hghtemplate % (number, symbol.lower(), str(Z))
                    elif pps == 'hgh.k':
                        selector = min  # Semicore - highest electron count
                        Z = selector(
                            [int(os.path.split(name)[1].split('-')[1][:-6][1:])
                             for name in filenames])
                        name = hghtemplate % (symbol, Z)
                    elif pps == 'tm':
                        selector = max  # Semicore - highest electron count
                        # currently only one version of psp per atom
                        name = hghtemplate % (number, symbol.lower(), '')
                    else:
                        assert pps == 'hgh.sc'
                        selector = max  # Semicore - highest electron count
                        Z = selector([int(os.path.split(name)[1].split('.')[1])
                                      for name in filenames])
                        name = hghtemplate % (number, symbol.lower(), str(Z))
                filename = join(path, name)
                if isfile(filename) or islink(filename):
                    found = True
                    self.ppp_list.append(filename)
                    break
            if not found:
                raise RuntimeError('No pseudopotential for %s!' % symbol)

    def get_number_of_iterations(self):
        return self.niter

    def read_number_of_iterations(self):
        niter = None
        for line in open(self.label + '.txt'):
            if line.find(' At SCF step') != -1: # find the last iteration number
                niter = int(line.split(',')[0].split()[-1].strip())
        return niter

    def get_electronic_temperature(self):
        return self.parameters.width

    def get_number_of_electrons(self):
        return self.nelect

    def read_number_of_electrons(self):
        nelect = None
        # only in log file!
        for line in open(self.label + '.log'):  # find last one
            if line.find('with nelect') != -1:
                nelect = float(line.split('=')[1].strip())
        return nelect

    def get_number_of_bands(self):
        return self.nband

    def read_number_of_bands(self):
        nband = None
        for line in open(self.label + '.txt'): # find last one
            if line.find('     nband') != -1: # nband, or nband1, nband*
                nband = int(line.split()[-1].strip())
        return nband

    def get_kpts_info(self, kpt=0, spin=0, mode='eigenvalues'):
        return self.read_kpts_info(kpt, spin, mode)

    def get_k_point_weights(self):
        return self.get_kpts_info(kpt=0, spin=0, mode='k_point_weights')

    def get_bz_k_points(self):
        raise NotImplementedError

    def get_ibz_k_points(self):
        return self.get_kpts_info(kpt=0, spin=0, mode='ibz_k_points')

    def get_spin_polarized(self):
        return self.spinpol

    def get_number_of_spins(self):
        return 1 + int(self.spinpol)

    def get_magnetic_moment(self, atoms):
        self.update(atoms)
        return self.magnetic_moment

    def read_magnetic_moment(self):
        magmom = None
        if not self.get_spin_polarized():
            magmom = 0.0
        else: # only for spinpolarized system Magnetisation is printed
            for line in open(self.label + '.txt'):
                if line.find('Magnetisation') != -1: # last one
                    magmom = float(line.split('=')[-1].strip())
        return magmom

    def get_fermi_level(self):
        return self.read_fermi()

    def get_eigenvalues(self, kpt=0, spin=0):
        return self.get_kpts_info(kpt, spin, 'eigenvalues')

    def get_occupations(self, kpt=0, spin=0):
        return self.get_kpts_info(kpt, spin, 'occupations')

    def read_fermi(self):
        """Method that reads Fermi energy in Hartree from the output file
        and returns it in eV"""
        E_f=None
        filename = self.label + '.txt'
        text = open(filename).read().lower()
        assert 'error' not in text
        for line in iter(text.split('\n')):
            if line.rfind('fermi (or homo) energy (hartree) =') > -1:
                E_f = float(line.split('=')[1].strip().split()[0])
        return E_f*Hartree

    def read_kpts_info(self, kpt=0, spin=0, mode='eigenvalues'):
        """ Returns list of last eigenvalues, occupations, kpts weights, or
        kpts coordinates for given kpt and spin.
        Due to the way of reading output the spins are exchanged in spin-polarized case.  """
        # output may look like this (or without occupation entries); 8 entries per line:
        #
        #  Eigenvalues (hartree) for nkpt=  20  k points:
        # kpt#   1, nband=  3, wtk=  0.01563, kpt=  0.0625  0.0625  0.0625 (reduced coord)
        #  -0.09911   0.15393   0.15393
        #      occupation numbers for kpt#   1
        #   2.00000   0.00000   0.00000
        # kpt#   2, nband=  3, wtk=  0.04688, kpt=  0.1875  0.0625  0.0625 (reduced coord)
        # ...
        #
        assert mode in ['eigenvalues' , 'occupations', 'ibz_k_points', 'k_point_weights'], 'mode not in [\'eigenvalues\' , \'occupations\', \'ibz_k_points\', \'k_point_weights\']'
        if self.get_spin_polarized():
            spin = {0: 1, 1: 0}[spin]
        if spin == 0:
           spinname = ''
        else:
           spinname = 'SPIN UP'.lower()
        # number of lines of eigenvalues/occupations for a kpt
        nband = self.get_number_of_bands()
        n_entry_lines = max(1, int((nband - 0.1)/self.n_entries_float) + 1)
        #
        filename = self.label + '.txt'
        text = open(filename).read().lower()
        assert 'error' not in text
        lines = text.split('\n')
        text_list = []
        # find the begining line of last eigenvalues
        contains_eigenvalues = 0
        for n, line in enumerate(lines):
            if spin == 0:
                if line.rfind('eigenvalues (hartree) for nkpt') > -1:
                #if line.rfind('eigenvalues (   ev  ) for nkpt') > -1: #MDTMP
                    contains_eigenvalues = n
            else:
                if (line.rfind('eigenvalues (hartree) for nkpt') > -1 and
                    line.rfind(spinname) > -1): # find the last 'SPIN UP'
                        contains_eigenvalues = n
        # find the end line of eigenvalues starting from contains_eigenvalues
        text_list = [lines[contains_eigenvalues]]
        for line in lines[contains_eigenvalues + 1:]:
            text_list.append(line)
            # find a blank line or eigenvalues of second spin
            if (not line.strip() or
                line.rfind('eigenvalues (hartree) for nkpt') > -1):
                break
        # remove last (blank) line
        text_list = text_list[:-1]

        assert contains_eigenvalues, 'No eigenvalues found in the output'

        n_kpts = int(text_list[0].split('nkpt=')[1].strip().split()[0])

        # get rid of the "eigenvalues line"
        text_list = text_list[1:]

        # join text eigenvalues description with eigenvalues
        # or occupation numbers for kpt# with occupations
        contains_occupations = False
        for line in text_list:
            if line.rfind('occupation numbers') > -1:
                contains_occupations = True
                break
        if mode == 'occupations':
            assert contains_occupations, 'No occupations found in the output'

        if contains_occupations:
            range_kpts = 2*n_kpts
        else:
            range_kpts = n_kpts
        #
        values_list = []
        offset = 0
        for kpt_entry in range(range_kpts):
            full_line = ''
            for entry_line in range(n_entry_lines+1):
                full_line = full_line+str(text_list[offset+entry_line])
            first_line = text_list[offset]
            if mode == 'occupations':
                if first_line.rfind('occupation numbers') > -1:
                    # extract numbers
                    full_line = [float(v) for v in full_line.split('#')[1].strip().split()[1:]]
                    values_list.append(full_line)
            elif mode in ['eigenvalues', 'ibz_k_points', 'k_point_weights']:
                if first_line.rfind('reduced coord') > -1:
                    # extract numbers
                    if mode == 'eigenvalues':
                        full_line = [Hartree*float(v) for v in full_line.split(')')[1].strip().split()[:]]
                        #full_line = [float(v) for v in full_line.split(')')[1].strip().split()[:]] #MDTMP
                    elif mode == 'ibz_k_points':
                        full_line = [float(v) for v in full_line.split('kpt=')[1].strip().split('(')[0].split()]
                    else:
                        full_line = float(full_line.split('wtk=')[1].strip().split(',')[0].split()[0])
                    values_list.append(full_line)
            offset = offset+n_entry_lines+1
        #
        if mode in ['occupations', 'eigenvalues']:
            return np.array(values_list[kpt])
        else:
            return np.array(values_list)



def inpify(key):
    return key.lower().replace('_', '').replace('.', '').replace('-', '')


keys_with_units = {
    'tsmear': 'eV',
    'paoenergyshift': 'eV',
    'zmunitslength': 'Bohr',
    'zmunitsangle': 'rad',
    'zmforcetollength': 'eV/Ang',
    'zmforcetolangle': 'eV/rad',
    'zmmaxdispllength': 'Ang',
    'zmmaxdisplangle': 'rad',
    'ecut': 'eV',
    'dmenergytolerance': 'eV',
    'electronictemperature': 'eV',
    'oneta': 'eV',
    'onetaalpha': 'eV',
    'onetabeta': 'eV',
    'onrclwf': 'Ang',
    'onchemicalpotentialrc': 'Ang',
    'onchemicalpotentialtemperature': 'eV',
    'mdmaxcgdispl': 'Ang',
    'mdmaxforcetol': 'eV/Ang',
    'mdmaxstresstol': 'eV/Ang**3',
    'mdlengthtimestep': 'fs',
    'mdinitialtemperature': 'eV',
    'mdtargettemperature': 'eV',
    'mdtargetpressure': 'eV/Ang**3',
    'mdnosemass': 'eV*fs**2',
    'mdparrinellorahmanmass': 'eV*fs**2',
    'mdtaurelax': 'fs',
    'mdbulkmodulus': 'eV/Ang**3',
    'mdfcdispl': 'Ang',
    'warningminimumatomicdistance': 'Ang',
    'rcspatial': 'Ang',
    'kgridcutoff': 'Ang',
    'latticeconstant': 'Ang'}

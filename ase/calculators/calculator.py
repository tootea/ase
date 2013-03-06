import os
import copy
import subprocess
from math import pi, sqrt

import numpy as np


class ReadError(Exception):
    pass


# Recognized names of calculators sorted alphabetically:
names = ['abinit', 'aims', 'asap', 'castep', 'dftb', 'elk', 'emt',
         'exciting', 'fleur', 'gpaw', 'gaussian', 'hotbit', 'jacapo',
         'lammps', 'lj', 'mopac', 'morse',
         'nwchem', 'siesta', 'turbomole', 'vasp']


special = {'elk': 'ELK',
           'emt': 'EMT',
           'fleur': 'FLEUR',
           'lammps': 'LAMMPS',
           'lj': 'LennardJones',
           'morse': 'MorsePotential',
           'nwchem': 'NWChem'}


def get_calculator(name):
    """Return calculator class."""
    if name == 'asap':
        from asap3 import EMT as Calculator
    elif name == 'gpaw':
        from gpaw import GPAW as Calculator
    elif name == 'hotbit':
        from hotbit import Calculator
    else:
        classname = special.get(name, name.title())
        module = __import__('ase.calculators.' + name, {}, None, [classname])
        Calculator = getattr(module, classname)
    return Calculator


def equal(a, b):
    """ndarray-enabled comparison function."""
    if isinstance(a, np.ndarray):
        b = np.array(b)
        return a.shape == b.shape and (a == b).all()
    if isinstance(b, np.ndarray):
        return equal(b, a)
    return a == b


def kptdensity2monkhorstpack(atoms, kptdensity=3.5, even=True):
    """Convert k-point density to Monkhorst-Pack grid size.

    atoms: Atoms object
        Contains unit cell and information about boundary conditions.
    kptdensity: float
        K-point density.  Default value is 3.5 point per Ang^-1.
    even: bool
        Round to even numbers.
    """

    recipcell = atoms.get_reciprocal_cell()
    kpts = []
    for i in range(3):
        if atoms.pbc[i]:
            k = 2 * pi * sqrt((recipcell[i]**2).sum()) * kptdensity
            if even:
                kpts.append(max(1, 2 * int(round(k / 2))))
            else:
                kpts.append(max(1, int(round(k))))
        else:
            kpts.append(1)
    return kpts


def kpts2mp(atoms, kpts, even=False):
    if kpts is None:
        return (1, 1, 1)
    if isinstance(kpts, (float, int)):
        return kptdensity2monkhorstpack(atoms, kpts, even)
    else:
        return kpts


def normalize_smearing_keyword(smearing):
    """Normalize smearing string to long names and lower case."""

    smearing = smearing.lower()
    if smearing == 'fd':
        smearing = 'fermi-dirac'
    elif smearing.startswith('mp'):
        smearing = 'mathfessel-paxton-' + smearing[2]
    if smearing == 'mathfessel-paxton-0':
        smearing == 'gaussian'
    return smearing


class Parameters(dict):
    """Dictionary for parameters.
    
    Special feature: If param is a Parameters instance, then param.xc
    is a shorthand for param['xc'].
    """
    
    def __getattr__(self, key):
        return self[key]

    def __setattr__(self, key, value):
        self[key] = value

    @classmethod
    def read(cls, filename):
        """Read parameters from file."""
        file = open(os.path.expanduser(filename))
        parameters = cls(eval(file.read()))
        file.close()
        return parameters

    def tostring(self):
        keys = sorted(self.keys())
        return 'dict(' + ',\n     '.join(
            '%s=%r' %(key, self[key]) for key in keys) + ')\n'
    
    def write(self, filename):
        file = open(filename, 'w')
        file.write(self.tostring())
        file.close()


class Calculator:
    """Base-class for all ASE calculators."""

    notimplemented = []
    "Properties calculator can't handle"

    default_parameters = {}
    'Default parameters'

    def __init__(self, restart=None, ignore_bad_restart_file=False, label=None,
                 atoms=None, **kwargs):
        """Basic calculator implementation.

        restart: str
            Prefix for restart file.  May contain a directory.  Default
            is None: don't restart.
        ignore_bad_restart_file: bool
            Ignore broken or missing restart file.  By defauls, it is an
            error if the restart file is missing or broken.
        label: str
            Name used for all files.  May contain a directory.
        atoms: Atoms object
            Optional Atoms object to which the calculator will be
            attached.  When restarting, atoms will get its positions and
            unit-cell updated from file.
        """

        self.state = None  # copy of atoms object from last calculation
        self.results = {}  # calculated properties (energy, forces, ...)
        self.parameters = None  # calculational parameters

        if restart is not None:
            try:
                self.read(restart)  # read parameters, state and results
            except ReadError:
                if ignore_bad_restart_file:
                    self.reset()
                else:
                    raise
        
        self.label = None
        self.directory = None
        self.prefix = None

        self.set_label(label)
        
        if self.parameters is None:
            # Use default parameters if they were not read from file: 
            self.parameters = self.get_default_parameters()

        if atoms is not None:
            atoms.calc = self
            if self.state is not None:
                # State was read from file.  Update atoms:
                if not (equal(atoms.numbers, self.state.numbers) and
                        (atoms.pbc == self.state.pbc).all()):
                    raise RuntimeError('Atoms not compatible with file')
                atoms.positions = self.state.positions
                atoms.cell = self.state.cell
                
        self.set(**kwargs)

        self.callbacks = {'before': [], 'after': []}  # call-back functions

        if not hasattr(self, 'name'):
            self.name = self.__class__.__name__

    def set_label(self, label):
        """Set label and convert label to directory and prefix.

        Examples:

        * label='abc': (directory='.', prefix='abc')
        * label='dir1/abc': (directory='dir1', prefix='abc')

        Calculators that must write results to files with fixed names
        can overwrite this method so that the directory is set to all
        of label."""

        self.label = label

        if label is None:
            self.directory = None
            self.prefix = None
        else:
            self.directory, self.prefix = os.path.split(label)
            if self.directory == '':
                self.directory = os.curdir

    def get_default_parameters(self):
        return Parameters(copy.deepcopy(self.default_parameters))

    def todict(self):
        data = copy.deepcopy(self.parameters)
        data['name'] = self.name
        return data

    def reset(self):
        """Clear all information from old calculation."""

        self.state = None
        self.results = {}

    def read(self, label):
        """Read atoms, parameters and calculated properties from output file.

        Read result from self.label file.  Raise ReadError if the file
        is not there.  If the file is corrupted or contains an error
        message from the calculation, a ReadError should also be
        raised.  In case of succes, these attributes must set:

        state: Atoms object
            The state of the atoms from last calculation.
        parameters: Parameters object
            The parameter dictionary.
        results: dict
            Calculated properties like energy and forces.

        The FileIOCalculator.read() method will typically read state
        and parameters and get the results dict by calling the
        read_results() method."""

        self.set_label(label)

    def get_atoms(self):
        if self.state is None:
            raise ValueError('Calculator has no atoms')
        atoms = self.state.copy()
        atoms.calc = self
        return atoms

    @classmethod
    def read_atoms(cls, restart, **kwargs):
        return cls(restart, **kwargs).get_atoms()

    def set(self, **kwargs):
        """Set parameters like set(key1=value1, key2=value2, ...).
        
        A dictionary containing the parameters that have been changed
        is returned.

        Subclasses must implement a set() method that will look at the
        chaneged parameters and decide if a call to reset() is needed.
        If the changed parameters are harmless, like a change in
        verbosity, then there is no need to call reset().

        The special keyword 'parameters' can be used to read
        parameters from a file."""

        if 'parameters' in kwargs:
            filename = kwargs.pop('parameters')
            parameters = Parameters.read(filename)
            parameters.update(kwargs)
            kwargs = parameters

        changed_parameters = {}

        for key, value in kwargs.items():
            oldvalue = self.parameters.get(key)
            if key not in self.parameters or not equal(value, oldvalue):
                if isinstance(oldvalue, dict):
                    # Special treatment for dictionary parameters:
                    for name in value:
                        if name not in oldvalue:
                            raise KeyError(
                                'Unknown subparameter "%s" in '
                                'dictionary parameter "%s"' % (name, key))
                    oldvalue.update(value)
                    value = oldvalue
                changed_parameters[key] = value
                self.parameters[key] = value

        return changed_parameters

    def check_state(self, atoms):
        """Check for system changes since last calculation."""
        if self.state is None:
            system_changes = ['positions', 'numbers', 'cell', 'pbc', 'magmoms']
        else:
            system_changes = []
            if not equal(self.state.positions, atoms.positions):
                system_changes.append('positions')
            if not equal(self.state.numbers, atoms.numbers):
                system_changes.append('numbers')
            if not equal(self.state.cell, atoms.cell):
                system_changes.append('cell')
            if not equal(self.state.pbc, atoms.pbc):
                system_changes.append('pbc')
            if not equal(self.state.get_initial_magnetic_moments(),
                         atoms.get_initial_magnetic_moments()):
                system_changes.append('magmoms')

        return system_changes

    def get_potential_energy(self, atoms, force_consistent=False):
        energy = self.get_property('energy', atoms)
        if force_consistent:
            return self.results.get('free energy', energy)
        else:
            return energy

    def get_forces(self, atoms):
        return self.get_property('forces', atoms).copy()

    def get_stress(self, atoms):
        return self.get_property('stress', atoms).copy()

    def get_dipole_moment(self, atoms):
        return self.get_property('dipole', atoms).copy()

    def get_magnetic_moment(self, atoms):
        return self.get_property('magmom', atoms)

    def get_magnetic_moments(self, atoms):
        return self.get_property('magmoms', atoms).copy()

    def get_property(self, name, atoms):
        if name in self.notimplemented:
            raise NotImplementedError

        system_changes = self.check_state(atoms)
        if system_changes:
            self.reset()

        if name not in self.results:
            self._calculate(atoms, [name], system_changes)
        return self.results[name]

    def calculation_required(self, atoms, properties):
        system_changes = self.check_state(atoms)
        if system_changes:
            return True
        for name in properties:
            if name not in self.results:
                return True
        return False
        
    def _calculate(self, atoms, properties, system_changes):
        """Call callbacks before and after actual calculation."""
        self.call_callbacks('before')
        self.state = atoms.copy()
        try:
            self.calculate(atoms, properties, system_changes)
        except Exception:
            self.reset()
            raise
        self.call_callbacks('after')

    def calculate(self, atoms, properties, system_changes):
        """Do the calculation.

        atoms: Atoms object
            Contains positions, unit-cell, ...
        properties: list of str
            List of what needs to be calculated.  Can be any combination
            of 'energy', 'forces', 'stress', 'dipole', 'magmom' and
            'magmoms'.
        system_changes: list of str
            List of what has changed since last calculation.  Can be
            any combination of these five: 'positons', 'numbers', 'cell',
            'pbc' and 'magmoms'.

        Subclasses need to implement this, but can ignore properties
        and system_changes if they want.
        """

        # Dummy calculation:
        self.results = {'energy': 0.0,
                        'forces': np.zeros((len(atoms), 3)),
                        'stress': np.zeros(6),
                        'dipole': np.zeros(3),
                        'magmom': 0.0,
                        'magmoms': np.zeros(len(atoms))}
                        
    def attach_callback(self, name, function, *args, **kwargs):
        """Attach function to be called before or after calculation.

        name: str
            One of 'before' or 'after'
        function: callable
            Function or callable to call.
        """

        self.callbacks[name].append((function, args, kwargs))

    def call_callbacks(self, name):
        for function, args, kwargs in self.callbacks[name]:
            function(*args, **kwargs)

    def get_spin_polarized(self):
        return False


class FileIOCalculator(Calculator):
    "Base class for calculators that write input files and read output files."

    command = None
    'Command used to start calculation'

    def __init__(self, restart=None, ignore_bad_restart_file=False,
                 label=None, atoms=None, command=None, **kwargs):
        """File-IO calculator.
        
        command: str
            Command used to start calculation.
        """

        Calculator.__init__(self, restart, ignore_bad_restart_file, label,
                            atoms, **kwargs)

        if command is not None:
            self.command = command
        else:
            name = 'ASE_' + self.name.upper() + '_COMMAND'
            self.command = os.environ.get(name, self.command)

    def calculate(self, atoms, properties=None, system_changes=None):
        self.write_input(atoms, properties, system_changes)
        if self.command is None:
            raise RuntimeError('Please set $%s environment variable ' %
                               ('ASE_' + self.name.upper() + '_COMMAND') +
                               'or supply the command keyword')
        command = self.command.replace('PREFIX', self.prefix)
        olddir = os.getcwd()
        try:
            os.chdir(self.directory)
            errorcode = subprocess.call(command, shell=True)
        finally:
            os.chdir(olddir)
        
        if errorcode:
            raise RuntimeError('%s returned an error: %d' %
                               (self.name, errorcode))
        self.read_results()

    def write_input(self, atoms, properties=None, system_changes=None):
        """Write input file(s).

        Call this method first in subclasses so that directories are
        created automatically."""

        if self.directory != os.curdir and not os.path.isdir(self.directory):
            os.makedirs(self.directory)

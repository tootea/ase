"""
A simple interface for ADF by Tomas Trnka <ttrnka@mail.muni.cz>
"""
import os
import os.path
import re
import shutil
import subprocess
import sys
import time

if "ADFHOME" in os.environ:
    sys.path.append(os.environ["ADFHOME"] + "/scripting")

import kf

import ase.io
from ase.units import kcal, mol, Hartree, Bohr

class ADF:
    name = 'ADF'
    def __init__(self, label='asejob', **kwargs):
        self.init_input = None
        self.workdir = None
        self.restart_file = None
        self.xyz_file = "current.xyz"
        self.command = "./{0:s} > {0:s}.stdout 2>&1 &"

        # save label
        self.label = label

        #set atoms
        self.atoms = None

        self.input_lines = None
        self.charge_level = "d"
        self.qm_atoms = set()
        self.active_atoms = set()
        self.first_run = True
        self.running = False
        self.tape21 = None
        self.energy = None
        self.forces = None
        self.atom_order_index = None
        self.iteration = 0

        # set user values
        self.set(**kwargs)

        if self.charge_level == "m":
            self.charge_level_id = 0
        elif self.charge_level == "d":
            self.charge_level_id = 1
        elif self.charge_level == "q":
            self.charge_level_id = 2
        else:
            raise RuntimeError("Invalid charge level")

        self.initialize()

    def set(self, **kwargs):
        """
        Sets the parameters on the according keywords
        Raises RuntimeError when wrong keyword is provided
        """
        for key in kwargs:
            if kwargs[key] is None:
                continue
            elif key == "input":
                self.init_input = kwargs[key]
            elif key == "restart":
                self.restart_file = os.path.abspath(kwargs[key])
            elif key == "xyz":
                self.xyz_file = kwargs[key]
            elif key == "workdir":
                self.workdir = os.path.abspath(kwargs[key])
            elif key == "active":
                self.active_atoms = set(kwargs[key])
            elif key == "command":
                self.command = kwargs[key]
            else:
                raise RuntimeError('ADF calculator: unknown keyword: ' + key)

        if self.workdir is None:
            self.workdir = "./"
        elif not self.workdir.endswith("/"):
            self.workdir += "/"

    def set_positions(self, pos):
        if (abs(pos - self.atoms.get_positions()) > 1e-6).any():
            self.atoms.set_positions(pos)
            self.run()

    def set_atoms(self, atoms):
        if self.atoms is None or (len(atoms) != len(self.atoms)):
            self.atoms = atoms.copy()
            self.run()
        else:
            self.set_positions(atoms.get_positions())

    def initialize(self):
        self.read_input(self.init_input)
        self.read_input_charges()
        self.read_qm_atoms()

    def read_input(self, fname):
        input_file = open(fname, "r")
        self.input_lines = input_file.readlines()
        input_file.close()

    def write_input(self, fname):
        new_input = open(fname, "w")
        new_input.writelines(self.input_lines)
        os.fchmod(new_input.fileno(), 0755)
        new_input.close()

    def fixup_input(self):
        i = 0
        first_line = None
        skipping = False
        while i < len(self.input_lines):
            line = self.input_lines[i].strip()
            if (skipping):
                if "END" in line:
                    skipping = False
                del self.input_lines[i]
            elif line.startswith("TITLE"):
                del self.input_lines[i]
            elif line.startswith("GEOMETRY"):
                skipping = True
                del self.input_lines[i]
            elif line.startswith("GRADIENT"):
                del self.input_lines[i]
            elif line.startswith("RESTART"):
                if line.endswith("&"):
                    skipping = True
                del self.input_lines[i]
            elif line.startswith("touch"):
                del self.input_lines[i]
            elif (first_line is None) and ("<<" in line):
                first_line = i + 1
            else:
                i += 1

        init_cmds = [
            "TITLE ASE ADF job: {:s}\n".format(self.label),
            "GRADIENT\n"
            ]

        if self.restart_file is not None:
            init_cmds += [
            "RESTART {:s} &\n".format(self.restart_file),
            " noGEO\n",
            "END\n"
            ]

        self.input_lines[first_line:first_line] = init_cmds

        self.input_lines += [
            "\n",
            "> FINISHED\n"
            ]

    def read_qm_atoms(self):
        self.qm_atoms.clear()
        in_conntable = False
        current_line_iter = iter(self.input_lines)
        for line in current_line_iter:
            if "MM_CONNECTION_TABLE" in line:
                in_conntable = True
                break

        if not in_conntable:
            return

        for line in current_line_iter:
            if "SUBEND" in line:
                break

            fields = line.split()
            if fields[2] == "QM":
                self.qm_atoms.add(int(fields[0]))

    def update_coordinates(self, atoms):
        start_line = 0
        atoms_found = False
        for start_line in xrange(0, len(self.input_lines)):
            if "ATOMS" in self.input_lines[start_line]:
                atoms_found = True
                break

        if not atoms_found:
            raise RuntimeError('ATOMS block not found')

        start_line += 1

        for end_line in xrange(start_line, len(self.input_lines)):
            if "END" in self.input_lines[end_line]:
                break

        xyz_atoms = ase.io.read(self.xyz_file, 0, "xyz")

        current_atom = iter(atoms)
        atom_lines = list()
        for iat in xrange(1, len(xyz_atoms) + 1):
            if iat in self.active_atoms:
                pos = current_atom.next().position
                xyz_atoms[iat - 1].position = pos
            else:
                pos = xyz_atoms[iat - 1].position

            atom_lines.append("{:6d} {:<3s} {:9.4f} {:9.4f} {:9.4f}\n".format
                              (iat, xyz_atoms[iat - 1].symbol, *pos))

        self.input_lines[start_line:end_line] = atom_lines

        ase.io.write(self.workdir + str(self.iteration) + ".xyz", xyz_atoms)
        self.iteration += 1

    def read_input_charges(self):
        self.charges = None
        in_charges = False
        current_line_iter = iter(self.input_lines)
        for line in current_line_iter:
            if "CHARGES" in line:
                in_charges = True
                break

        if not in_charges:
            return

        self.charges = dict()

        for line in current_line_iter:
            if "END" in line:
                break

            fields = line.split()

            self.charges[int(fields[0])] = float(fields[1])

    def read_charges(self, fname):
        f = open(fname, "r")

        for line in f:
            if "Multipole derived atomic charges" in line:
                break
        for line in f:
            if re.search("MDC-m +MDC-d +MDC-q *$", line):
                break
        f.next()

        for line in f:
            fields = line.split()
            if len(fields) != 5:
                break
            iat = int(fields[0])
            if iat in self.qm_atoms:
                self.charges[iat] = float(fields[2 + self.charge_level_id])

        for line in f:
            m = re.match("^ *Added capping atom *([0-9]+) charge *(-*[0-9\.]+) to atom *([0-9]+) with new value *(-*[0-9\.]+)", line)
            if m:
                icap = int(m.group(1))
                ireal = int(m.group(3))
                qreal = float(m.group(4))

                if icap in self.charges:
                    del self.charges[icap]
                self.charges[ireal] = qreal

    def update_charges(self):
        if self.charges is None:
            return

        charge_lines = list()
        for i in self.charges:
            charge_lines.append("{:9d} {:8.4f}\n".format(i, self.charges[i]))

        start_line = 0
        for start_line in xrange(0, len(self.input_lines)):
            if "CHARGES" in self.input_lines[start_line]:
                break

        start_line += 1

        for end_line in xrange(start_line, len(self.input_lines)):
            if "END" in self.input_lines[end_line]:
                break

        self.input_lines[start_line:end_line] = charge_lines

    def run(self):
        self.prepare_workdir()
        if not self.first_run:
            self.xyz_file = self.workdir + "prev/current.xyz"
            self.restart_file = self.workdir + "prev/TAPE21"

        self.update_coordinates(self.atoms)
        self.update_charges()
        self.fixup_input()

        self.write_input(self.workdir + "current/" + self.label)
        self.first_run = False

        os.system("cd " + self.workdir + "current/ && " + self.command.format(self.label))
        self.running = True

    def prepare_workdir(self):
        if self.tape21 is not None:
            self.tape21.close()
            self.tape21 = None

        if not os.path.exists(self.workdir):
            os.mkdir(self.workdir, 0755)

        if os.path.exists(self.workdir + "prev"):
            shutil.rmtree(self.workdir + "prev/")

        if os.path.exists(self.workdir + "current"):
            os.rename(self.workdir + "current/", self.workdir + "prev/")

        os.mkdir(self.workdir + "current", 0755)


    def finish_run(self):
        if self.tape21 is not None:
            return

        if not self.running:
            self.run()

        stdout_path = self.workdir + "current/{:s}.stdout".format(self.label)
        t21_path = self.workdir + "current/TAPE21"

        while not (os.path.exists(self.workdir + "current/FINISHED") and
                   os.path.exists(stdout_path) and
                   os.path.exists(t21_path)
                   ):
            time.sleep(1)

        self.running = False

        saved_stdout = self.workdir + str(self.iteration) + ".out"
        if (os.path.exists(saved_stdout)):
            os.remove(saved_stdout)
        os.link(stdout_path, saved_stdout)

        self.read_charges(stdout_path)
        self.tape21 = kf.kffile(t21_path)
        self.energy = None
        self.forces = None

    def get_potential_energy(self, atoms=None):
        if atoms is not None:
            self.set_atoms(atoms)
        self.finish_run()

        if self.energy is None:
            self.energy = float(self.tape21.read("Energy", "Bond Energy")) * Hartree
            saved_energy = open(self.workdir + str(self.iteration) + ".E", "w")
            saved_energy.write("{:f}\n".format(self.energy))
            saved_energy.close()
        return self.energy

    def from_internal_order(self, a):
        if self.atom_order_index is None:
            self.atom_order_index = list()
            index_tmp = self.tape21.read("Geometry", "atom order index")
            for i in xrange(1, len(index_tmp) / 2 + 1):
                if i in self.active_atoms:
                    self.atom_order_index.append(index_tmp[i - 1] - 1)

        return a.take(self.atom_order_index, 0)

    def get_forces(self, atoms=None):
        if atoms is not None:
            self.set_atoms(atoms)
        self.finish_run()

        if self.forces is None:
            grad_internal = self.tape21.read("GeoOpt", "Gradients_CART").reshape((-1, 3))
            self.forces = -self.from_internal_order(grad_internal) * Hartree / Bohr

        return self.forces

    def get_stress(self, atoms=None):
        raise NotImplementedError

    def calculation_required(atoms, quantities):
        if self.positions_changed(atoms) or self.tape21 is None:
            return True

        for qty in quantities:
            if qty != "energy" and qty != "forces":
                return True

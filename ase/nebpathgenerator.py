#!/usr/bin/python

import ase.atoms
import ase.calculators.adf
import math
import numpy
import scipy.interpolate

class NEBPathGenerator(object):
    def __init__(self):
        self.input_points = list()

    def add_point(self, atoms, input_, xyz, restart):
        self.input_points.append((atoms, input_, xyz, restart))

    def generate(self, N, **kwargs):
        images = list()
        M = len(self.input_points)

        for i in xrange(0, N):
            relpos = float(i) * (M - 1) / (N - 1)
            iclosest = int(round(relpos))

            iprev = int(math.floor(relpos))
            inext = int(math.ceil(relpos))

            print i, relpos, iclosest, iprev, inext

            (src_at, src_inp, src_xyz, src_rst) = self.input_points[iclosest]

            im = src_at.copy()
            if iprev != inext:
                diff = self.input_points[inext][0].positions - self.input_points[iprev][0].positions
                im.set_positions(self.input_points[iprev][0].positions + (relpos - iprev) * diff)

            if i != 0 and i != N - 1:
                im.set_calculator(ase.calculators.adf.ADF(input=src_inp,
                    xyz=src_xyz, restart=src_rst, workdir=str(i), **kwargs))

            images.append(im)

        return images

class SplineNEBPathGenerator(NEBPathGenerator):
    def __init__(self, tol=0.0):
        super(SplineNEBPathGenerator, self).__init__()
        self.tol = tol

    def interpolate(self, x, u, u_interpolated):
        (m, dim) = x.shape
        n = len(u_interpolated)
        s = self.tol**2 * m
        x_interpolated = numpy.empty((n, dim))

        if m > 4:
            k = 3
        elif m >= 2:
            k = m - 1
        else:
            raise ValueError('at least two points needed')

        for i in xrange(0, dim):
            spl = scipy.interpolate.UnivariateSpline(u, x[:, i], k=k, s=s)
            x_interpolated[:, i] = spl(u_interpolated)

        return x_interpolated


    def generate(self, N, **kwargs):
        images = list()
        M = len(self.input_points)

        x = numpy.array([pt[0].positions.flatten() for pt in self.input_points])

        u = numpy.empty(M)
        u[0] = 0.0
        for i in xrange(1, M):
            u[i] = u[i - 1] + numpy.linalg.norm(x[i] - x[i - 1])
        u /= u[-1]

        u_interpolated = numpy.linspace(0.0, 1.0, num=N, endpoint=True)
        x_interpolated = self.interpolate(x, u, u_interpolated)
        x_interpolated.shape = (N, -1, 3)

        for i in xrange(0, N):
            i_closest = numpy.abs(u - u_interpolated[i]).argmin()

            print i, u_interpolated[i], i_closest

            (src_at, src_inp, src_xyz, src_rst) = self.input_points[i_closest]

            im = src_at.copy()
            im.set_positions(x_interpolated[i])

            if i != 0 and i != N - 1:
                im.set_calculator(ase.calculators.adf.ADF(input=src_inp,
                    xyz=src_xyz, restart=src_rst, workdir=str(i), **kwargs))

            images.append(im)

        return images


class NonuniformLinear(NEBPathGenerator):
    def generate(self, maxmov, **kwargs):
        images = list()
        M = len(self.input_points)

        k = 0
        for i in xrange(0, M - 1):
            diff = self.input_points[i + 1][0].positions - self.input_points[i][0].positions
            d = max([numpy.linalg.norm(v) for v in diff])
            num_interpol = int(math.floor(d / maxmov))
            print i, d, num_interpol, ":"

            for j in xrange(0, num_interpol + 1):
                relpos = float(j) / (num_interpol + 1)
                iclosest = i + int(round(relpos))

                (src_at, src_inp, src_xyz, src_rst) = self.input_points[iclosest]

                print " ", k, relpos, iclosest

                im = src_at.copy()

                if j != 0:
                    im.set_positions(self.input_points[i][0].positions + relpos * diff)

                if k != 0:
                    im.set_calculator(ase.calculators.adf.ADF(input=src_inp,
                        xyz=src_xyz, restart=src_rst, workdir=str(k), **kwargs))

                images.append(im)
                k += 1

        images.append(self.input_points[-1][0].copy())

        return images
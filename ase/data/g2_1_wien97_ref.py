# atomization energies in kcal / mol (= 43.364 meV)
# All values evaluated with PBE xc-orbitals and densities at
# experimental geometries. Zero-point vibration has been removed
# from experimental energies (from [1]).
atomization = {
# Molec   expt    LSD    PBE    RPBE   BLYP
'H2'  : ( 109.5, 113.2, 104.6, 105.5, 109.4),
'LiH' : (  57.8,  61.0,  53.5,  53.4,  58.1),
'CH4' : ( 419.3, 462.3, 419.8, 410.6, 416.6),
'NH3' : ( 297.4, 337.3, 301.7, 293.2, 301.4),
'OH'  : ( 106.4, 124.1, 109.8, 106.3, 109.6),
'H2O' : ( 232.2, 266.5, 234.2, 226.6, 232.5),
'HF'  : ( 140.8, 162.2, 142.0, 137.5, 141.0),
'Li2' : (  24.4,  23.9,  19.9,  20.2,  20.5),
'LiF' : ( 138.9, 156.1, 138.6, 132.9, 140.1),
'Be2' : (   3.0,  12.8,   9.8,   7.9,   6.1),
'C2H2': ( 405.4, 460.3, 414.9, 400.4, 405.3),
'C2H4': ( 562.6, 632.6, 571.5, 554.5, 560.7),
'HCN' : ( 311.9, 361.0, 326.1, 313.6, 320.3),
'CO'  : ( 259.3, 299.1, 268.8, 257.9, 261.8),
'N2'  : ( 228.5, 267.4, 243.2, 232.7, 239.8),
'NO'  : ( 152.9, 198.7, 171.9, 161.6, 166.0),
'O2'  : ( 120.5, 175.0, 143.7, 133.3, 135.3),
'F2'  : (  38.5,  78.2,  53.4,  45.6,  49.4),
'P2'  : ( 117.3, 143.8, 121.1, 114.1, 121.0),
'Cl2' : (  58.0,  83.0,  65.1,  58.9,  57.2)
}

# exchange-only atomization energies in kcal / mol (= 43.364 meV)
# All values evaluated with PBE xc-orbitals and densities at
# experimental geometries. (from [1]).
ex_atomization = {
# Molec   exact   LSD    PBE    RPBE   BLYP
'H2'  : (  84.0,  81.5,  84.8,  85.8,  85.4),
'LiH' : (  33.9,  33.6,  36.9,  36.8,  36.2),
'CH4' : ( 327.2, 369.9, 336.0, 326.9, 331.2),
'NH3' : ( 199.5, 255.0, 227.4, 218.9, 222.6),
'OH'  : (  67.3,  96.2,  84.5,  80.9,  82.7),
'H2O' : ( 154.6, 212.9, 183.9, 176.4, 180.5),
'HF'  : (  96.1, 136.1, 117.1, 112.6, 115.4),
'Li2' : (   3.5,   6.5,   6.4,   6.7,   3.9),
'LiF' : (  86.8, 129.7, 116.5, 110.8, 113.6),
'Be2' : ( -11.0,   9.4,   3.1,   1.2,   1.6),
'C2H2': ( 290.6, 382.7, 333.0, 318.5, 325.5),
'C2H4': ( 423.9, 517.7, 456.5, 439.6, 447.4),
'HCN' : ( 194.5, 294.0, 256.1, 243.5, 249.1),
'CO'  : ( 169.2, 261.9, 224.0, 213.1, 218.7),
'N2'  : ( 110.2, 211.4, 184.1, 173.6, 177.6),
'NO'  : (  45.6, 156.9, 122.8, 112.5, 117.0),
'O2'  : (  24.9, 147.5, 104.4,  94.1,  99.3),
'F2'  : ( -43.3,  64.0,  32.5,  24.7,  28.8),
'P2'  : (  31.8,  98.4,  73.1,  66.1,  70.1),
'Cl2' : (  15.5,  68.2,  39.8,  33.7,  37.0)
}

# Exchange energy of some spherical atoms in Hartrees (= 27.211 eV)
# All functionals were evaluated with self-consistent exchange-only
# OEP orbitals and densities (from [1]).
ex_energy = {
# Atom       exact     LSD       PBE       RPBE      BLYP
'H'   : (   0.3125,   0.2680,   0.3059,   0.3112,   0.3098),
'He'  : (   1.0258,   0.8840,   1.0136,   1.0313,   1.0255),
'Li'  : (   1.7807,   1.5379,   1.7572,   1.7876,   1.7753),
'Be'  : (   2.6658,   2.3124,   2.6358,   2.6801,   2.6578),
'N'   : (   6.6044,   5.9008,   6.5521,   6.6252,   6.5961),
'Ne'  : (  12.1050,  11.0335,  12.0667,  12.1593,  12.1378),
'Na'  : (  14.0131,  12.7859,  13.9506,  14.0528,  14.0304),
'Mg'  : (  15.9884,  14.6117,  15.9147,  16.0260,  16.0005),
'P'   : (  22.6341,  20.7931,  22.5028,  22.6369,  22.6221),
'Ar'  : (  30.1747,  27.8632,  29.9961,  30.1494,  30.1535),
'Kr'  : (  93.8330,  88.6245,  93.4257,  93.6645,  93.8721),
'Xe'  : ( 179.0635, 170.5660, 178.2450, 178.5649, 179.0427)
}

# Correlation energy of some spherical atoms in Hartrees (= 27.211 eV)
# All functionals were evaluated with self-consistent exchange-only
# OEP orbitals and densities (from [1]).
# 'Exact' values are from reference [2]
ec_energy = {
# Atom     exact   LSD     PBE     BLYP
'H'   : ( 0.0000, 0.0222, 0.0060, 0.0000),
'He'  : ( 0.0420, 0.1125, 0.0420, 0.0438),
'Li'  : ( 0.0455, 0.1508, 0.0514, 0.0534),
'Be'  : ( 0.0950, 0.2240, 0.0856, 0.0945),
'N'   : ( 0.1858, 0.4268, 0.1799, 0.1919),
'Ne'  : ( 0.3929, 0.7428, 0.3513, 0.3835),
'Na'  : ( 0.3988, 0.8010, 0.3715, 0.4083),
'Mg'  : ( 0.4424, 0.8874, 0.4110, 0.4594),
'P'   : ( 0.5446, 1.1127, 0.5265, 0.5664),
'Ar'  : ( 0.7314, 1.4242, 0.7067, 0.7508),
'Kr'  : ( 3.2693, 1.7672, 1.7486, 2.0788),
'Xe'  : ( 5.1773, 2.9184, 2.7440, 3.1789)
}

def convert(input, column):
    keys = input.keys()
    keys.sort()
    data = {}
    for k in keys:
        data[k] = input[k][column]
    return data

info = {}

info['atomization energy'] = {}

info['atomization energy'].update({'reference': convert(atomization, 0)})
info['atomization energy'].update({'LSD': convert(atomization, 1)})
info['atomization energy'].update({'PBE': convert(atomization, 2)})
info['atomization energy'].update({'RPBE': convert(atomization, 3)})
info['atomization energy'].update({'BLYP': convert(atomization, 4)})

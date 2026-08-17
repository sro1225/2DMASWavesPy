#
#    2DMASWavesPy, a 2D implemenmtation of mwaswavespy, a Python package
#    for processing and inverting MASW data
#    Original code Copyright (C) 2023  Elin Asta Olafsdottir (elinasta(at)hi.is)
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.

from setuptools import setup, Extension, find_packages
from Cython.Build import cythonize
import numpy as np

extensions = [
    Extension(
        "maswavespy.cy_dispersion_imaging",
        ["src/maswavespy/cy_dispersion_imaging.pyx"],
        include_dirs=[np.get_include()],
    ),
    Extension(
        "maswavespy.cy_theoretical_dc",
        ["src/maswavespy/cy_theoretical_dc.pyx"],
        include_dirs=[np.get_include()],
    ),
]

setup(
    name="maswavespy",
    packages=find_packages(where="src"),
    package_dir={"": "src"},

    install_requires=[
        "numpy>=1.24",
        "pandas>=2.0",
        "matplotlib>=3.8",
        "scipy>=1.11",
        "cython",
        "obspy>=1.4",
        "pyyaml>=6.0",
    ],

    ext_modules=cythonize(
        extensions,
        language_level="3"
    ),
    include_dirs=[np.get_include()],
)


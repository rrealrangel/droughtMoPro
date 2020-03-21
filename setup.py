from setuptools import setup
from os import path

this_directory = path.abspath(path.dirname(__file__))

with open(path.join(this_directory, 'README.md')) as f:
    long_description = f.read()

def read(*rnames):
    return open(path.join(path.dirname(__file__), *rnames)).read()

setup(
     name='droughtMoPro',
     version='v0.1.4',
     author='R. A. Real-Rangel',
     author_email='rrealr@iingen.unam.mx',
     description='Generator of drought monitoring products.',
     license='GPL-3.0',
     keywords="drought monitoring mapping series",
     url='https://github.com/rrealrangel/droughtMoPro',
     packages=['droughtMoPro'],
     long_description=long_description,
     long_description_content_type='text/markdown',
     classifiers=[
        "Development Status :: 4 - Beta",
        "Topic :: Scientific/Engineering :: Hydrology",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)"
        ],
     install_requires=[
         'datetime',
         'dask',
         'netcdf4',
         'numpy',
         'openpyxl',
         'pandas',
         'pathlib2',
         'scipy',
         'toml',
         'toolz',
         'xarray',
         ],
     zip_safe=False
     )

from setuptools import setup

setup(
     name='droughtMoPro',
     version='v0.1.3.2',
     author='R. A. Real-Rangel',
     author_email='rrealr@iingen.unam.mx',
     description='Generator of drought monitoring products.',
     license='GPL-3.0',
     keywords="drought monitoring mapping series",
     url='https://github.com/rrealrangel/droughtMoPro',
     packages=['droughtMoPro'],
     long_description=read('README'),
     classifiers=[
        "Development Status :: 3 - Alpha",
        "Topic :: Scientific/Engineering :: Hydrology",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)"]
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

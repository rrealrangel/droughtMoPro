from setuptools import setup

setup(
     name='droughtMoPro',
     version='v0.1.3.2',
     description='Generator of drought monitoring products.',
     url='https://github.com/rrealrangel/droughtMoPro',
     author='R. A. Real-Rangel',
     author_email='rrealr@iingen.unam.mx',
     license='GPL-3.0',
     packages=[
         'droughtMoPro'
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

from setuptools import setup, find_packages

#https://godatadriven.com/blog/a-practical-guide-to-using-setup-py/
setup(
    name='MJOcast',
    version='0.1.0',
    author='William E. Chapman',
    author_email='wchapman@ucar.edu',
    description='A package to create the forecasting of MJO indices',
    license='GNU General Public License v2',
    packages=find_packages(include='MJOcast*'),
    install_requires=[
        'pyyaml',
        'pandas',#pandas 2.0? 
        'numpy',
        'matplotlib',
        'xarray',
        'eofs',
        'datetime',
        'scipy',
        'netCDF4'
    ],
    extras_require={
        'full_func': ['jupyter','dask'],
        'dev': ['build', 'pytest', 'pytest-pep8','coverage'],
      },
    package_data={'MJOcast': ['*/*.nc']},
)

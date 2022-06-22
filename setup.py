import setuptools
from setuptools import setup
from glob import glob


setup(
    name='seq_cluster',
    version='0.8',
    packages=['seq_cluster', 'seq_cluster/clustering', 'seq_cluster/data', 'seq_cluster/distances',
              'seq_cluster/distances/matrices', 'seq_cluster/main', 'seq_cluster/visualization'],
    package_data={'seq_cluster/distances': ['matrices/*']},
    install_requires = ['numpy==1.22.0'],
    url='',
    license='GNU General Public License v3.0',
    author='Ralf Stefan Neumann',
    author_email='',
    description=''
)

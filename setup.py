import os
import sys
import shutil
from setuptools import setup, find_packages, find_namespace_packages
from setuptools.command.install import install
#from distutils.core import setup

#NAME = 'flexce'
# do not use x.x.x-dev.  things complain.  instead use x.x.xdev
#VERSION = '1.0.1dev'
#RELEASE = 'dev' not in VERSION

setup(name='flexce',
      version='1.0.1',
      description='Flexible Galactic Chemical Evolution Model',
      author='Brett Andrews',
      author_email='brett.h.andrews@gmail.com',
      url='https://github.com/dnidever/flexce',
      requires=['numpy','astropy(>=4.0)','scipy','dlnpyutils','doppler'],
      zip_safe = False,
      include_package_data=True,
      packages=find_namespace_packages(where="python"),
      package_dir={"": "python"}
)

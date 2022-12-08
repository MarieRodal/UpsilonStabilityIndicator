# -*- coding: utf-8 -*-
"""
Created on Sun Nov 15 14:39:28 2020

@author: Marie Rodal 
"""

from setuptools import setup


def readme():
    with open('README.md') as f:
        return f.read()


setup(name='UpsilonStabilityIndicator',
      version='0.1',
      description='Uses Arima modelling to assess the stability of a time series',
      long_description=readme(),
      classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.7',
        'Topic :: Time series :: ARIMA',
      ],
      keywords='',
      url='',
      author='Marie Rodal',
      author_email='marie.rodal@uantwerpen.be',
      license='MIT',
      packages=['UpsilonStabilityIndicator'],
      install_requires=[
          'rpy2',
          'time',
          'multiprocessing',
          'tqdm',
          'matplotlib', 
          'sdeint',
          'numpy',
      ],
      )




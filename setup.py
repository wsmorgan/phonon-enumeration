#!/usr/bin/env python
try:
    from setuptools import setup
    args = {}
except ImportError:
    from distutils.core import setup
    print("""\
*** WARNING: setuptools is not found.  Using distutils...
""")

from setuptools import setup
try:
    from pypandoc import convert
    def read_md(f):
        "Read md"
        return convert(f, 'rst')
except ImportError:
    print("warning: pypandoc module not found, could not convert Markdown to RST")
    def read_md(f):
        "Read md"
        return open(f, 'r').read()

setup(name='phenum',
      version='2.0.2',
      description='Enumeration of symmetrically unique derivative superstructures of crystals.',
      long_description=read_md('README.md'),
      author='Wiley S Morgan',
      author_email='wiley.s.morgan@gmail.com',
      url='https://github.com/wsmorgan/phonon-enumeration',
      license='MIT',
      install_requires=[
          "numpy",
          "termcolor",
          "pyparsing",
          "matplotlib",
      ],
      packages=['phenum'],
      scripts=['phenum/enumeration.py','phenum/makeStr.py'],
      classifiers=[
          'Development Status :: 5 - Production/Stable',
          'Intended Audience :: Science/Research',
          'Natural Language :: English',
          'License :: OSI Approved :: MIT License',          
          'Operating System :: MacOS',
          'Programming Language :: Python',
          'Programming Language :: Python :: 2',
          'Programming Language :: Python :: 2.7',
          'Programming Language :: Python :: 3',
          'Programming Language :: Python :: 3.4',
          'Topic :: Scientific/Engineering :: Physics',
      ],
     )

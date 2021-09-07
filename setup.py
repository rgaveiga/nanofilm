import os.path
from setuptools import setup

def read(fname):
    return open(os.path.join(os.path.dirname(__file__),fname)).read()

setup(name='nanofilm',
      version='0.1.5',
      description="Manipulation of molecules adsorbed onto a substrate.",
      long_description_content_type='text/markdown',
      long_description=read('Readme.MD'),
      author='Roberto Gomes de Aguiar Veiga',
      url="https://github.com/rgaveiga/nanofilm",
      packages=['nanofilm'])



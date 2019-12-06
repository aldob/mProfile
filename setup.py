from setuptools import setup, find_packages

setup(name='mprofile',
      version='1.0',
      description='nucleotid-resolution mutation calling',
      url='https://github.com/aldob/mProfile',
      packages=find_packages(),
      author='Aldo Bader',
      entry_points={
        'console_scripts': [
            'mprofile = mprofile.cli:main',
        ],
      },
      zip_safe=False)
 
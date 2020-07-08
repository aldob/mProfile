from setuptools import setup, find_packages

setup(name='mProfile',
      version='1.1',
      description='nucleotid-resolution mutation calling',
      url='https://github.com/aldob/mProfile',
      packages=find_packages(),
      author='Aldo Bader',
      entry_points={
        'console_scripts': [
            'callMUT = mProfile.callMUT:main',
            'TransloCapture = mProfile.TransloCapture:main'            
        ],
      },
      zip_safe=False)
 
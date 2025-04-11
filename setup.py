from setuptools import setup, find_packages

setup(name='mProfile-mut',
      version='2.1',
      description='nucleotide-resolution mutation calling',
      url='https://github.com/aldob/mProfile',
      download_url = 'https://github.com/aldob/mProfile/archive/v1.4.tar.gz',
      packages=find_packages(),
      author='Aldo Saleh Bader',
      install_requires=["pysam"],
      entry_points={
        'console_scripts': [
            'callMUT = mProfile.callMUT:main',
            'TransloCapture = mProfile.TransloCapture:main',
            'StructureMap = mProfile.StructureMap:main'
        ],
      },
      zip_safe=False)
 

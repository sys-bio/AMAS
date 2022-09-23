from setuptools import setup, find_packages

INSTALL_REQUIRES = [
    "compress-pickle",
    "coverage==5.5",
    "coveralls",
    "editdistance",
    "nose",
    "numpy",
    "pandas",
    "python-libsbml",
    "pyyaml",
    "requests",
    "scikit-learn"
    ]

def doSetup(install_requires):
  setup(
      name='AMASsb',
      version='0.0.3',
      author='Woosub Shin',
      author_email='woosubs@umich.edu',
      packages=find_packages(exclude=['tests', 'notebooks']),
      url='https://github.com/woosubs/AMAS',
      description='AMAS (Automatic Model Annotation System)',
      long_description=open('README.md').read(),
      long_description_content_type='text/markdown',
      package_dir={'AMAS': 'AMAS'},
      install_requires=install_requires,
      include_package_data=True,
      classifiers=[
          'Development Status :: 3 - Alpha',
          'Intended Audience :: Developers',      # Define that your audience are developers
          'Topic :: Scientific/Engineering :: Bio-Informatics',
          'License :: OSI Approved :: MIT License',   # Again, pick a license
          'Programming Language :: Python :: 3.8',
          'Programming Language :: Python :: 3.9',
        ],
      )

if __name__ == '__main__':
  doSetup(INSTALL_REQUIRES)

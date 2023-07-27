from setuptools import setup, find_packages

INSTALL_REQUIRES = [
    "compress-pickle",
    "editdistance",
    "nose",
    "numpy",
    "pandas",
    "python-libsbml",
    "pyyaml",
    ]

def doSetup(install_requires):
  setup(
      name='AMAS-sb',
      version='1.0.1',
      author='Woosub Shin',
      author_email='woosubs@umich.edu',
      packages=find_packages(exclude=['tests', 'notebooks']),
      url='https://github.com/sys-bio/AMAS',
      description='AMAS (Automatic Model Annotation System)',
      long_description=open('README.md').read(),
      long_description_content_type='text/markdown',
      package_dir={'AMAS': 'AMAS'},
      install_requires=install_requires,
      include_package_data=True,
      scripts=['AMAS/recommend_species.py',
               'AMAS/recommend_reactions.py',
               'AMAS/recommend_annotation.py',
               'AMAS/update_annotation.py',
               'AMAS/recommend_species',
               'AMAS/recommend_reactions',
               'AMAS/recommend_annotation',
               'AMAS/update_annotation',
               'AMAS/recommend_species.bat',
               'AMAS/recommend_reactions.bat',
               'AMAS/recommend_annotation.bat',
               'AMAS/update_annotation.bat',
               ],
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

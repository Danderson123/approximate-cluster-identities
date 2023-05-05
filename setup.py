from setuptools import setup, find_packages

setup(
    name='minimiser-distance',
    version='0.1',
    description='A package to calculate and visualise approximate jaccard distances for short nucleotide sequences using minimisers.',
    author='Daniel Anderson',
    author_email='danp.anderson@outlook.com',
    packages=find_packages(),
    install_requires=[
        'biopython',
        'pandas',
        'networkx',
        'numpy',
        'tqdm',
        'joblib'
    ],
    entry_points={
        'console_scripts': [
            'minimiser-distance=minimiser_distance:main',
        ],
    },
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.6',
    ],
)
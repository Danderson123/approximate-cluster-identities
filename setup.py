from setuptools import setup, find_packages

setup(
    name='approximate_cluster_identities',
    version='0.1',
    long_description_content_type = "text/markdown",
    description='A package to calculate and visualise approximate cluster identities for a large number of short nucleotide sequences using minimisers.',
    author='Daniel Anderson',
    author_email='danp.anderson@outlook.com',
    packages=find_packages(),
    install_requires=[
        'biopython',
        'pandas',
        'matplotlib',
        'networkx',
        'numpy',
        'seaborn',
        'tqdm',
        'joblib'
    ],
    entry_points={
        'console_scripts': [
            'aci=approximate_cluster_identities.__main__:main',
        ],
    },
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bioinformatics',
        'License :: Apache',
        'Programming Language :: Python :: 3.6',
    ]
)

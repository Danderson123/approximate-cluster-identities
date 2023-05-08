from setuptools import setup, find_packages

# Read the contents of your README file
with open('README.md', 'r', encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='approximate_cluster_identities',
    version='0.1.3',
    description='A package to calculate and visualise approximate cluster identities for a large number of short nucleotide sequences using minimisers.',
    long_description=long_description,
    long_description_content_type='text/markdown',
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
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.6',
    ]
)

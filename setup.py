from setuptools import setup


setup(
    name='d3pendr',
    version='0.1',
    description=(
        'Differential 3 Prime End analysis of Nanopore Direct RNAseq'
    ),
    author='Matthew Parker',
    entry_points={
        'console_scripts': [
            'd3pendr = d3pendr.main:d3pendr',
            'ndapars = d3pendr.main:ndapars',
        ]
    },
    scripts=[
        'scripts/filter_internal_priming.py',
        'scripts/filter_nanopore_oversplitting.py'
    ],
    packages=[
        'd3pendr',
    ],
    install_requires=[
        'numpy',
        'scipy',
        'pandas',
        'click',
        'pysam',
        #'pyBigWig',
        'statsmodels',
        'joblib'
    ],
)

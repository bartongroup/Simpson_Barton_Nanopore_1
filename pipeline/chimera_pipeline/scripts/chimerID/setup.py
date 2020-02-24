from setuptools import setup


setup(
    name='chimerID',
    version='0.1',
    description=(
        'find chimeric reads in nanopore data'
    ),
    author='Matthew Parker',
    entry_points={
        'console_scripts': [
            'chimerID = chimerID.main:chimerID',
        ]
    },
    packages=[
        'chimerID',
    ],
    install_requires=[
        'numpy',
        'pandas',
        'click',
        'joblib',
        'pysam',
    ],
)
from setuptools import setup


setup(
    name='pAks',
    version='0.1',
    description=(
        'Use output from nanopolish to find loci with different polyA tail length distributions'
    ),
    author='Matthew Parker',
    entry_points={
        'console_scripts': [
            'pAks = pAks.diff_pA:test_differential_polya',
            'pAmed = pAks.diff_pA:median_polya_lengths',
            'pAlabel = pAks.label:label_bam'
        ]
    },
    packages=[
        'pAks',
    ],
    install_requires=[
        'numpy',
        'pandas',
        'scipy',
        'click',
        'pysam',
        'statsmodels',
        'chimerID'
    ],
)
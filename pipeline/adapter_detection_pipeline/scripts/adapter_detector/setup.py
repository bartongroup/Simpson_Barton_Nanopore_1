from setuptools import setup


setup(
    name='adapter_detector',
    version='0.1',
    description=(
        'detect 5 prime adapters'
    ),
    author='Matthew Parker',
    entry_points={
        'console_scripts': [
            'adapter_detector = adapter_detector.filter:cli',
            'train_adapter_detector = adapter_detector.train:cli'
        ]
    },
    packages=[
        'adapter_detector',
    ],
    include_package_data=True,
    install_requires=[
        'numpy',
        'keras',
        'click',
        'click_log',
        'joblib',
        'pysam',
        'ont_fast5_api'
    ],
)

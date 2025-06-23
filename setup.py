from setuptools import setup
from setuptools import find_packages
from src.shumi.version import __version__

setup(
    name="shumi",
    description="Simulations of high-throughput UMI-based single-genome sequencing data.",
    version=__version__,
    url="https://github.com/niaid/umisim",
    author="Pierce Radecki",
    author_email="pierce.radecki@nih.gov",
    packages=find_packages('src'),
    package_dir={'': 'src'},
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License"
    ],
    entry_points={
        'console_scripts': [
            # 'umisim = umisim.cli:main',
            'shumi = shumi.cli:main'
        ]
    },
    install_requires=['biopython',
                      'matplotlib',
                      'numpy',
                      'pysam',
                      'pyyaml',
                      'tqdm',
                      'scipy',
                      'setuptools', ]
)

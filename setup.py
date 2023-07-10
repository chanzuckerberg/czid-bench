"""Packaging for czid-bench module.

References:
https://packaging.python.org/guides/distributing-packages-using-setuptools/
https://github.com/pypa/sampleproject
"""

from os import path
from setuptools import setup, find_packages
import czid_bench

local = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(local, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='czid-bench',
    # Versions should comply with PEP 440:
    # https://www.python.org/dev/peps/pep-0440/
    version=czid_bench.__version__,
    description='Tools to create and scores MGS benchmarks',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/chanzuckerberg/czid-bench',
    author='Chan Zuckerberg Initiative',
    author_email='help@czid.org',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Metagenomic Researchers',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ],
    keywords='czid benchmark metrics',
    py_modules=['czid-bench'],
    python_requires='>=3.6, <4',
    install_requires=[
        'InSilicoSeq>=1.4.2',
        'ncbi-acc-download',
        'numpy',
        'scikit-learn',
        'smart_open',
        'pyyaml'
    ],
    entry_points={
        'console_scripts': [
            'czid-bench-compare=czid_bench.compare:main',
            'czid-bench-generate=czid_bench.generate:main',
            'czid-bench-score=czid_bench.score:main',
        ],
    },
    packages=find_packages(),
    project_urls={
        'Docs': 'https://github.com/chanzuckerberg/czid-bench',
        'Bug Reports': 'https://github.com/chanzuckerberg/czid-bench/issues',
        'Sign up for IDseq': 'https://czid.org/',
    },
)

"""Packaging for idseq-bench module.

References:
https://packaging.python.org/guides/distributing-packages-using-setuptools/
https://github.com/pypa/sampleproject
"""

from os import path
from setuptools import setup
import idseq_bench

local = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(local, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='idseq-bench',
    # Versions should comply with PEP 440:
    # https://www.python.org/dev/peps/pep-0440/
    version=idseq_bench.__version__,
    description='Tools to create and scores MGS benchmarks',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/chanzuckerberg/idseq-bench',
    author='Chan Zuckerberg Initiative',
    author_email='help@idseq.net',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Metagenomic Researchers',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ],
    keywords='idseq benchmark metrics',
    py_modules=['idseq-bench'],
    python_requires='>=3.6, <4',
    install_requires=[
        'InSilicoSeq',
        'ncbi-acc-download'
    ],
    entry_points={
        'console_scripts': [
            'idseq-bench-compare=idseq_bench.compare:main',
            'idseq-bench-generate=idseq_bench.generate:main',
            'idseq-bench-score=idseq_bench.score:main',
        ],
    },
    packages=['idseq_bench'],
    project_urls={
        'Docs': 'https://github.com/chanzuckerberg/idseq-bench',
        'Bug Reports': 'https://github.com/chanzuckerberg/idseq-bench/issues',
        'Sign up for IDseq': 'https://idseq.net/',
    },
)
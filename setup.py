from setuptools import setup, find_packages

setup(
    name='pbjelly',
    author='acenglishbcm',
    version="16.0beta1",
    description='Software for Long-Read Sequencing Data from PacBio',
    python_requires='>=3',
    packages=find_packages(),
    install_requires=[
        'networkx==1.7',
        'pyparsing'
    ],
    scripts=[
        'pbsuite/jelly/Assembly.py',
        'pbsuite/jelly/Collection.py',
        'pbsuite/jelly/Extraction.py',
        'pbsuite/jelly/Jelly.py',
        'pbsuite/jelly/m4pie.py',
        'pbsuite/jelly/Setup.py',
        'pbsuite/jelly/Support.py',
        'pbsuite/utils/fakeQuals.py',
        'pbsuite/utils/summarizeAssembly.py',
        'pbsuite/utils/readSummary.py'
    ]
)

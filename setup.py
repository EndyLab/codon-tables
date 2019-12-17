import setuptools

with open('README.md', 'r') as handle:
    long_description = handle.read()

setuptools.setup(
    name='CodonTables-jecalles',
    version='0.0.1',
    author='Jonathan Calles',
    author_email='jecalles@stanford.edu',
    description='Tools for engineering genetic codes.',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/EndyLab/CodonTables'
)

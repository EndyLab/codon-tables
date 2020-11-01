import setuptools

with open('README.md', 'r') as handle:
    long_description = handle.read()

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

setuptools.setup(
    name='codontables',
    version='0.0.1',
    author='Jonathan Calles',
    author_email='jecalles@stanford.edu',
    description='Tools for engineering genetic codes.',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/EndyLab/codon-tables/tree/package',
    packages=setuptools.find_packages(),
    classifiers=[
        'Programming Language :: Python :: 3',
        'Licence :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>3.6',
    install_requires=requirements,
    package_data={'codontables':['res/*']}
)

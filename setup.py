from setuptools import setup

with open("README", 'r') as f:
    long_description = f.read()

with open("done/__init__.py", 'r') as f:
    for l in f:
        if l.startswith('__version__'):
            version = l.split('=')[1].strip().strip('"')
setup(
   name='SchistoAmpliPan',
   version=version,
   description='Design amplicon panels for targeted sequencing',
   author='Ritah Nabunje',
   author_email='ritahnabunje@gmail.com',
   long_description=long_description,
   packages=['Schisto AmpliPan'],
   install_requires=['vcf', 'csv', 'intervaltree', 'pandas', 'itertools', 'itertools', 'subprocess', 'Bio' , 'math', 'datetime', 'glob', 'os'],
)

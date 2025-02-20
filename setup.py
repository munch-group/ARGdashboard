import setuptools, os, sys

# read the contents of your README file
from os import path
this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

# SEARCH/REPLACE package_name WITH THE ACTUAL PACKAGE NAME BELOW:

setuptools.setup(
    name='popgen-dashboards',
    version="1.1.5",
    author="Kasper Munch",
    author_email="kaspermunch@birc.au.dk",
    description="Short description",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/kaspermunch/argdashboard',
    packages=setuptools.find_packages(),
    package_data={"argdashboard": ["assets/*.css"]},
    python_requires='>=3.6',
    entry_points = { 'console_scripts': [ 'arg-dashboard = argdashboard.arg_dashboard:run', ] },    
    install_requires=[
        'jupyterlab',
        'dash',
        'dash-bootstrap-components',
        'pandas',
        'networkx'
    ])


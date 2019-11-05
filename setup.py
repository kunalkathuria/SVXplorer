from setuptools import find_packages, setup
import re


def long_description():
    with open("README.md", "r") as f:
        d = f.read()
    return d


def version():
    return re.search('^__version__\s*=\s*"(.*)"', open('src/svxplorer/_version.py').read(), re.M).group(1)


setup(
    name='svxplorer',
    url='https://github.com/kunalkathuria/SVXplorer',
    author='Kunal Kathuria',
    # author_email='',  # TODO
    version=version(),
    # description='', #TODO
    long_description=long_description(),
    package_dir={"": "src"},
    packages=find_packages("src"),
    # install_requires=[], # TODO: move requirements to setup.py
    entry_points={
        'console_scripts': [
            'SVXplorer = svxplorer.SVXplorer:main',
        ],
    },
    classifiers=[
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',

    ]
)

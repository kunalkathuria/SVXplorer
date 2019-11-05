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
    author_email='kunal.kathuria@gmail.com',  # TODO check if correct
    version=version(),
    # description='', #TODO
    long_description=long_description(),
    package_dir={"": "src"},
    packages=find_packages("src"),
    install_requires=[
        "bitarray",
        "interlap",
        "logging",
        "networkx",
        "numpy",
        "pandas",
        "pybedtools",
        "pysam",
        "scikit-learn"
    ],
    entry_points={
        'console_scripts': [
            'SVXplorer = svxplorer.SVXplorer:main',
            'writeDiscordantFragments = svxplorer.writeDiscordantFragments:main',
            'writeBEDs = svxplorer.writeBEDs:main',
            # 'VCFtoBEDPE = svxplorer.vcftoBedpe:main', #TODO: uncomment when vcftoBedpe.py is updated
            # 'VCFtoBED = svxplorer.vcftoBed:main', #TODO: uncomment when vcftoBed.py is updated
            'uniqueSuppFilter = svxplorer.uniqueSuppFilter:main',  # TODO: suggested change of alias for users
            'preserveSmallClusters = svxplorer.preserveSmallClusters:main',
            'pickBestCluster = svxplorer.pickBestCluster:main',  # TODO: test if works
            # 'markLikelyConflicts = svxplorer.markLikelyConflicts:main'  #TODO: uncomment when markLikelyConflicts.py is updated
            'markDuplicateClusterRegions = svxplorer.markDuplicateClusterRegions:main',
            'formPEClusters = svxplorer.formPEClusters:main',
            'covPUFilter = svxplorer.covPUFilter:main',
            'consolidatePEClusters = svxplorer.consolidatePEClusters:main',
            # 'calcConflict = svxplorer.calcConflict:main', #TODO: uncomment when calcConflict.py is updated
            # 'addSplitReads = svxplorer.addSplitReads:main' #TODO: uncomment when addSplitReads.py is updated
        ],
    },
    classifiers=[
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',

    ]
)

from setuptools import setup, find_packages

setup(
    name='DICE',
    version='1.1.0',
    author='Samson Weiner',
    author_email='samson.weiner@uconn.edu',
    description='DICE (short for Distance-based Inference of Copy-number Evolution) is a collection of fast and accurate methods for reconstructing cell lineage trees from single-cell copy number aberration data.',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/samsonweiner/DICE',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent'
    ],
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'dice=core.dice:main',
        ],
    }
)
"""
Copyright (c) 2015 bigWigFeaturize
"""
from setuptools import setup, Extension

bigWigFeaturize_C_module = Extension(
    'bigWigFeaturize.bigWigFeaturize',
    define_macros = [('MAJOR_VERSION', '0'),
                     ('MINOR_VERSION', '1')],
    include_dirs = [
        '/usr/local/include/',
        '/users/nasa/nandi/include/python2.7/',
        '~/nandi/lib/python2.7/site-packages/numpy/core/include/',
        '/mnt/lab_data/kundaje/projects/encodeenhancerpredict/igr/kentUtils/src/inc'
    ],
    libraries = [
        'm', 
        'z', 
        'ssl', 
        '/mnt/lab_data/kundaje/projects/encodeenhancerpredict/igr/kentutils/src/lib/local/jkweb.a',
        '/mnt/lab_data/kundaje/projects/encodeenhancerpredict/igr/kentutils/samtabix/lib*.a'
    ],
    library_dirs = ['/usr/local/lib'],
    sources = ['src/bigWigFeaturize.c']
)

config = {
    'include_package_data': True,
    'ext_modules': [bigWigFeaturize_C_module,], 
    'description': 'bigWigFeaturize',
    'author': 'Nasa Sinnot-Armstrong',
    'url': 'NA',
    'download_url': 'https://github.com/nboley/bigWigFeaturize/',
    'author_email': 'nasa@stanford.edu',
    'version': '0.1',
    'packages': ['bigWigFeaturize', ],
    'setup_requires': [],
    'install_requires': ['pybedtools',],
    'scripts': [],
    'name': 'bigWigFeaturize'
}

if __name__== '__main__':
    setup(**config)

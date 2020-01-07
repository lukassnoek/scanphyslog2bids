import os
from setuptools import setup, find_packages
PACKAGES = find_packages()

# Get version and release info, which is all stored in shablona/version.py
ver_file = os.path.join('scanphyslog2bids', 'version.py')

with open(ver_file) as f:
    exec(f.read())

# Long description will go up on the pypi page
with open('README.md') as f:
    LONG_DESCRIPTION = f.read()

with open('requirements.txt') as f:
    REQUIRES = f.readlines()

opts = dict(name=NAME,
            maintainer=MAINTAINER,
            maintainer_email=MAINTAINER_EMAIL,
            description=DESCRIPTION,
            long_description_content_type='text/markdown',
            long_description=LONG_DESCRIPTION,
            url=URL,
            download_url=DOWNLOAD_URL,
            license=LICENSE,
            classifiers=CLASSIFIERS,
            author=AUTHOR,
            author_email=AUTHOR_EMAIL,
            platforms=PLATFORMS,
            version=VERSION,
            packages=PACKAGES,
            package_data=PACKAGE_DATA,
            install_requires=REQUIRES,
            requires=REQUIRES,
            entry_points={
                'console_scripts': [
                    'scanphyslog2bids = scanphyslog2bids.core:cmd_interface',
                    'add_header_to_physio_tsv = scanphyslog2bids.utils:add_header_to_physio_tsv'
                    ]
                }
            )

if __name__ == '__main__':
    setup(**opts)

# -*- coding: utf-8 -*-
#!/usr/bin/env python
#

from setuptools import setup, find_packages

import pathlib

try:
    from wheel.bdist_wheel import bdist_wheel as _bdist_wheel
    class bdist_wheel(_bdist_wheel):
        def finalize_options(self):
            _bdist_wheel.finalize_options(self)
            self.root_is_pure = False
except ImportError:
    bdist_wheel = None

here = pathlib.Path(__file__).parent.resolve()

long_description = (here / "README.md").read_text(encoding="utf-8")

setup(
    name="sitesyncro",
    version="1.0",
    description="SiteSyncro - Site-specific chronological modeling and synchronization",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/demjanp/SiteSyncro",
    author="Peter Demjan",
    author_email="peter.demjan@gmail.com",
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Programming Language :: Python :: 3.10",
    ],
    keywords="archaeology, radiocarbon, chronology, stratigraphy",
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    include_package_data=True,
	package_data={"":[
	]},
    python_requires=">=3.10, <3.11",
	install_requires=[
        'numpy>=1.26.4, <2',
        'scipy>=1.13.0, <2',
        'matplotlib>=3.8.4, <4',
        'scikit-learn>=1.4.2, <2',
        'tqdm>=4.66.0, <5',
        'requests>=2.31.0, <3',
        'networkx>=3.3.0, <4',
	],
	cmdclass={'bdist_wheel': bdist_wheel},
)

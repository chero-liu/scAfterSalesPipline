# setup.py

from setuptools import setup, find_packages

from scAfterSalesPipline.__init__ import __VERSION__


setup(
    name="scAfterSalesPipline",
    version=__VERSION__,
    description="scAfterSalesPipline for oe single cell RNA After Sales Pipline",
    author="liuchenglong",
    author_email="chenglong.liu@oebiotech.com",
    install_requires=[
        "pyyaml",
        "argparse",
        "numpy",
        "pandas",
        "matplotlib",
        "ruamel.yaml",
    ],
    packages=find_packages(),
    entry_points={
        "console_scripts": [
            "scAfterSalesPipline = scAfterSalesPipline.scAfterSalesPipline:main",
        ],
    },
    include_package_data=True,
)

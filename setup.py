from setuptools import setup, find_packages

setup(
    name="Geo2Gmsh",
    version="1.0.0",
    packages=find_packages(),
    install_requires=[
        "numpy",
        "gmsh"
    ],
)
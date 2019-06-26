import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="GRLP",
    version="1.0.0",
    author="Andrew D. Wickert",
    author_email="awickert@umn.edu",
    description="Evolves gravel-bed river long profiles",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/awickert/GRLP",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 2",
        "License :: OSI Approved :: GNU GPL v3.0",
        "Operating System :: OS Independent",
    ],
)

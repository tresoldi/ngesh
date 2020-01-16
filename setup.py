import pathlib
from setuptools import setup

# The directory containing this file
LOCAL_PATH = pathlib.Path(__file__).parent

# The text of the README file
README_FILE = (LOCAL_PATH / "README.md").read_text()

# This call to setup() does all the work
setup(
    name="ngesh",
    version="0.3.1",
    description="Simulate random phylogenetic trees",
    long_description=README_FILE,
    long_description_content_type="text/markdown",
    url="https://github.com/tresoldi/ngesh",
    author="Tiago Tresoldi",
    author_email="tresoldi@shh.mpg.de",
    license="MIT",
    classifiers=[
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
        "Topic :: Software Development :: Libraries",
    ],
    packages=["ngesh"],
    keywords=["random phylogenetic tree", "phylogenetics"],
    include_package_data=True,
    install_requires=["abzu", "ete3", "numpy", "six"],
    entry_points={"console_scripts": ["ngesh=ngesh.__main__:main"]},
    test_suite="tests",
    tests_require=[],
)

from setuptools import setup, find_packages

setup(
    name="PySurv",
    version="0.1.1",
    author="Pushkar Mutha",
    author_email="pushkar.mutha@outlook.com",
    description="A Python tool for generating and plotting Kaplan-Meier survival curves.",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/pushkarmutha/PySurv",  # Update with your GitHub repo
    packages=find_packages(),
    install_requires=[
        "numpy",
        "pandas",
        "lifelines",
        "matplotlib"
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ],
    python_requires=">=3.6",
)
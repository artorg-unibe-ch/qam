import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="quantitative-ablation-margin",
    version="0.1.0",
    author="Raluca M. Sandu",
    author_email="raluca-sandu@rwth-aachen.de",
    description="Calculation and visualization of ablation margins",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/artorg-unibe-ch/qam",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GPL :: 3 ",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)

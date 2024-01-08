import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(name='stepRNA',
    version='1.0.0',
    author='Ben Murcott',
    author_email='bmm41@bath.ac.uk',
    description='A package to look at inserting DNA sequences into plasmids for cloning or PCR for SWBio DTP Data Science and Machine Learning Module.',
    long_description = long_description,
    long_description_content_type = "text/markdown",
    url='https://github.com/bmm514/PlasmidInsertChecker.git',
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: Unix"
    ],
    packages=["plasmidin"],
    python_requires=">=3.8",
    install_requires=[
        "biopython==1.81",
        "numpy==1.26.2",
        "pandas==2.1.3",
        "reportlab==4.0.8"
        ]
)
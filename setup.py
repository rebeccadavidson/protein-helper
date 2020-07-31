import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="protein-helper",
    version="0.0.1",
    author="Rebecca Davidson",
    author_email="becca.davidson@gmail.com",
    description="A package to help with protein function predictions tasks.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/rebeccadavidson/protein-helper",
    package_dir={"": "src"},
    packages=setuptools.find_packages(where="src"),
    python_requires='>=3.7',
    install_requires=[
        'biopython',
        'click',
        'matplotlib',
        'networkx',
    ],
    entry_points={
        'console_scripts': [
            'protein-helper=protein_helper.scripts.protein_helper:cli',
        ]
    }

)

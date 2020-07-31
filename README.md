# protein-helper
A library to help with protein function prediction tasks.

## Getting started

Follow these instructions to get this project running on your local machine for developement, testing, and running analyses locally.

###Prerequisites

In order to install this project you will need to install bioinformatics software to your system.

#### Install dependencies on a Mac with Homebrew
    brew install diamond
    brew install pipenv

Install other dependencies from source: cd-hit

#### Create the environment
    pipenv install --dev

#### Activate the virtual environment
    pipenv shell

## Running commands
One may run help to get a comprehensive list of the commands available:

    protein-helper -h

This is an an example of a command that creates a cytoscape network representing protein similarity 
among a set of proteins:

    protein-helper network \
    --input-fasta curated_seqs.fa \
    --cytoscape-cyjs curated_seqs.cyjs \
    --min-percent-identity 65 \
    --temp-dir cdhit_results \
    --threads 4

This file can then be imported into Cytoscape for network analysis.
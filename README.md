# DrugBank Knowledge Graph Pipeline

This repository contains a Python script for transforming data from the DrugBank knowledge base into a knowledge graph. The script extracts data, converts it into RDF triples using the Schema.org vocabulary, and updates the knowledge graph in an AllegroGraph database.

## What the script does:
- Automatically extracts data from DrugBank XML file (acquired trough an academic license - limited data).
- Maps DrugBank data to RDF triples using Schema.org.
- Uploads RDF triples to AllegroGraph for storage and visualization.

## Requirements
- Python 3.12
- rdflib
- lxml
- agraph-python
- python-dotenv (optional)
- DrugBank XML file

## Usage
1. Clone the repo: `git clone https://github.com/Ackolakitus/DrugBank-knowledge-graph.git`

2. Install dependencies: `pip install rdflib lxml agraph-python python-dotenv`

3. Download the latest release of [DrugBank](https://go.drugbank.com/releases) XML (all drugs) file.

4. Run the script `python extract_and_create.py`

## Note

Check the GitHub Actions workflow for automation of the whole process.

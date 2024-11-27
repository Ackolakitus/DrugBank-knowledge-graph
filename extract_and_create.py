import json
import re
import sys
import time
from datetime import date

from dotenv import load_dotenv
from rdflib import Graph, URIRef, Literal, BNode
from rdflib.namespace import SDO, RDF
import lxml.etree as ET

from franz.openrdf.connect import ag_connect

jena_update_url = "http://localhost:3030/drugbank/update"
jena_sparql_url = "http://localhost:3030/drugbank/sparql"

ns = {'drugbank': 'http://www.drugbank.ca'}
classyfire_base_uri = 'http://classyfire.wishartlab.com/tax_nodes'

drug_ext_ids_uris = {
    "UniProtKB": "https://www.uniprot.org/uniprotkb/IDENTIFIER/entry",
    "Wikipedia": "https://en.wikipedia.org/wiki/IDENTIFIER",
    "ChEBI": "https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:IDENTIFIER",
    "ChEMBL": "https://www.ebi.ac.uk/chembl/web_components/explore/compound/IDENTIFIER",
    "PubChem Compound": "https://pubchem.ncbi.nlm.nih.gov/compound/IDENTIFIER",
    "PubChem Substance": "https://pubchem.ncbi.nlm.nih.gov/substance/IDENTIFIER",
    "KEGG Compound": "https://www.kegg.jp/entry/IDENTIFIER",
    "KEGG Drug": "https://www.kegg.jp/entry/IDENTIFIER",
    "ChemSpider": "https://www.chemspider.com/Chemical-Structure.IDENTIFIER.html",
    "BindingDB": "https://www.bindingdb.org/rwd/bind/chemsearch/marvin/MolStructure.jsp?monomerid=IDENTIFIER",
    "UniProt Accession": "https://www.uniprot.org/uniprotkb/IDENTIFIER/entry",
    "GenBank Protein Database": "https://www.ncbi.nlm.nih.gov/protein/IDENTIFIER",
    "HUGO Gene Nomenclature Committee (HGNC)": "https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/IDENTIFIER",
    "GenAtlas": "http://genatlas.medecine.univ-paris5.fr/fiche.php?symbol=IDENTIFIER",
    "PharmGKB": "https://www.pharmgkb.org/chemical/IDENTIFIER",
    "PDB": "https://www.rcsb.org/ligand/IDENTIFIER",
    "IUPHAR": "https://www.guidetopharmacology.org/GRAC/LigandDisplayForward?ligandId=IDENTIFIER",
    "Guide to Pharmacology": "https://www.guidetopharmacology.org/GRAC/ObjectDisplayForward?objectId=IDENTIFIER",
    "ZINC": "https://zinc15.docking.org/substances/IDENTIFIER/",
    "RxCUI": "https://mor.nlm.nih.gov/RxNav/search?searchBy=RXCUI&searchTerm=IDENTIFIER",
    "Drugs Product Database (DPD)": None,
    "National Drug Code Directory": None,
    "GenBank Gene Database": None,
    "Therapeutic Targets Database": None,
    "polypeptides": "https://go.drugbank.com/polypeptides/IDENTIFIER",
    "drugs": "https://go.drugbank.com/drugs/IDENTIFIER",
    "bioentities": "https://go.drugbank.com/bio_entities/IDENTIFIER",
    "salts": "https://go.drugbank.com/salts/IDENTIFIER",
    "metabolites": "https://go.drugbank.com/metabolites/",
    "cas": "https://commonchemistry.cas.org/detail?cas_rn=IDENTIFIER",
    "unii": "https://precision.fda.gov/uniisearch/srs/unii/IDENTIFIER",
    "rs": "https://www.ncbi.nlm.nih.gov/snp/IDENTIFIER"
}

group_nodes = {
    'approved': "?approved=1",
    'nutraceutical': "?approved=0&nutraceutical=1",
    'illicit': "?approved=0&illicit=1",
    'withdrawn': "?approved=0&withdrawn=1",
    'investigational': "?approved=0&investigational=1",
    'experimental': "?approved=0&experimental=1",
}

group_descriptions = {
    'approved': 'A drug that has been approved in at least one jurisdiction, at some point in time. This does not mean the drug is currently approved or available, just that it has been approved and marketed at some point, somewhere. Different jurisdictions also have a different concept for "approval". For example drugs that are available over-the-counter in the U.S. may not be technically approved, whereas over-the-counter drugs in Canada are considered approved.',
    'nutraceutical': 'In DrugBank, nutraceuticals are drugs which are regulated and processed at a pharmaceutical grade and have a demonstrable nutritional effect. Another major characteristic of nutraceuticals is that they do not always have a patent protection even when used as therapeutic agents.',
    'illicit': 'In DrugBank this status is given to drugs whereby their use has been prohibited due to either the stimulation or inhibition of the central nervous system or due to the production of hallucinogenic effects. Illicit drugs are often scheduled drugs, whereby they are approved but are limited in their distribution or have additional constraints on how they can be prescribed.',
    'withdrawn': 'In DrugBank this state is given to the drugs that have been discontinued.  Although the reason as to why such drugs are withdrawn can range from patient safety and toxicity issues to limited commercial viability, the formal rationale behind such decisions lies typically with the regional or national public health administration that issued the withdrawal.\nA drug maintains its drug approval status even once withdrawn, as approval indicates that it has been approved at some point in time.',
    'investigational': 'In DrugBank an investigational drug refers to the drug development status in which the drug is being researched for a determinate condition and has reached clinical trials. A drug can have multiple statuses filled if for example, it has been approved for a determinate  indication but it is currently on clinical trials for a different indication.',
    'experimental': 'A compound that has been shown experimentally to bind specific proteins in mammals, bacteria, viruses, fungi, or parasites. This includes compounds that are Pre-Investigational New Drug Applications (Pre-IND, or Discovery Phase compounds).',
}

biotech_uri = "https://go.drugbank.com/biotech_drugs/"
small_molecule_uri = "https://go.drugbank.com/drugs"

biotech_term_set = URIRef(biotech_uri)
small_molecule_term_set = URIRef(small_molecule_uri)
pathways_term_set = URIRef("https://smpdb.ca/view?subject=")
snp_effects = URIRef("https://www.ncbi.nlm.nih.gov/snp/")


def save_locally(graph, filename, data_format='turtle'):
    with open(filename, "wb") as f:
        graph.serialize(f, format=data_format)

    print("Saved to " + filename)


def get_element_text(element, path):
    xml = element.find(path, ns)

    if xml is None:
        return None

    if xml.text in (None, '') or xml.text.lower() == 'none':
        return None

    return xml.text.strip()


def extract_text_with_stars(text):
    pattern = r'(\*\*.*?\*\*)(.*?)(?=\*\*|$)'
    matches = re.findall(pattern, text, re.DOTALL)

    result = [header + "\n\n" + content.strip() for header, content in matches]

    if not matches:
        result = [text]

    return result


def map_name(drug_xml, graph, drug):
    name = get_element_text(drug_xml, 'drugbank:name')

    if not name:
        return

    graph.add((drug, SDO.name, Literal(name)))


def map_description(drug_xml, graph, drug):
    description = get_element_text(drug_xml, 'drugbank:description')

    if description is None:
        return

    description = description.replace("\r", "")
    graph.add((drug, SDO.description, Literal(description)))


def map_groups(drug_xml, graph, drug):
    drug_groups = [group.text for group in drug_xml.findall('drugbank:groups/drugbank:group', ns)]
    if not drug_groups:
        return

    drug_type = drug_xml.attrib.get('type')

    for gr in drug_groups:
        if drug_type is not None and gr != 'vet_approved':
            if drug_type == 'biotech':
                group_uri = URIRef(biotech_uri + group_nodes.get(gr))
                graph.add((drug, SDO.taxonomicRange, group_uri))
            elif drug_type == "small molecule":
                group_uri = URIRef(small_molecule_uri + group_nodes.get(gr))
                graph.add((drug, SDO.taxonomicRange, group_uri))


def map_label_details(drug_xml, graph, drug):
    label_details_url = get_element_text(drug_xml, 'drugbank:fda-label')

    if label_details_url is None:
        return

    if str(label_details_url).startswith("//s3"):
        label_details_url = "https:" + label_details_url

    label_uri = URIRef(label_details_url)
    graph.add((drug, SDO.labelDetails, label_uri))


def map_indications(drug_xml, graph, drug):
    indication_text = get_element_text(drug_xml, 'drugbank:indication')

    if indication_text is None:
        return

    term = BNode()
    graph.add((drug, SDO.potentialUse, term))
    graph.add((term, RDF.type, SDO.DefinedTerm))
    graph.add((term, SDO.name, Literal('Indication')))
    graph.add((term, SDO.description, Literal(indication_text)))


def map_toxicity(drug_xml, graph, drug):
    toxicity_text = get_element_text(drug_xml, 'drugbank:toxicity')
    if toxicity_text is None:
        return

    toxicities = extract_text_with_stars(toxicity_text)

    for t in toxicities:
        graph.add((drug, SDO.overdosage, Literal(t.replace('\r', ''))))


def map_monoisotopic_mass(drug_xml, graph, drug):
    monoisotopic_mass = get_element_text(drug_xml, 'drugbank:monoisotopic-mass')

    if not monoisotopic_mass:
        return

    graph.add((drug, SDO.monoisotopicMolecularWeight, Literal(monoisotopic_mass)))


def map_synonyms(drug_xml, graph, drug):
    synonyms = [synonym.text for synonym in drug_xml.findall('drugbank:synonyms/drugbank:synonym', ns)]

    if not synonyms:
        return

    for s in synonyms:
        graph.add((drug, SDO.alternateName, Literal(s)))


def map_routes(drug_xml, graph, drug):
    routes_list = [get_element_text(route, 'drugbank:route')
                   for route in
                   drug_xml.findall('drugbank:dosages/drugbank:dosage', ns)]

    if not routes_list:
        return

    routes = {
        item.strip()
        for list_item in routes_list if list_item
        for item in list_item.split(";") if item.strip()
    }

    for route in routes:
        graph.add((drug, SDO.administrationRoute, Literal(route)))


def map_strengths(drug_xml, graph, drug):
    strengths = {text for strength in drug_xml.findall('drugbank:dosages/drugbank:dosage', ns)
                 if (text := strength.findtext('drugbank:strength', None, ns).lower().strip())}

    if not strengths:
        return


    for strength in strengths:
        graph.add((drug, SDO.availableStrength, Literal(strength)))

        # value, unit = extract_value_unit(strength)
        # identifier = BNode()
        # graph.add((identifier, RDF.type, SDO.DrugStrength))
        # graph.add((drug, SDO.strengthValue, value))
        # graph.add((drug, SDO.strengthUnit, Literal(unit)))


def map_food_warnings(drug_xml, graph, drug):
    food_warning = [interaction.text for interaction in
                    drug_xml.findall('drugbank:food-interactions/drugbank:food-interaction', ns) if
                    "alcohol" not in interaction.text.lower()]

    if not food_warning:
        return

    for fw in food_warning:
        graph.add((drug, SDO.foodWarning, Literal(fw)))


def map_alcohol_warnings(drug_xml, graph, drug):
    alcohol_warning = [interaction.text for interaction in
                       drug_xml.findall('drugbank:food-interactions/drugbank:food-interaction', ns) if
                       "alcohol" in interaction.text.lower()]

    if not alcohol_warning:
        return

    for aw in alcohol_warning:
        graph.add((drug, SDO.alcoholWarning, Literal(aw)))


def map_drug_interactions(drug_xml, graph, drug):
    drug_interactions = [interaction.find('drugbank:drugbank-id', ns).text for interaction in
                         drug_xml.findall('drugbank:drug-interactions/drugbank:drug-interaction', ns)]

    if not drug_interactions:
        return

    for di in drug_interactions:
        graph.add((drug, SDO.interactingDrug, URIRef(f'https://go.drugbank.com/drugs/{di}')))


def map_pharmacodynamics(drug_xml, graph, drug):
    pharmacodynamics = get_element_text(drug_xml, 'drugbank:pharmacodynamics')

    if pharmacodynamics is None:
        return

    pharmacodynamics = pharmacodynamics.replace("*", "").replace("\r", "").strip()
    graph.add((drug, SDO.clinicalPharmacology, Literal("Pharmacodynamics:\n\n" + pharmacodynamics)))


def map_route_of_elimination(drug_xml, graph, drug):
    route_of_elimination = get_element_text(drug_xml, 'drugbank:route-of-elimination')

    if route_of_elimination is None:
        return

    route_of_elimination = route_of_elimination.replace("*", "").replace("\r", "").strip()
    graph.add((drug, SDO.clinicalPharmacology, Literal("Route of elimination:\n\n" + route_of_elimination)))


def map_absorption(drug_xml, graph, drug):
    absorption = get_element_text(drug_xml, 'drugbank:absorption')

    if absorption is None:
        return

    absorption = absorption.replace("*", "").replace("\r", "").strip()
    graph.add((drug, SDO.clinicalPharmacology, Literal("Absorption:\n\n" + absorption)))


def map_half_life(drug_xml, graph, drug):
    half_life = get_element_text(drug_xml, 'drugbank:half-life')
    if half_life is None:
        return

    half_life = half_life.replace("*", "").replace("\r", "").strip()
    graph.add((drug, SDO.clinicalPharmacology, Literal("Half-life:\n\n" + half_life)))


def map_volume_of_distribution(drug_xml, graph, drug):
    volume_of_distribution = get_element_text(drug_xml, 'drugbank:volume-of-distribution')
    if volume_of_distribution is None:
        return

    volume_of_distribution = volume_of_distribution.replace("*", "").replace("\r", "").strip()
    graph.add((drug, SDO.clinicalPharmacology, Literal("Volume of distribution:\n\n" + volume_of_distribution)))


def map_clearance(drug_xml, graph, drug):
    clearance = get_element_text(drug_xml, 'drugbank:clearance')
    if clearance is None:
        return

    clearance = clearance.replace("*", "").replace("\r", "").strip()
    graph.add((drug, SDO.clinicalPharmacology, Literal("Clearance:\n\n" + clearance)))


def map_protein_binding(drug_xml, graph, drug):
    protein_binding = get_element_text(drug_xml, 'drugbank:protein-binding')

    # if protein_binding is None:
    #     return
    #
    # protein_binding = protein_binding.strip()
    # protein_binding_node = BNode()
    #
    # graph.add((drug, SDO.hasMolecularFunction, protein_binding_node))
    # graph.add((protein_binding_node, RDF.type, SDO.PropertyValue))
    # graph.add((protein_binding_node, SDO.description, Literal(
    #     "A description of the drug’s affinity for plama proteins and the proportion of the drug that is bound to them when in circulation within the body.")))
    # graph.add((protein_binding_node, SDO.value, Literal(protein_binding)))


def map_mechanism_of_action(drug_xml, graph, drug):
    mechanism_of_action = get_element_text(drug_xml, 'drugbank:mechanism-of-action')

    if mechanism_of_action is None:
        return

    mechanism_of_action = mechanism_of_action.replace("*", "").replace("\r", "").strip()
    graph.add((drug, SDO.mechanismOfAction, Literal(mechanism_of_action)))


def map_over_the_counter(drug_xml, graph, drug):
    over_the_counter = False

    products = [product for product in drug_xml.findall('drugbank:products/drugbank:product', ns)]

    if not products:
        graph.add((drug, SDO.prescriptionStatus, SDO.PrescriptionOnly))
        return

    # NDC LINKS <https://ndclist.com/ndc/{}>
    for product in products:
        otc = product.find('drugbank:over-the-counter', ns)
        if otc is not None and otc.text == "true":
            over_the_counter = True
            break

    status = SDO.OTC if over_the_counter else SDO.PrescriptionOnly
    graph.add((drug, SDO.prescriptionStatus, Literal(status)))


def map_proprietary_names(drug_xml, graph, drug):
    proprietary_names = {product_name.strip()
                         for product in drug_xml.findall('drugbank:products/drugbank:product', ns)
                         if (product_name := get_element_text(product, 'drugbank:name'))}

    if not proprietary_names:
        return

    for name in proprietary_names:
        if name:
            graph.add((drug, SDO.proprietaryName, Literal(name)))


def map_external_identifiers(drug_xml, graph, drug):
    external_identifiers = [{
        'resource': get_element_text(identifier, 'drugbank:resource'),
        'identifier': get_element_text(identifier, 'drugbank:identifier')
    } for identifier in drug_xml.findall('drugbank:external-identifiers/drugbank:external-identifier', ns)
    ]
    if not external_identifiers:
        return

    drug_info = {
        'unii': get_element_text(drug_xml, 'drugbank:unii'),
        'cas-number': get_element_text(drug_xml, 'drugbank:cas-number')}

    if drug_info.get('unii') is not None:
        graph.add(
            (drug, SDO.sameAs, URIRef(drug_ext_ids_uris.get('unii').replace("IDENTIFIER:", drug_info.get('unii')))))

    if drug_info.get('cas-number') is not None:
        graph.add((drug, SDO.sameAs,
                   URIRef(drug_ext_ids_uris.get('cas').replace("IDENTIFIER:", drug_info.get('cas-number')))))

    for item in external_identifiers:

        ext_id = item.get('identifier')
        resource = item.get('resource')
        if ext_id is None or resource is None:
            return

        ext_id_uri = drug_ext_ids_uris.get(resource)

        if not ext_id_uri:
            continue

        uri = URIRef(ext_id_uri.replace("IDENTIFIER", ext_id))
        if resource != "Wikipedia":
            graph.add((drug, SDO.code, uri))
            graph.add((uri, RDF.type, SDO.MedicalCode))
            graph.add((uri, SDO.codeValue, Literal(ext_id)))
            graph.add((uri, SDO.codingSystem, Literal(resource)))
        else:
            same_as_uri = URIRef(ext_id_uri.replace("IDENTIFIER", ext_id))
            graph.add((drug, SDO.sameAs, same_as_uri))
            # graph.add((same_as_uri, RDF.type, SDO.URL))


def map_pdb_entries(drug_xml, graph, drug):
    pdb_entries = [entry.text for entry in drug_xml.findall('drugbank:pdb-entries/drugbank:pdb-entry', ns) if
                   entry is not None]

    if not pdb_entries:
        return

    for p in pdb_entries:
        entry_uri = URIRef(f'https://www.wwpdb.org/pdb?id={p}')
        graph.add((drug, SDO.mainEntityOfPage, entry_uri))
        # graph.add((entry_uri, RDF.type, SDO.URL))


def map_salts(drug_xml, graph, drug):
    salts = [{
        'drugbank-id': get_element_text(salt, 'drugbank:drugbank-id'),
        'name': get_element_text(salt, 'drugbank:name'),
        'unii': get_element_text(salt, 'drugbank:unii'),
        'cas-number': get_element_text(salt, 'drugbank:cas-number'),
        'inchikey': get_element_text(salt, 'drugbank:inchikey'),
        'average-mass': get_element_text(salt, 'drugbank:average-mass'),
        'monoisotopic-mass': get_element_text(salt, 'drugbank:monoisotopic-mass')}
        for salt in drug_xml.findall('drugbank:salts/drugbank:salt', ns)
    ]

    if not salts:
        return

    for s in salts:
        if s.get('drugbank-id') is not None:
            salt = URIRef(f'https://go.drugbank.com/salts/{s.get("drugbank-id")}')
            graph.add((drug, SDO.hasBioChemEntityPart, salt))

            graph.add((salt, RDF.type, SDO.MolecularEntity))
            graph.add((salt, SDO.identifier, Literal(s.get('drugbank-id'))))
            graph.add((salt, SDO.additionalType, Literal('Salt')))

            if s.get('name'):
                graph.add((salt, SDO.name, Literal(s.get('name'))))
            if s.get('inchikey'):
                graph.add((salt, SDO.inChIKey, Literal(s.get('inchikey'))))
            if s.get('average-mass'):
                graph.add((salt, SDO.molecularWeight, Literal(s.get('average-mass'))))
            if s.get('monoisotopic-mass'):
                graph.add((salt, SDO.monoisotopicMolecularWeight, Literal(s.get('monoisotopic-mass'))))

            if s.get('unii') is not None:
                graph.add((salt, SDO.sameAs,
                           URIRef(drug_ext_ids_uris.get('unii').replace("IDENTIFIER:", s.get('unii')))))

            if s.get('cas') is not None:
                graph.add((salt, SDO.sameAs,
                           URIRef(drug_ext_ids_uris.get('cas').replace("IDENTIFIER:", s.get('cas-number')))))

            graph.add((salt, SDO.isPartOfBioChemEntity, drug))
            graph.add((drug, SDO.hasBioChemEntityPart, salt))


def map_calculated_properties(drug_xml, graph, drug):
    calculated_properties = {
        get_element_text(prop, 'drugbank:kind'): get_element_text(prop, 'drugbank:value')
        for prop in drug_xml.findall('drugbank:calculated-properties/drugbank:property', ns)
    }
    experimental_properties = {
        get_element_text(prop, 'drugbank:kind'): get_element_text(prop, 'drugbank:value')
        for prop in drug_xml.findall('drugbank:experimental-properties/drugbank:property', ns)
    }

    calculated_properties.update(experimental_properties)

    if not calculated_properties:
        return

    for prop, value in calculated_properties.items():
        if prop is None or value is None:
            continue

        match prop:
            case 'InChI':
                graph.add((drug, SDO.inChI, Literal(calculated_properties.get('InChI'))))
            case'InChIKey':
                graph.add((drug, SDO.inChIKey, Literal(calculated_properties.get('InChIKey'))))
            case 'IUPAC Name':
                graph.add((drug, SDO.iupacName, Literal(calculated_properties.get('IUPAC Name'))))
            case'Molecular Formula':
                graph.add((drug, SDO.molecularFormula, Literal(calculated_properties.get('Molecular Formula'))))
            case 'Molecular Weight':
                graph.add((drug, SDO.molecularWeight, Literal(calculated_properties.get('Molecular Weight'))))
            case 'SMILES':
                graph.add((drug, SDO.smiles, Literal(calculated_properties.get('SMILES'))))
            case _:
                continue
                # additional_prop = BNode()
                # graph.add((drug, SDO.additionalProperty, additional_prop))
                # graph.add((additional_prop, SDO.name, Literal(prop)))
                # graph.add((additional_prop, SDO.value, Literal(value)))


def map_atc_codes(drug_xml, graph, drug):
    atc_codes = []

    for code in drug_xml.findall('drugbank:atc-codes/drugbank:atc-code', ns):
        levels = code.findall('drugbank:level', ns)
        n = len(levels)

        local_codes = []
        for i, level in enumerate(reversed(levels)):
            local_codes.append((i, level.attrib.get('code'), level.text))

        local_codes.append((n, code.attrib.get('code'), None))

        atc_codes.append(local_codes)

    if not atc_codes:
        return

    for code in atc_codes:
        for level in code:
            _, level_code, level_name = level
            graph.add((drug, SDO.code, URIRef(f'https://atcddd.fhi.no/atc_ddd_index/?code={level_code}')))
            graph.add((URIRef(f'https://atcddd.fhi.no/atc_ddd_index/?code={level_code}'), RDF.type, SDO.MedicalCode))

            graph.add(
                (URIRef(f'https://atcddd.fhi.no/atc_ddd_index/?code={level_code}'), SDO.codeValue, Literal(level_code)))
            graph.add((URIRef(f'https://atcddd.fhi.no/atc_ddd_index/?code={level_code}'), SDO.codingSystem,
                       Literal('Anatomical Therapeutic Chemical (ATC) classification system')))
            if level_name is not None:
                graph.add(
                    (URIRef(f'https://atcddd.fhi.no/atc_ddd_index/?code={level_code}'), SDO.name, Literal(level_name)))


def map_drug_pathways(drug_xml, graph, drug):
    drug_pathways = [{
        'smpdb-id': get_element_text(pathway, 'drugbank:smpdb-id'),
        'name': get_element_text(pathway, 'drugbank:name'),
        'category': get_element_text(pathway, 'drugbank:category'),
        'drugs': [get_element_text(d, 'drugbank:drugbank-id') for d in
                  pathway.findall('drugbank:drugs/drugbank:drug', ns)],
        'enzymes': [enzyme.text.strip() for enzyme in pathway.findall('drugbank:enzymes/drugbank:uniprot-id', ns)]
    } for pathway in drug_xml.findall('drugbank:pathways/drugbank:pathway', ns)]

    if not drug_pathways:
        return

    for pathway in drug_pathways:
        # category = pathway.get('category').replace("_","+").title()
        if pathway.get('smpdb-id') is not None:
            pathway_term = URIRef(f'https://smpdb.ca/view/{pathway.get("smpdb-id")}')
            graph.add((pathway_term, RDF.type, SDO.DefinedTerm))
            graph.add((pathway_term, SDO.inDefinedTermSet, pathways_term_set))
            graph.add((pathway_term, SDO.name, Literal(pathway.get('name'))))

            for d in pathway.get('drugs'):
                drug_uri = URIRef(f'https://go.drugbank.com/drugs/{d}')
                graph.add((drug_uri, RDF.type, SDO.MolecularEntity))
                graph.add((drug_uri, RDF.type, SDO.Drug))
                graph.add((drug_uri, SDO.isInvolvedInBiologicalProcess, pathway_term))

            for e in pathway.get('enzymes'):
                # enzyme_uri = URIRef(f'https://www.uniprot.org/uniprotkb/{e}/entry')
                enzyme_uri = URIRef(f'https://go.drugbank.com/polypeptides/{e}')
                graph.add((enzyme_uri, RDF.type, SDO.Protein))
                graph.add((enzyme_uri, SDO.additionalType, Literal("Enzyme")))
                # graph.add((enzyme_uri, SDO.sameAs, URIRef(f'https://www.uniprot.org/uniprotkb/{e}/entry')))
                graph.add((enzyme_uri, SDO.isInvolvedInBiologicalProcess, pathway_term))


def map_drug_targets_enzymes_carriers_transporters(drug_xml, graph, drug):
    # A drug target is a specific molecule, often a protein, in the body that is closely linked to a particular disease process and can be influenced by a drug to produce a desired therapeutic outcome.
    # Enzymes that are inhibited/induced or involved in metabolism
    # Carrier or transporter proteins involved in movement of the drug across biological membranes
    all_elements = (drug_xml.findall('drugbank:targets/drugbank:target', ns)
                    + drug_xml.findall('drugbank:enzymes/drugbank:enzyme', ns)
                    + drug_xml.findall('drugbank:carriers/drugbank:carrier', ns)
                    + drug_xml.findall('drugbank:transporters/drugbank:transporter', ns))

    drug_targets_enzymes_carriers_transporters = [{
        'id': get_element_text(target, 'drugbank:id'),
        'name': get_element_text(target, 'drugbank:name'),
        'actions': [action.text for action in target.findall('drugbank:actions/drugbank:action', ns)],
        'polypeptide-id': target.find('drugbank:polypeptide', ns).attrib.get('id') if target.find(
            'drugbank:polypeptide', ns) is not None else None,
    } for target in all_elements]

    for element in all_elements:
        if (p_id := element.find('drugbank:polypeptide', ns).attrib.get('id') if element.find(
                'drugbank:polypeptide', ns) is not None else None) is not None:
            map_general_references(element, graph, URIRef(
                drug_ext_ids_uris.get('polypeptides').replace("IDENTIFIER", p_id)), 'drugbank:references')

    if not drug_targets_enzymes_carriers_transporters:
        return

    for entry in drug_targets_enzymes_carriers_transporters:
        if entry.get('polypeptide-id') is not None:
            polypeptide = URIRef(
                drug_ext_ids_uris.get('polypeptides').replace("IDENTIFIER", entry.get('polypeptide-id')))

            graph.add((drug, SDO.bioChemInteraction, polypeptide))
            graph.add((polypeptide, SDO.bioChemInteraction, drug))

            if (e_id := entry.get('id')) is not None:
                # same_as_uri = URIRef(f'{drug_ext_ids_uris.get('bioentities').replace("IDENTIFIER", entry.get('id'))}')
                graph.add((polypeptide, SDO.identifier, Literal(e_id)))
                # graph.add((same_as_uri, RDF.type, SDO.URL))

            if (e_name := entry.get('name')) is not None:
                graph.add((polypeptide, SDO.name, Literal(e_name)))


def map_general_references(drug_xml, graph, entity, starting_tag='drugbank:general-references'):
    references = {
        'articles': [{
            'ref-id': get_element_text(article, 'drugbank:ref-id'),
            'pubmed-id': get_element_text(article, 'drugbank:pubmed-id'),
            'citation': get_element_text(article, 'drugbank:citation'),
        } for article in drug_xml.findall(starting_tag + '/drugbank:articles/drugbank:article', ns)],
        'textbooks': [{
            'ref-id': get_element_text(textbook, 'drugbank:ref-id'),
            'isbn': get_element_text(textbook, 'drugbank:isbn'),
            'citation': get_element_text(textbook, 'drugbank:citation'),
        } for textbook in drug_xml.findall(starting_tag + '/drugbank:textbooks/drugbank:textbook', ns)],
        'links': [{
            'ref-id': get_element_text(link, 'drugbank:ref-id'),
            'title': get_element_text(link, 'drugbank:title'),
            'url': get_element_text(link, 'drugbank:url'),
        } for link in drug_xml.findall(starting_tag + '/drugbank:links/drugbank:link', ns)],
        'attachments': [{
            'ref-id': get_element_text(attachment, 'drugbank:ref-id'),
            'title': get_element_text(attachment, 'drugbank:title'),
            'url': get_element_text(attachment, 'drugbank:url'),
        } for attachment in drug_xml.findall(starting_tag + '/drugbank:attachments/drugbank:attachment', ns)]
    }

    for article in references.get('articles'):
        if article.get('pubmed-id') is not None:
            uri = URIRef(f'https://pubmed.ncbi.nlm.nih.gov/{article.get("pubmed-id")}/')
            graph.add((uri, RDF.type, SDO.CreativeWork))
            graph.add((uri, SDO.about, entity))
            graph.add((entity, SDO.subjectOf, uri))

            # if article.get('ref-id') is not None:
            #     graph.add((uri, SDO.identifier, Literal(article.get('ref-id'))))
            if article.get('citation') is not None:
                graph.add((uri, SDO.citation, Literal(article.get('citation'))))

    for textbook in references.get('textbooks'):
        if textbook.get('isbn') is not None:
            uri = URIRef(f'{textbook.get("isbn").replace(" ", "-")}/')
            graph.add((uri, RDF.type, SDO.CreativeWork))
            graph.add((uri, SDO.about, entity))
            graph.add((entity, SDO.subjectOf, uri))

            # if textbook.get('ref-id') is not None:
            #     graph.add((uri, SDO.identifier, Literal(textbook.get('ref-id'))))
            if textbook.get('citation') is not None:
                graph.add((uri, SDO.citation, Literal(textbook.get('citation'))))

    for link in references.get('links'):
        if (url := link.get('url')) is not None:
            if str(url).startswith("//s3"):
                url = "https:" + url
            uri = URIRef(url)
            graph.add((entity, SDO.sameAs, uri))
            # graph.add((uri, RDF.type, SDO.URL))
            # if link.get('ref-id') is not None:
            #     graph.add((uri, SDO.identifier, Literal(link.get('ref-id'))))

    for attachment in references.get('attachments'):
        if (url := attachment.get('url')) is not None:
            if str(url).startswith("//s3"):
                url = "https:" + url
            uri = URIRef(url)
            graph.add((entity, SDO.sameAs, uri))
            # graph.add((uri, RDF.type, SDO.URL))
            # if attachment.get('ref-id') is not None:
            #     graph.add((uri, SDO.identifier, Literal(attachment.get('ref-id'))))


def map_polypeptides(drug_xml, graph, drug):
    all_polypeptides = drug_xml.findall('drugbank:targets/drugbank:target/drugbank:polypeptide', ns) + drug_xml.findall(
        'drugbank:enzymes/drugbank:enzyme/drugbank:polypeptide', ns) + drug_xml.findall(
        'drugbank:carriers/drugbank:carrier/drugbank:polypeptide', ns) + drug_xml.findall(
        'drugbank:transporters/drugbank:transporter/drugbank:polypeptide', ns)

    polypeptides = [{
        'id': polypeptide.attrib.get('id'),
        'name': get_element_text(polypeptide, 'drugbank:name'),
        'general-function': get_element_text(polypeptide, 'drugbank:general-function'),
        'specific-function': get_element_text(polypeptide, 'drugbank:specific-function'),
        'gene-name': get_element_text(polypeptide, 'drugbank:gene-name'),
        'gene-id': next((get_element_text(polypeptide, 'drugbank:identifier')
                         for identifier in
                         polypeptide.findall('drugbank:external-identifiers/drugbank:external-identifier', ns)
                         if identifier.find('drugbank:resource', ns).text == "HUGO Gene Nomenclature Committee (HGNC)"),
                        None),
        'cellular-location': get_element_text(polypeptide, 'drugbank:cellular-location'),
        'molecular-weight': get_element_text(polypeptide, 'drugbank:molecular-weight'),
        'sequence': get_element_text(polypeptide, 'drugbank:amino-acid-sequence'),
        'synonyms': [synonym.text for synonym in polypeptide.findall('drugbank:synonyms/drugbank:synonym', ns) if
                     synonym is not None],
        'identifiers': [{
            'resource': get_element_text(identifier, 'drugbank:resource'),
            'id': get_element_text(identifier, 'drugbank:identifier'),
        } for identifier in polypeptide.findall('drugbank:external-identifiers/drugbank:external-identifier', ns)]
    } for polypeptide in all_polypeptides]

    if not polypeptides:
        return

    for p in polypeptides:
        if p.get('id') is not None:
            polypeptide = URIRef(drug_ext_ids_uris.get('polypeptides').replace("IDENTIFIER", p.get('id')))
            graph.add((polypeptide, RDF.type, SDO.Protein))
            graph.add((polypeptide, SDO.additionalType, Literal("Polypeptide")))
            graph.add((polypeptide, SDO.identifier, Literal(p.get('id'))))
            graph.add((drug, SDO.bioChemInteraction, polypeptide))
            graph.add((polypeptide, SDO.bioChemInteraction, drug))

            if p.get('name') is not None:
                graph.add((polypeptide, SDO.name, Literal(p.get('name'))))
            if p.get('general-function') is not None:
                graph.add((polypeptide, SDO.description, Literal(p.get('general-function'))))
            if p.get('specific-function') is not None:
                graph.add((polypeptide, SDO.description, Literal(p.get('specific-function'))))

            for s in p.get('synonyms'):
                if s is not None:
                    graph.add((polypeptide, SDO.alternateName, Literal(s)))

            for id_pack in p.get('identifiers'):
                if (resource := id_pack.get('resource')) is None:
                    continue


                if link := drug_ext_ids_uris.get(resource):
                    p_id = id_pack.get('id')

                    if p_id is not None and resource != "HUGO Gene Nomenclature Committee (HGNC)":
                        graph.add((polypeptide, SDO.sameAs, URIRef(link.replace("IDENTIFIER", p_id))))

            if p.get('sequence') is not None:
                graph.add((polypeptide, SDO.hasBioPolymerSequence, Literal(p.get('sequence'))))

            if p.get('gene-id') is not None:
                gene_uri = URIRef(drug_ext_ids_uris.get("HUGO Gene Nomenclature Committee (HGNC)").replace("IDENTIFIER", p.get('gene-id')))
                graph.add((polypeptide, SDO.isEncodedByBioChemEntity, gene_uri))
                graph.add((gene_uri, SDO.encodesBioChemEntity, polypeptide))
                graph.add((gene_uri, RDF.type, SDO.Gene))
                if p.get('gene-name') is not None:
                    graph.add((gene_uri, SDO.name, Literal(p.get('gene-name'))))

            # graph.add((polypeptide, SDO.isLocatedInSubcellularLocation, p.get('cellular-location'))) TREBA DEFINED TERM


def map_genes(drug_xml, graph, drug):
    all_polypeptides = drug_xml.findall('drugbank:targets/drugbank:target/drugbank:polypeptide', ns) + drug_xml.findall(
        'drugbank:enzymes/drugbank:enzyme/drugbank:polypeptide', ns) + drug_xml.findall(
        'drugbank:carriers/drugbank:carrier/drugbank:polypeptide', ns) + drug_xml.findall(
        'drugbank:transporters/drugbank:transporter/drugbank:polypeptide', ns)

    genes = [{
        'name': get_element_text(polypeptide, 'drugbank:gene-name'),
        'id': next((identifier.find('drugbank:identifier', ns).text
                    for identifier in
                    polypeptide.findall('drugbank:external-identifiers/drugbank:external-identifier', ns)
                    if identifier.find('drugbank:resource', ns).text == "HUGO Gene Nomenclature Committee (HGNC)"),
                   None),
        'locus': get_element_text(polypeptide, 'drugbank:locus'),
        'polypeptide': polypeptide.attrib.get('id'),
        'sequence': get_element_text(polypeptide, 'drugbank:gene-sequence'),
    } for polypeptide in all_polypeptides]

    if not genes:
        return

    for g in genes:
        g_id = g.get('id')

        if g_id is not None:
            g_uri = drug_ext_ids_uris.get("HUGO Gene Nomenclature Committee (HGNC)").replace("IDENTIFIER", g_id)
            gene = URIRef(g_uri)

            graph.add((gene, RDF.type, SDO.Gene))
            graph.add((gene, SDO.identifier, Literal(g_id)))
            graph.add((gene, SDO.additionalType, Literal("Gene")))
            if g.get('name') is not None:
                graph.add((gene, SDO.name, Literal(g.get('name'))))
            if g.get('sequence') is not None:
                graph.add((gene, SDO.hasBioPolymerSequence, Literal(g.get('sequence'))))
            if g.get('polypeptide') is not None:
                polypeptide_uri = URIRef(
                    drug_ext_ids_uris.get('polypeptides').replace("IDENTIFIER", g.get('polypeptide')))
                graph.add((gene, SDO.encodesBioChemEntity, polypeptide_uri))
                graph.add((polypeptide_uri, SDO.isEncodedByBioChemEntity, gene))


def map_snp_effects(drug_xml, graph, drug, drug_id):
    all_effects = [{
        'protein-name': get_element_text(effect, 'drugbank:protein-name'),
        'gene-name': get_element_text(effect, 'drugbank:gene-symbol'),
        'uniprot-id': get_element_text(effect, 'drugbank:uniprot-id'),
        'rs-id': get_element_text(effect, 'drugbank:rs-id'),
        'allele': get_element_text(effect, 'drugbank:allele'),
        'defining-change': get_element_text(effect, 'drugbank:defining-change'),
        'description': get_element_text(effect, 'drugbank:description'),
        'pubmed-id': get_element_text(effect, 'drugbank:pubmed-id'),
    } for effect in drug_xml.findall('drugbank:snp-effects/drugbank:effect', ns)]

    if not all_effects:
        return

    for effect in all_effects:
        if (p_id := effect.get('uniprot-id')) is not None:
            protein = URIRef(drug_ext_ids_uris.get('polypeptides').replace("IDENTIFIER", p_id))

            graph.add((protein, RDF.type, SDO.Protein))
            graph.add((protein, SDO.name, Literal(effect.get('protein-name'))))
            graph.add((drug, SDO.bioChemInteraction, protein))
            graph.add((protein, SDO.bioChemInteraction, drug))

            if (g_name := effect.get('gene-name')) is not None:

                if (rs_id := effect.get('rs-id')) is not None:
                    snp_uri = URIRef(drug_ext_ids_uris.get('rs').replace("IDENTIFIER", rs_id))

                    graph.add((snp_uri, RDF.type, SDO.DefinedTerm))
                    graph.add((snp_uri, SDO.inDefinedTermSet, snp_effects))
                    graph.add((snp_uri, SDO.additionalType, Literal("SNP")))
                    graph.add((snp_effects, SDO.hasDefinedTerm, snp_uri))
                    graph.add((drug, SDO.isInvolvedInBiologicalProcess, snp_uri))
                    graph.add((protein, SDO.isInvolvedInBiologicalProcess, snp_uri))

                    if (desc := effect.get('description')) is not None:
                        graph.add((snp_uri, SDO.description, Literal(desc)))
                    if (pubmed_id := effect.get('pubmed-id')) is not None:
                        graph.add((URIRef(f'https://pubmed.ncbi.nlm.nih.gov/{pubmed_id}/'), RDF.type, SDO.CreativeWork ))
                        graph.add((snp_uri, SDO.subjectOf, URIRef(f'https://pubmed.ncbi.nlm.nih.gov/{pubmed_id}/')))

                    disambiguating_description = f"{drug_id} {p_id} {g_name}\n"
                    disambiguating_description += f"Allele: {allele}\n" if (allele := effect.get(
                        'allele')) is not None else ''
                    disambiguating_description += f"Defining change: {change}" if (change := effect.get(
                        'defining-change')) is not None else ''

                    graph.add((snp_uri, SDO.disambiguatingDescription, Literal(disambiguating_description)))


def map_all_classifications(graph):
    classyfire = open('./classyfire.json')
    json_array = json.load(classyfire)

    classification_tree = {}
    for item in json_array:
        item_name = str(item['name']).lower().capitalize()
        chemont_id = str(item['chemont_id']).replace("HEMONTID:", "")
        parent_id = str(item['parent_chemont_id']).replace("HEMONTID:", "")

        classification_tree[item_name] = chemont_id

        parent_uri = f'{classyfire_base_uri}/{parent_id}' if parent_id != 'None' else classyfire_base_uri

        item_uri = f'{classyfire_base_uri}/{chemont_id}'

        graph.add((URIRef(parent_uri), RDF.type, SDO.Taxon))
        if parent_id == 'None':
            graph.add((URIRef(parent_uri), SDO.taxonRank, Literal('Root')))

        graph.add((URIRef(item_uri), RDF.type, SDO.Taxon))
        graph.add((URIRef(item_uri), SDO.name, Literal(item_name)))

        graph.add((URIRef(parent_uri), SDO.childTaxon, URIRef(item_uri)))
        graph.add((URIRef(item_uri), SDO.parentTaxon, URIRef(parent_uri)))

    return classification_tree


def map_classifications(drug_xml, graph, drug, classification_tree):
    classification = drug_xml.find('drugbank:classification', ns)

    kingdoms = set()
    superclasses = set()
    classes = set()
    subclasses = set()
    direct_parents = set()
    alternative_parents = set()

    if classification is not None:
        description = get_element_text(classification, 'description')

        if description is None:
            return

        graph.add((drug, SDO.taxonomicRange, Literal(description)))

        kingdom = str(classification.find('drugbank:kingdom', ns).text.lower().capitalize()) if get_element_text(
            classification, 'drugbank:kingdom') is not None else None

        superclass = str(classification.find('drugbank:superclass', ns).text.lower().capitalize()) if get_element_text(
            classification, 'drugbank:superclass') is not None else None

        c = classification.find('drugbank:class', ns).text.lower().capitalize() if get_element_text(classification,
                                                                                                    'drugbank:class') is not None else 'None'

        subclass = str(classification.find('drugbank:subclass', ns).text.lower().capitalize()) if get_element_text(
            classification, 'drugbank:subclass') is not None else None

        direct_parent = str(
            classification.find('drugbank:direct-parent', ns).text.lower().capitalize()) if get_element_text(
            classification, 'drugbank:direct-parent') is not None else None

        alt_parents = [alt.text.lower().capitalize() for alt in
                       classification.findall('drugbank:alternative-parent', ns)]

        # substituents = [substituent.text for substituent in classification.findall('drugbank:substituent', ns)]

        if (kingdom is not None) and (kingdom != "None"):
            kingdoms.add(kingdom.lower().capitalize())
            this_id = classification_tree.get(kingdom)
            graph.add((URIRef(f'{classyfire_base_uri}/{this_id}'), SDO.taxonRank, Literal('Kingdom')))
            graph.add((drug, SDO.taxonomicRange, (URIRef(f'{classyfire_base_uri}/{this_id}'))))

        if (superclass is not None) and (superclass != "None"):
            superclasses.add(superclass.lower().capitalize())
            this_id = classification_tree.get(superclass)
            graph.add((URIRef(f'{classyfire_base_uri}/{this_id}'), SDO.taxonRank, Literal('Superclass')))
            graph.add((drug, SDO.taxonomicRange, (URIRef(f'{classyfire_base_uri}/{this_id}'))))

        if (c is not None) and (c != "None"):
            classes.add(c.lower().capitalize())
            this_id = classification_tree.get(c)
            graph.add((URIRef(f'{classyfire_base_uri}/{this_id}'), SDO.taxonRank, Literal('Class')))
            graph.add((drug, SDO.taxonomicRange, (URIRef(f'{classyfire_base_uri}/{this_id}'))))

        if (subclass is not None) and (subclass != "None"):
            subclasses.add(subclass.lower().capitalize())
            this_id = classification_tree.get(subclass)
            graph.add((URIRef(f'{classyfire_base_uri}/{this_id}'), SDO.taxonRank, Literal('Subclass')))
            graph.add((drug, SDO.taxonomicRange, (URIRef(f'{classyfire_base_uri}/{this_id}'))))

        if (direct_parent is not None) and (direct_parent != "None"):
            direct_parents.add(direct_parent.lower().capitalize())
            this_id = classification_tree.get(direct_parent)
            graph.add((drug, SDO.taxonomicRange, (URIRef(f'{classyfire_base_uri}/{this_id}'))))

        if not alt_parents:
            return

        for alt_parent in alt_parents:
            if (alt_parent is not None) and (alt_parent != "None"):
                alternative_parents.add(alt_parent.lower().capitalize())
                this_id = classification_tree.get(alt_parent)
                graph.add((drug, SDO.taxonomicRange, (URIRef(f'{classyfire_base_uri}/{this_id}'))))


# def map_subcellular_locations(filepath, graph):
#     with open(filepath, 'r') as file:
#         reader = csv.DictReader(file, delimiter='\t')
#
#         subcellular_locations = {row['Name'] : {
#             'location_id': row['Subcellular location ID'],
#             'go-id': row['Gene Ontologies'],
#             'synonyms': row['Synonyms'],
#             'description': row['Description'],
#         }
#         for row in reader if row['Category'] == "Cellular component" and row['Gene Ontologies'] is not None}
#
#     go_url = URIRef('https://amigo.geneontology.org/amigo/search/ontology')
#     graph.add((go_url, RDF.type, SDO.DefinedTermSet))
#
#     return subcellular_locations

def map_groups_term_set(graph):
    graph.add((biotech_term_set, RDF.type, SDO.DefinedTermSet))
    graph.add((biotech_term_set, SDO.name, Literal("Biotech")))
    graph.add((biotech_term_set, SDO.description, Literal(
        "Biotech is used for any drug that is derived from living systems or organisms, usually composed of high molecular weight mixtures of protein. Contains elements which denote the groups this biotech drug belongs to. Groups include approved, nutraceutical, illicit, withdrawn, investigational, and experimental.")))

    graph.add((small_molecule_term_set, RDF.type, SDO.DefinedTermSet))
    graph.add((small_molecule_term_set, SDO.name, Literal("Small molecule")))
    graph.add((small_molecule_term_set, SDO.description, Literal(
        "Small molecule describes a low molecular weight organic compound, a drug that can enter cells easily because it has a low molecular weight. Contains elements which denote the groups this small drug belongs to. Groups include approved, nutraceutical, illicit, withdrawn, investigational, and experimental.")))

    for group, value in group_nodes.items():
        graph.add((URIRef(small_molecule_uri + value), RDF.type, SDO.DefinedTerm))
        graph.add((URIRef(small_molecule_uri + value), SDO.name, Literal(group.capitalize())))
        graph.add((URIRef(small_molecule_uri + value), SDO.inDefinedTermSet, small_molecule_term_set))
        graph.add((URIRef(small_molecule_uri + value), SDO.description, Literal(group_descriptions.get(group))))

        graph.add((URIRef(biotech_uri + value), RDF.type, SDO.DefinedTerm))
        graph.add((URIRef(biotech_uri + value), SDO.name, Literal(group.capitalize())))
        graph.add((URIRef(biotech_uri + value), SDO.inDefinedTermSet, biotech_term_set))
        graph.add((URIRef(biotech_uri + value), SDO.description, Literal(group_descriptions.get(group))))

def map_snps_term_set(graph):
    graph.add((snp_effects, RDF.type, SDO.DefinedTermSet))
    graph.add((snp_effects, SDO.name, Literal("Single nucleotide polymorphism effects")))
    graph.add((snp_effects, SDO.alternateName, Literal("SNP effects")))
    graph.add((snp_effects, SDO.description, Literal(
        "A set of single nucleotide polymorphisms (SNPs) relevent to drug activity or metabolism, and the effects these may have on pharmacological activity. A single nucleotide polymorphism (abbreviated SNP, pronounced snip) is a genomic variant at a single base position in the DNA. Scientists study if and how SNPs in a genome influence health, disease, drug response and other traits. SNP effects in the patient may require close monitoring, an increase or decrease in dose, or a change in therapy.")))

def map_pathways_term_set(graph):
    graph.add((pathways_term_set, RDF.type, SDO.DefinedTermSet))
    graph.add((pathways_term_set, SDO.name, Literal("Pathways")))
    graph.add((pathways_term_set, SDO.description, Literal(
        "Metabolic, disease, and biological pathways that the drug is involved in, as identified by the Small Molecule Pathway Database (SMPDB).")))

def add_metadata(connection):
    drugbank = connection.createURI("https://finki.ukim.mk/drugbank")
    metadata = connection.createURI("https://finki.ukim.mk/metadata")

    dataset = connection.createURI("https://ag1twkresvsgv5ez.allegrograph.cloud")
    creator = connection.createURI("https://www.linkedin.com/in/aleksandar-ristovski-b13030253/")
    publisher = connection.createURI("https://www.finki.ukim.mk/")
    contributor = connection.createURI("https://mjovanovik.com/")
    source = connection.createURI("https://www.drugbank.com/")

    triple_count = connection.size(contexts=[drugbank])
    current_date = date.today()

    prefixes = f"""
        PREFIX dcterms: <http://purl.org/dc/terms/>
        PREFIX void: <http://rdfs.org/ns/void#>
        PREFIX foaf: <http://xmlns.com/foaf/0.1/>
        PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#> 
        PREFIX finki: <https://finki.ukim.mk/> 
        PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>

    """

    delete_query = f"""
            WITH finki:metadata
            DELETE {{ 
                    {dataset} void:triples ?oldTriples;
                              dcterms:modified ?modified.
            }}
            WHERE {{
                    {dataset} void:triples ?oldTriples;
                              dcterms:modified ?modified.
            }}
        """

    connection.executeUpdate(f'{prefixes}{delete_query}')

    check_created_query = f"""        
        WITH finki:metadata
        INSERT {{
          {dataset} a void:Dataset ;            
                    dcterms:created "{current_date}"^^xsd:date .
        }}
        WHERE {{ 
              FILTER( NOT EXISTS {{ {dataset} dcterms:created ?created . }})
        }}
        """
    connection.executeUpdate(f'{prefixes}{check_created_query}')

    update_query = f"""        
        WITH finki:metadata
        INSERT {{
                {dataset} a void:Dataset ;
                    dcterms:title "Drugbank RDF dataset" ;
                    dcterms:description "RDF data extracted from Drugbank dataset" ;
                    dcterms:creator {creator} ;
                    dcterms:contributor {contributor} ;
                    dcterms:publisher {publisher} ;
                    dcterms:source {source} ;
                    void:triples {triple_count} ;
                    dcterms:modified "{current_date}"^^xsd:date.

                {creator} a foaf:Person ;
                    foaf:name "Aleksandar Ristovski" .

                {contributor} a foaf:Person ;
                    foaf:name "Milosh Jovanovikj" .

                {publisher} a foaf:Organization ;
                    rdfs:label "FCSE" ;
                    foaf:name "Faculty of Computer Science and Engineering" .
        }}
        WHERE{{
        }}
    """
    connection.executeUpdate(f'{prefixes}{update_query}')

def number_of_properties(drug_xml, experimental_set, calculated_set):
    calculated_properties = {
        get_element_text(prop, 'drugbank:kind'): get_element_text(prop, 'drugbank:value')
        for prop in drug_xml.findall('drugbank:calculated-properties/drugbank:property', ns)
    }

    experimental_properties = {
        get_element_text(prop, 'drugbank:kind'): get_element_text(prop, 'drugbank:value')
        for prop in drug_xml.findall('drugbank:experimental-properties/drugbank:property', ns)
    }

    for prop, value in calculated_properties.items():
        calculated_set.add(prop)

    for prop, value in experimental_properties.items():
        experimental_set.add(prop)

def main():
    # load_dotenv()
    # connection = ag_connect('diplomska-graph')
    connection = ag_connect('diplomska-graph', create=True, clear=True)

    drugbank = connection.createURI("https://finki.ukim.mk/drugbank")

    start = time.perf_counter()
    context = ET.iterparse('./data/full database.xml')
    _, root = next(context)

    g = Graph()
    g.bind("schema", SDO)

    map_groups_term_set(g)
    map_pathways_term_set(g)
    map_snps_term_set(g)
    classification_tree = map_all_classifications(g)

    # subcellular_locations = map_subcellular_locations("./data/subcellular_locations.tsv", g)

    connection.addData(f"""{g.serialize()}""", context=drugbank)

    i = 0
    counter = 0
    graph = Graph()
    graph.bind("schema", SDO)

    for event, drug_xml in context:
        if (drug_xml.tag == "{http://www.drugbank.ca}drug"
                and drug_xml.attrib.get('type') in ['biotech', 'small molecule'] and event == 'end'):

            drug_id = next(identifier.text
                           for identifier in drug_xml.findall("drugbank:drugbank-id", ns)
                           if identifier.attrib.get('primary') == "true")

            drug = URIRef(f'https://go.drugbank.com/drugs/{drug_id}')



            graph.add((drug, RDF.type, SDO.Drug))
            graph.add((drug, RDF.type, SDO.MolecularEntity))
            graph.add((drug, SDO.identifier, Literal(drug_id)))

            map_name(drug_xml, graph, drug)
            map_description(drug_xml, graph, drug)
            map_groups(drug_xml, graph, drug)

            map_label_details(drug_xml, graph, drug)
            map_toxicity(drug_xml, graph, drug)
            map_monoisotopic_mass(drug_xml, graph, drug)

            map_pharmacodynamics(drug_xml, graph, drug)
            map_route_of_elimination(drug_xml, graph, drug)
            map_absorption(drug_xml, graph, drug)
            map_half_life(drug_xml, graph, drug)
            map_volume_of_distribution(drug_xml, graph, drug)
            map_clearance(drug_xml, graph, drug)
            map_mechanism_of_action(drug_xml, graph, drug)

            map_synonyms(drug_xml, graph, drug)
            map_general_references(drug_xml, graph, drug)
            map_routes(drug_xml, graph, drug)
            map_strengths(drug_xml, graph, drug)
            map_food_warnings(drug_xml, graph, drug)
            map_alcohol_warnings(drug_xml, graph, drug)
            map_drug_interactions(drug_xml, graph, drug)

            map_over_the_counter(drug_xml, graph, drug)
            map_proprietary_names(drug_xml, graph, drug)
            map_external_identifiers(drug_xml, graph, drug)
            map_pdb_entries(drug_xml, graph, drug)
            map_salts(drug_xml, graph, drug)
            map_calculated_properties(drug_xml, graph, drug)

            map_atc_codes(drug_xml, graph, drug)
            map_drug_pathways(drug_xml, graph, drug)
            map_drug_targets_enzymes_carriers_transporters(drug_xml, graph, drug)
            map_polypeptides(drug_xml, graph, drug)
            map_genes(drug_xml, graph, drug)

            map_snp_effects(drug_xml, graph, drug, drug_id)
            map_classifications(drug_xml, graph, drug, classification_tree)

            map_indications(drug_xml, graph, drug)

            i += 1
            if len(graph) > 10000:
                print(f'From {counter} to {counter + i} drugs: {len(graph)}')
                counter += i
                i = 0
                connection.addData(f"""{graph.serialize()}""", context=drugbank)
                graph.remove((None, None, None))

            drug_xml.clear()
            root.clear()
    del context

    if len(graph) != 0:
        print(f'From {counter} to {counter + i} drugs: {len(graph)}')
        connection.addData(f"""{graph.serialize()}""", context=drugbank)

    time.sleep(2)
    add_metadata(connection)
    print(f'Ended in {time.perf_counter() - start}')


if __name__ == '__main__':
    main()

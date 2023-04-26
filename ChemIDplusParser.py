# This file parses the NIH ChemIDplus Subset Data https://www.nlm.nih.gov/databases/download/chemidplus.html
# It works with the zipfile (like chemid-20230222.zip) that could be obtained from here https://ftp.nlm.nih.gov/projects/chemidlease

import os
import re

from lxml import etree
from zipfile import ZipFile



# These functions are class-independent
# Separate placement simplifies subsequent usage of map
def extract_info_list(chemical_xml, infolist_tag):

    """
    This function extracts XML tags from individual chemicals in the NIH ChemIDplus datase
    Input: XML with individual chemical; tag to be extracted (string)
    Output: list of characters for a .tsv file or None

    The function is primarily intended for the internal use
    """
    mol_name = chemical_xml.attrib["displayName"]
    info_list = chemical_xml.findall(".//" + infolist_tag)

    if len(info_list) < 1:
        return None
    
    info_sublist = []

    for element in info_list:
        element_val = element.text
        source_list = element.findall(".//Source")
        sources = []

        for source in source_list:
            sources.append(source.text)
        
        sources_str = ";".join(sources)
        syn_source = mol_name + "\t" + element_val + "\t" + sources_str + "\n"
        info_sublist.append(syn_source)

    return info_sublist

def extract_synonyms_modif(chemical_xml, bracket_terms):

    mol_name = chemical_xml.attrib["displayName"]
    synonyms_list = chemical_xml.findall(".//Synonyms")

    # Generating initial list of synonyms
    if len(synonyms_list) > 1:
        synonyms_list_init = []
        for element_synonym in synonyms_list:
            synonym = element_synonym.text
            synonyms_list_init.append(synonym)
    else:
        synonyms_list_init = []

    if mol_name == "":
        return None
    
    molecules_synonyms = []

    mol_name_clean = mol_name
    for term in bracket_terms:
        mol_name_clean = mol_name_clean.replace(term, "")
    mol_name_clean = mol_name_clean.strip()

    # Inspecting if the synonym is new
    if mol_name_clean != mol_name and mol_name_clean not in synonyms_list_init:
        generated_1st_syn = mol_name + "\t" + mol_name_clean + "\t" + "GENERATED\n"
        molecules_synonyms.append(generated_1st_syn)
        synonyms_list_init.append(mol_name_clean)

    if len(synonyms_list) < 1 and len(molecules_synonyms) < 0:
        return None
    
    if len(synonyms_list) >= 1:
        for element_synonym in synonyms_list:
            synonym = element_synonym.text
            source_list = element_synonym.findall(".//Source")
            sources = []

            for source in source_list:
                sources.append(source.text)
            
            sources_str = ";".join(sources)
            syn_source = mol_name + "\t" + synonym + "\t" + sources_str + "\n"
            
            molecules_synonyms.append(syn_source)

            # Inspecting if the synonym can be generated
            if re.search("\[", string=synonym):
                for term in bracket_terms:
                    synonym_clean = synonym.replace(term, "")
                synonym_clean = synonym.strip()

                # Inspecting if the synonym is new
                if synonym_clean != synonym and synonym_clean not in synonyms_list_init:
                    syn_source_clean = mol_name + "\t" + synonym_clean + "\t" + "GENERATED\n"
                    molecules_synonyms.append(syn_source_clean)
                    synonyms_list_init.append(synonym_clean)

    return molecules_synonyms

class ChemIDplusParser:

    """
    This class contains several functions to parse  NIH ChemIDplus Subset Data https://www.nlm.nih.gov/databases/download/chemidplus.html
    It works with the zipfile (like chemid-20230222.zip) that could be obtained from here https://ftp.nlm.nih.gov/projects/chemidlease
    """

    def __init__(self, database_zip_file):

        # Unzip the downloaded file
        ZipFile(database_zip_file).extractall("extracted_zip")

        # Import XML tree
        for file in os.listdir("extracted_zip"):
            if re.search(pattern=".xml", string=file):
                path = "extracted_zip" + "/" + file
                self.xml_file = etree.parse(path)
                self.chemical_list = self.xml_file.findall("Chemical")
                print("The list contains ", len(self.chemical_list), " compounds")
        
        # Importing bracket terms
        bracket_terms = open("bracket_terms.txt", "r")
        bracket_terms = bracket_terms.read()
        bracket_terms = str.split(bracket_terms, sep="\n")
        self.bracket_terms = bracket_terms

    def extract_all_classes(self, file):

        class_list_param = ["ClassificationCode"]*len(self.chemical_list)
        extacted_classes = list(map(extract_info_list, self.chemical_list, class_list_param))

        # Writing output to a .tsv file
        output_file = open(file, "w")
        header = "Molecule" + "\t" + "Class" + "\t" + "Class_source\n"
        output_file.write(header)
        [output_file.write(drug_class) for element in extacted_classes if element is not None for drug_class in element]
        output_file.close()


    def extract_all_synonyms(self, file, generate_synonyms=True):

        if generate_synonyms:
            bracket_term_lists = [self.bracket_terms]*len(self.chemical_list)
            chem_list = self.chemical_list
            extacted_syn = list(map(extract_synonyms_modif, chem_list, bracket_term_lists))

        else:
            syn_list_param = ["Synonyms"]*len(self.chemical_list)
            chem_list = self.chemical_list
            extacted_syn = list(map(extract_info_list, chem_list, syn_list_param))

        # Writing output to a .tsv file
        output_file = open(file, "w")
        header = "Molecule" + "\t" + "Synonym" + "\t" + "Synonym_source\n"
        output_file.write(header)
        [output_file.write(synonym) for element in extacted_syn if element is not None for synonym in element]
        output_file.close()
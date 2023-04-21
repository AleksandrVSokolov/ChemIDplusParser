# Get datasets from ChemIDplus

from ChemIDplusParser import ChemIDplusParser

parser = ChemIDplusParser("/home/aleksandr/Desktop/WORK/ChemIDPlusDataParserPY/chemid-20230222.zip")

# Extracting lists
parser.extract_all_classes("classes_chemIDplus.txt")
parser.extract_all_synonyms("synonyms_chemIDplus.txt")
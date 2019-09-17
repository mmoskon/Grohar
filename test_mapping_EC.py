"""
# http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc145
from Bio import ExPASy
from Bio import SwissProt

accessions = ["O23729", "O23730", "O23731"]
records = []

for accession in accessions:
    handle = ExPASy.get_sprot_raw(accession)
    record = SwissProt.read(handle)
    EC = record.description.split(";")[1].split("=")[1]
    records.append(EC)
    
"""

"""
import mygene

mg = mygene.MyGeneInfo()
g = mg.getgene(100759423)
"""

#https://www.biostars.org/p/104733/

"""
# dela, ampak zelo pocasi...
kegg.find("genes", '100759423')
#Out[386]: 'cge:100759423\tK18311 N-acetylaspartylglutamate/N-acetylaspartylglutamylglutamate synthase [EC:6.3.2.41 6.3.2.42] | (RefSeq) Rimkla; ribosomal modification protein rimK like family member A\n'
"""


#http://nbviewer.jupyter.org/url/pythonhosted.org//bioservices/_downloads/Entrez_EUtils.ipynb
from bioservices import EUtils
e = EUtils()
db = 'gene'
id_list = '100759423'
results = e.EPost(db, id_list)



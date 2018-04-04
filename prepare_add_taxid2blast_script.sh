#! /bin/bash -x
#mkdir accession2taxid
#cd accession2taxid
#rm nucl*
#wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl*gz
#gunzip *gz
#split -l 50000000 nucl_gb.accession2taxid nucl_gb.accession2taxid.
#split -l 50000000 nucl_est.accession2taxid nucl_est.accession2taxid.
#split -l 50000000 nucl_gss.accession2taxid nucl_gss.accession2taxid.
#split -l 50000000 nucl_wgs.accession2taxid nucl_wgs.accession2taxid.
#cd ..

path=$(pwd)"/accession2taxid/gi_taxid_nucl.dmp.*"

sed -i.bak 's|/PATH/to/FILE/accession2taxid/gi_taxid_nucl.dmp.\*|'$path'|g' add_taxid2blast.py 

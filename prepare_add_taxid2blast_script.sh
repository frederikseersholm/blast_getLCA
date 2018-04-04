#! /bin/bash -x
mkdir Accession2taxid
cd Accession2taxid
rm nucl*
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl*gz
gunzip *gz
split -l 50000000 nucl_gb.accession2taxid nucl_gb.accession2taxid.
split -l 50000000 nucl_est.accession2taxid nucl_est.accession2taxid.
split -l 50000000 nucl_gss.accession2taxid nucl_gss.accession2taxid.
split -l 50000000 nucl_wgs.accession2taxid nucl_wgs.accession2taxid.
cd ..

path=$(pwd)"/Accession2taxid/gi_taxid_nucl.dmp.*"

sed -i.bak 's|/PATH/to/FILE/Accession2taxid/nucl_\*.accession2taxid.??|'$path'|g' add_taxid2blast.py 

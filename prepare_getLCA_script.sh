#! /bin/bash -x
mkdir taxdump 
wget ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz 
tar -xzf taxdump.tar.gz -C taxdump 
rm taxdump.tar.gz

nodes=$(pwd)"/taxdump/nodes.dmp"
names=$(pwd)"/taxdump/names.dmp"

sed -i.bak 's|/PATH/to/FILE/taxdump/names.dmp|'$names'|g' blast_getLCA.py
sed -i.bak 's|/PATH/to/FILE/taxdump/nodes.dmp|'$nodes'|g' blast_getLCA.py

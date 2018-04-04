# blast_getLCA

## Guide to run the program blast_getLCA.py

#### 1) Download and unpack the blast_getLCA master
```
wget https://github.com/frederikseersholm/blast_getLCA/archive/master.tar.gz
tar xpvf master.tar.gz
cd blast_getLCA-master
```
#### 2) Run the script prepare_getLCA_script.sh
This script downloads the NCBI taxonomy files nodes.dmp and names.dmp to the folder 'taxdump', and adds the paths to the files to the python script.
```
bash prepare_getLCA_script.sh
```
#### 3) Run the script on a test blast-file, modified for the blast_getLCA algorithm
This script downloads the NCBI taxonomy files nodes.dmp and names.dmp to the folder 'taxdump', and adds the paths to the files to the python script.
```
python blast_getLCA.py test.taxid.blast
```


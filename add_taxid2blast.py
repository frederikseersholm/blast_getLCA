#download databases

#wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl*gz
import optparse
import sys

import glob
def main():
   
    p = optparse.OptionParser()
    p.add_option('--outfile', '-o')
    options, arguments = p.parse_args()
    blastfile=arguments[0]
    if options.outfile is None:
        outfile=blastfile.replace('.blast','.taxid.blast')
        print 'No output file specified, writing to: '+outfile
    else:
        outfile=options.outfile
    
    blastfile2=open(blastfile,'r')
    outfile=open(outfile,'w')
    
###################
    count=0
    gi_list=[]
    

    for line in blastfile2.readlines():
        text=line.split()
        gi_tab=text[2]
        if 'LOCALDB' in gi_tab:
            gi_nr=gi_tab
        else:
            gi_nr=gi_tab.split('|')[1]
        if gi_nr not in gi_list:
            gi_list.append(gi_nr)
    blastfile2.close()
            
###################
    gi_dict={}
    for gi_nr in gi_list:
        gi_dict[gi_nr]='TAXID_NOT_FOUND'
####################
    count=0
    filelist=glob.glob('/Users/frederikseersholm1/Frederik/PhD/tools/forgithub/test/blast_getLCA-master/accession2taxid/gi_taxid_nucl.dmp.*')

    for name in filelist:
        print name
    gi2taxid={}
    for gi_name in filelist:
        count+=1
        if count>1:
            break
        print '\nprocessing file '+str(count)+' of '+str(len(filelist))+'\n'
        gi=open(gi_name,'r')
        
        lines_processed = 0
        for line in gi.readlines():
            lines_processed = lines_processed + 1
            if (lines_processed % 5000000 == 0):
                sys.stdout.write('-')
            
            text=line.split()
            gi2taxid[text[0]]=text[1]
        
        gi.close()
        gi2taxid['110190A']='110190A'
        for gi_nr in gi_list:
            if gi_dict[gi_nr]=='TAXID_NOT_FOUND':
                if 'LOCALDB' in gi_nr:
                    gi_dict[gi_nr]=gi_nr.replace('LOCALDB','')
                    continue
                try:
                    gi_dict[gi_nr]=gi2taxid[gi_nr]
                
                except:
                    gi_dict[gi_nr]='TAXID_NOT_FOUND'
        
    print '\n"gi_taxid_nucl.dmp"-files loaded into memory \n'    
    count=0
     

##########################PARSE BLAST FILE 2#########################
    blastfile3=open(blastfile,'r')
    for line in blastfile3.readlines():
        #print line
        #print text[1]
        text=line.split()
        gi_tab=text[2]
        
        if 'LOCALDB' in gi_tab:
            gi_nr=gi_tab
        else:
            gi_nr=gi_tab.split('|')[1]    
        try:
            taxid=gi_dict[gi_nr]
        except:
            taxid='TAXID_NOT_FOUND'
        print gi_nr,taxid
        text[1]=taxid+'|'+text[1]
       
        outfile.write('\t'.join(text)+'\n')
    blastfile3.close()
    outfile.close()

#####################
if __name__ == '__main__':
    main()
    







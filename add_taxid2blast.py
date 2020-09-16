
import optparse
import sys

import glob
def main():
    filelist=glob.glob('/PATH/to/FILE/Accession2taxid/nucl_*.accession2taxid.??')
    p = optparse.OptionParser()
    p.add_option('--outfile', '-o')
    options, arguments = p.parse_args()
    blastfile=arguments[0]
    if options.outfile is None:
        outfile=blastfile.replace('.blast','.taxid.blast')
        print('No output file specified, writing to: '+outfile)
    else:
        outfile=options.outfile
    
    blastfile2=open(blastfile,'r')
    outfile=open(outfile,'w')
    
###################
    count=0
    acc_list=[]
    

    for line in blastfile2.readlines():
        
        text=line.split()
        acc_nr=text[1]
        if acc_nr not in acc_list:
            acc_list.append(acc_nr)
    blastfile2.close()
            
###################
    acc_dict={}
    for acc_nr in acc_list:
        acc_dict[acc_nr]='TAXID_NOT_FOUND'
####################


    count=0
    
    
    acc2taxid={}
    for acc_name in filelist:
        count+=1
        #if count>1:
        #    break
        print('\nprocessing file '+str(count)+' of '+str(len(filelist))+'\n')
        acc=open(acc_name,'r')
        acc2taxid={}
        
        lines_processed = 0
        for line in acc.readlines():
            lines_processed = lines_processed + 1
            if (lines_processed % 1000000 == 0):
                sys.stdout.write('-')
            
            text=line.split()
            acc2taxid[text[0]]=text[2]
        
        acc.close()
        for acc_nr in acc_list:

            if acc_dict[acc_nr]=='TAXID_NOT_FOUND':
                if 'LOCALDB' in acc_nr:
                    acc_dict[acc_nr]=acc_nr.split('LOCALDB')[len(acc_nr.split('LOCALDB'))-1]
                    continue
                try:
                    acc_dict[acc_nr]=acc2taxid[acc_nr]
                
                except:
                    acc_dict[acc_nr]='TAXID_NOT_FOUND'
##########################PARSE BLAST FILE 2#########################
    blastfile3=open(blastfile,'r')
    for line in blastfile3.readlines():
        
        text=line.split()
        try:
            taxid=acc_dict[text[1]]
        except:
            taxid='TAXID_NOT_FOUND'
        text[1]=taxid+'|'+text[1]
        outfile.write('\t'.join(text)+'\n')
    blastfile3.close()
    outfile.close()

#####################
if __name__ == '__main__':
    main()

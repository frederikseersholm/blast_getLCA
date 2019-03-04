import optparse
import sys
sys.path.insert(0, '/Users/frederikseersholm1/Frederik/PhD/tools/blast_getLCA-master') 

from blast_getLCA import *

def find_tree(taxid_list):
    
    int_nodes=['1']
    output='root;'
    cont=0
    
    while cont<2:
        groups=int_nodes
        int_nodes=[]
        #print '\nnew while round'
        
        for group in groups:
            
            taxids_subset=[taxid for taxid in taxid_list if group in find_parents(taxid)]
            #print 'current group: '+group+'  n species in group:'+str(len(taxids_subset))
            
            node=taxidlist2LCA(taxids_subset)
            #print 'LCA: ',node
            nodes2=[[find_parents(taxid)[i-1] for i,parent in enumerate(find_parents(taxid)) if parent==node][0] for taxid in taxids_subset]  
            u_nodes2=list(set(nodes2))
            
            
            if len(u_nodes2)==1:
                u_nodes2=taxids_subset
            
            #print 'unique nodes ',u_nodes2
            
            
            [int_nodes.append(i) for i in u_nodes2]
            #print 'int_nodes',int_nodes
    
            join_str='('+','.join([name[i] for i in u_nodes2])+')'
            #print join_str, name[group]
            output=output.replace(name[group],join_str)
            #print output
        if len(int_nodes)>=len(taxid_list):
            cont+=1 #this step ensures that, once all nodes have been resolved a last loop is run to change the names of all nodes to the lowest possible
    
    return(output)

################################################################################################


#########################    MAIN PROGRAM    #######################


################################################################################################
def main():
   
    p = optparse.OptionParser()
    p.add_option('--outfile', '-o')
    p.add_option('--ignoretaxid', '-i',dest="wrongtax",default='no_ignore',help="list of taxids to ignore, separated by comma")
    p.add_option('--ignoreoffspring', '-g',dest="wronggenus",default='no_ignore',help="list of taxids in which all lower order taxa are ignored")
    p.add_option("-r", "--rank", dest="rank", default='no_drop', help="Lowest rank to include in tree")
    options, arguments = p.parse_args()
    infile=arguments[0]
    if options.outfile is None:
        outfile=infile.replace('.getLCA.tsv','.tre')
        outfile=outfile.split('/')[len(outfile.split('/'))-1]
        print 'No output file specified, writing to: '+outfile
    else:
        outfile=options.outfile
    
    infile=open(infile,'r')
    outfile=open(outfile,'w')
    prev_taxid=''
    taxids=[]
    for line in infile.readlines():
        taxid=line.split()[9]
        
        if taxid=='NOT_FOUND':
            continue

        if options.wrongtax!='no_ignore':
            if taxid in options.wrongtax.split(','):
                continue

        if options.rank!='no_drop':
            taxid=drop_to_level2(taxid,options.rank)
            if taxid=='NOT_FOUND':
                continue
        
        if options.wronggenus!='no_ignore':
            if sum([i in find_parents(taxid) for i in options.wronggenus.split(',')])>0:
                print 'ignored'
                continue
        
        if taxid==prev_taxid:
            continue
        
        taxids.append(taxid)
        prev_taxid=taxid
    ##############FILTER OUT TAXIDS CONTAINED IN OTHER TAXIDS (higher order parents)#############
    parents=[]
    
    [[parents.append(taxid_list) for taxid_list in find_parents(taxid)] for taxid in taxids]
    taxids=[taxid for taxid in taxids if parents.count(taxid)==1]
    #taxids=[taxid for taxid in taxids if taxid !='9823']
    
    #print [name[i] for i in taxids]

    ###################
    tree=find_tree(taxids)
    #print tree
    outfile.write(tree+'\n')

#####################
if __name__ == '__main__':
    main()
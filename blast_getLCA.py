import fileinput
import optparse

names=open('/PATH/to/FILE/taxdump/names.dmp','r')
nodes=open('/PATH/to/FILE/taxdump/nodes.dmp','r')


def get_LCA_from_blast(blastlines,idthreshold,limits):
    nms=[]
    threshold=[]
    ids=[]
    idps=[]
    lseqs=[]
    names=[]
    gap_mm=[]
##################sort lines by highest Nm score###########
    for line in blastlines:
        text=line.split()
        length=text[4]
        al_len=text[5]
        mm=text[6]
        gaps=text[8]
        iden=text[11]
        
        nm=max([int(length),int(al_len)])-int(iden)
        nms.append(nm)
     
    blastlines.sort(key=dict(zip(blastlines, nms)).get)
    nms.sort()
    threshold=nms[0]
#########################append tax ids with Nm score over 'threshold'#################
    for line,nm in zip(blastlines,nms):
        text=line.split()
        length=text[4]
        mm=text[6]
        gaps=text[8]
        idp=round((float(length)-(float(nm)))*100/float(length),2)
        taxid=text[1].split('|')[0]
        
        if nm>threshold:
            break
        ids.append(taxid)
        names.append(text[1])
        gap_mm.append('_'.join([gaps,mm,length]))
        idps.append(str(idp))
        
    if not all([i==idps[0] for i in idps]):
        idps.append('NOT_ALL_MATCH')
####################find LCA if more than 1 id has been accepted################ 
    try:
        lca_id=taxidlist2LCA(ids)
    except:
        lca_id='NOT_FOUND'

################## Assign to higher order taxa based on ID-percent thresholds (limits) ################
    drop='Not_dropped'
    idp=float(idps[0])
    
    if limits[0]>idp and idp>float(idthreshold):
        
        which_drop=[limits[0]>idp>=limits[1],limits[1]>idp>=limits[2],limits[2]>idp]
        drop_level=[i for i,j in zip(['genus','family','order'],which_drop) if j]
        drop='Dropped2'+drop_level[0]+name[lca_id].replace(' ','_')
        lca_id=drop_to_level2(lca_id,drop_level[0])

####################Filter out based on ID-threshold################
    lca=':'.join(find_parents_w_rank_short(lca_id)).replace(' ','_')
    nm=nms[0]/float(length) 
    if (nm>(1-float(idthreshold)/100)):#Similarity threshold set to "threshold" [default=95%]
        lca='NOMATCH_similarity_below_'+str(float(idthreshold)).replace('.0','')+'%'+lca
        lca_id='NOT_FOUND'

   
##################output line###############  
    stats='tothits:'+str(len(blastlines))+'_accepted-hits:'+str(len(ids))+'_Min-Nm:'+str(nms[0])+'_IDp:'+str(idps[0])
    return('\t'.join([text[0],lca,get_rank(lca_id).replace(' ','_'),':'.join(set(ids)),stats,length,':'.join(set(idps)),':'.join(set(gap_mm)),drop,lca_id])+'\n')

################################################################################################


#########################    LOAD names.dmp and nodes.dmp INTO MEMORY    #######################


################################################################################################
name={}
id_from_name={}
parent={}
rank={}
gi2taxid={}
embl={}

lines_processed = 0
for line in names.readlines():
    lines_processed = lines_processed + 1
    if (lines_processed % 500000 == 0):
         print('names.dmp: processing line ' + str(lines_processed))
         
    text=line.replace('\t','').replace('\n','').split('|')
    text=text[0:4]
    
    id_from_name[text[1]]=text[0]
    if text[3]=='scientific name':
        name[text[0]]=text[1]

        

lines_processed = 0        
for line in nodes.readlines():
    lines_processed = lines_processed + 1
    if (lines_processed % 500000 == 0):
         print('nodes.dmp: processing line ' + str(lines_processed))
    text=line.replace('\t','').replace('\n','').split('|')
    parent[text[0]]=text[1]
    rank[text[0]]=text[2]
    embl[text[0]]=text[3]
    

################################################################################################


#########################    LCA FUNCTIONS    #######################


################################################################################################

def find_rankofparents(current_taxid):
    parents=[] 
    found = False
    while found == False:
        parents.append(rank[current_taxid])
        if (current_taxid == '1'):
            return(parents)
            found = True      
        else:
            current_taxid = parent[current_taxid]
def name1(taxid):
    textname=name[taxid]
    return(textname)
            
def get_rank(taxid):
    if taxid!='NOT_FOUND':
                
        try:
            rank1=rank[taxid]
        except:
            rank1='rank_not_found'
    else:
        rank1=taxid
    return(rank1)


def find_parents(current_taxid):
    parents=[] 
    #continue searching until loop reaches 'root' in tree of life
    found = False
    while found == False:
        parents.append(current_taxid)
        if (current_taxid == '1'):
            return(parents)
            found = True      
        else:
            current_taxid = parent[current_taxid]

def find_parents_w_rank(current_taxid):
    a=find_parents(current_taxid)
    b=find_rankofparents(current_taxid)
    return([name[i]+';'+j for i,j in zip(a,b)])

def find_parents_w_rank_short(current_taxid):
    a=find_parents(current_taxid)
    b=find_rankofparents(current_taxid)
    output=[name[i]+';'+j for i,j in zip(a,b) if j in ['subspecies','species','genus','family','suborder','order','superorder','class']]
    
    if name[a[0]]+';'+b[0] not in output:
        output=[name[a[0]]+';'+b[0]]+output
    return(output)

def find_parents_smartsort(current_taxid,org_name):
    a=find_parents(current_taxid)
    b=find_rankofparents(current_taxid)
    output=[]
    for rank in ['subspecies','species','genus','subfamily','family','suborder','order','superorder','class']:
        if rank in b:
            outputname=[name[i]+';'+j for i,j in zip(a,b) if j==rank]
            output.append(outputname[0])
        else:
            output.append('AA;'+rank)

    output.append(org_name)
    return(output)

def find_genus(taxid):
    return(drop_to_level2(taxid,'genus'))
    
def drop_to_level2(taxid,level):
    
    try:
        newtaxid=[taxid2 for taxid2 in find_parents(taxid) if rank[taxid2]==level]
        if len(newtaxid)==0:
            newtaxid=[taxid]
        return(newtaxid[0])
    except:
        return('NOT_FOUND')
  
    

def find_LCA(taxid1,taxid2):
    for i in find_parents(taxid1):
        if i in find_parents(taxid2):
            return(i)
            break

def taxidlist2LCA(taxid_list2):
    count=0
    prev_LCA_id=[]
    taxid_list=[i for i in taxid_list2 if 'TAXID_NOT_FOUND' not in i]
    for taxid in taxid_list:
        count+=1
        
        if count==1:
            prev_LCA_id=taxid
            continue

        else:
            
            try:
                prev_LCA_id=find_LCA(taxid,prev_LCA_id)
            except:
                continue
    return(prev_LCA_id)

def smartsort(getLCA_lines):
    parents=[]
    for line in getLCA_lines:

        taxid=line.split()[9]
        
        try:
            parents.append(find_parents_smartsort(taxid,line))
            
        except:
            parents.append(['AA','AA','AA','AA','AA','AA','AA','AA','AA',line])


    l3=sorted(parents, key = lambda x: (x[8],x[7],x[6],x[5],x[4],x[3],x[2],x[1],x[0]))
    l4=[i[9] for i in l3]
    return(l4)

################################################################################################


#########################    MAIN PROGRAM    #######################


################################################################################################
def main():
############################ Import arguments ##############    
    p = optparse.OptionParser()
    p.add_option('--ignoretaxid', '-i',dest="wrongtax",help="csv file of taxids to ignore, first column should contain taxids")
    p.add_option("-t", "--threshold", dest="threshold", default=95, help="Ignores reads where the best alignment has less than a certain percentage identity to the reference [default=95]")
    p.add_option("-l", "--limits", dest="limits", default='98-95-90', help="List of identity thresholds which will inform the algorithm to drop lowest common ancestor to genus (default: 95-98%), family level (default: 90-95%) or order level (default: below 90%). Separated by dash '-' [default=98-95-90]")
    options, arguments = p.parse_args()
    infiles=arguments
    limits=[int(i) for i in options.limits.split('-')]	
    print 'Identity threshold: '+str(options.threshold)
    print 'Identity limits: '+','.join([str(i) for i in limits])
    #Open file defining taxids to ingore and store them as a list
    wrong_tax=[]
    if options.wrongtax!=None:
        with open(options.wrongtax) as file:
            wrong_tax = [line.split(';')[0].split(',')[0].strip() for line in file]
        wrong_tax=[i for i in wrong_tax if not any(c.isalpha() for c in i)]
    

    #add taxids that should always be ignored
    [wrong_tax.append(taxid) for taxid in ['1749399','155900','37029','32644','419950','547489']]
    
    outlines=[]
############################ Loop over input files ##############
    
    prev_name=[]
    prev_text=[]
    
    for infile in infiles:

        outfile=infile.replace('.blast','')+'.getLCA.tsv'
        print '\nWriting to: '+outfile
        infile=open(infile,'r')
        outfile=open(outfile,'w')
        count_total=0
        
############################ loop over line in samfile ##############        
        for line in infile.readlines():
            
            text=line.split()
############################ find LCA from all lines with the same sequence identifier (field #1 in samfile) ##############   
            if text[0]!=prev_name and count_total!=0:
                lines2=[line2 for line2 in lines if line2.split()[1].split('|')[0] not in wrong_tax]

                if len(lines2)>0:
                    outlines.append(get_LCA_from_blast(lines2,options.threshold,limits))
                else:
                    outlines.append(lines[0].split()[0]+'\tNOMATCH_all_taxids_ignored\t\t\t\t\t\t\t\n')
                
                count_total=0
                #break
            
            if count_total==0:
    
                prev_name=text[0]
                prev_text=text
                lines=[]
                lines2=[]
            lines.append(line)
            count_total+=1
        if 'lines' in locals():
            
            lines2=[line2 for line2 in lines if line2.split()[1].split('|')[0] not in wrong_tax]
            outlines.append(get_LCA_from_blast(lines2,options.threshold,limits))

############################ Sort output lines and print to file ############            
        [outfile.write(i) for i in smartsort(outlines)]
        outfile.close()
           
        
    
if __name__ == '__main__':
    main() 


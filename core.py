import msprime
import tskit
import random
import time
import numpy as np

#Function from https://tskit.dev/tutorials/introgression.html
def get_migrating_tracts(ts):
    # You should look at how the data is organized in tskit, it uses tables, as described here:
    # https://tskit.dev/tutorials/tables_and_editing.html
    # The following line is parsing through to population table in order to find the id of population ND.
    neanderthal_id = [p.id for p in ts.populations() if p.metadata['name']=='ND'][0]
    migrating_tracts = []
    # Get all tracts that migrated into the neanderthal population
    for migration in ts.migrations():
        if migration.dest == neanderthal_id:
            migrating_tracts.append((migration.left, migration.right))
    return np.array(migrating_tracts)


# This function will cut the genomic sequence in segments of length 'cut' and for each segment will count the number of mutations that individual 'ind' has, but which cannot be found in any individual belonging to population 'pop'.
# p1 is the position of first individual belonging to 'pop' and p2 the number of individual belonging to 'pop'. It is only there for the function to run faster, it can probably done in a cleaner one. TLDR: It is not important.

def createSeqObs(ts,cut,ind,p1,p2,pop) -> list:
    tables = ts.dump_tables()
    nodes = tables.nodes
    seq = np.zeros(int(ts.sequence_length/cut),dtype=int)  #The list which will contain the result of our function.
    pop_id = [p.id for p in ts.populations() if p.metadata['name']==pop][0] #The id of the population we are comparing our individual to.
    r=0 
    for v in ts.variants(): #This will loop through all the variants in the genome. Meaning all the positions where a mutation occured. 
        i=int(v.site.position/cut) # We want to cut our genome in small segments of length 'cut'. So we look at the position of the current variant and divide by 'cut'
        r=r+1
        b = 0
        c = False
        x=0
        while x < len(v.genotypes): # v.genotypes is the array containing the genotype of all the sampled individuals for the variant 'v'
            if(nodes[x].individual==ind): # Check if this position corresponds to individual 'ind'
                if(v.genotypes[x]==1): #If individual 'ind' has mutation '1' at this position, we put it in 'b', otherwise 'b' has value 0. At a site with a mutation, the genotype of an individual can take the value 0 or 1. In real life it will take a value on the alphabet A C G T, but this information is not important. 
                    b = 1
                x=p1
            elif(nodes[x].population==pop_id): # Check if this position corresponds to an individual in 'pop'
                if(v.genotypes[x]==b): 
                    c = True # The mutation found in 'ind' is also present in the population 'pop' 
                    x=p1+p2 # We can skip all the remaining individuals
    
            x=x+1
        if not c : # If c is False it means that the mutation in 'ind' has not been found in the population 'pop' 
            seq[i]=seq[i]+1
    return seq

# Mainly sorts the tracts in the must inneficient way, not important.
def clean_tracts(tractInit):
    tract = np.copy(tractInit)
    tract=tract.astype(int)
    flag = True
    while(flag):
        flag=False
        for i in range(len(tract)):
            for j in range(len(tract)):
                if not flag and tract[i,0]==tract[j,1]:
                    tract[j,1]=tract[i,1]
                    tract = np.delete(tract,i,0)
                    flag=True
    flag = True
    while(flag):
        flag=False
        for i in range(len(tract)):
            for j in range(i+1,len(tract)):
                if tract[i,0]>tract[j,0]:
                    save0=tract[i,0]
                    save1=tract[i,1]
                    tract[i,0]=tract[j,0]
                    tract[i,1]=tract[j,1]
                    tract[j,0]=save0
                    tract[j,1]=save1
                    flag=True
    return tract
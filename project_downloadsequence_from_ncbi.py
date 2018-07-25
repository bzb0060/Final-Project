##use BioPython package to download the required gene from NCBI/PubMed.
#Bio.Entrez visit http://eutils.ncbi.nlm.nih.gov
## Bio.Blast to blast these genes with reference gene


from Bio import Entrez #NCBI entrez database
Entrez.email = 'bzb0060@auburn.edu'  ##tell NCBI who you are
#handle = Entrez.einfo()
#result = handle.read()
#print handle
#print result

handle = Entrez.esearch(db='nucleotide', term = "Poaceae[Orgn] AND als[Gene]", retmax = 50) ##esearch(db, term, **keywds) #Run an Entrez search and return a handle to the results.
record = Entrez.read(handle)
handle.close() #make sure to close the file... its just an important convention to prevent you from doing anything you don't want to do
print record['Count']  #97
print len(record["IdList"]) #50
print record["IdList"] #['1376185352', '1376185350', '1376185348', '1376185346', '1269612522', '1269612520', '1269612518', ...]

                        
                                                                        
handle = Entrez.efetch(db="nucleotide", id="1376185352", rettype="gb", retmode="text") ##rettype and retmode to make the download data type is GenBank
#record = SeqIO.read(handle,"genbank")
print handle.read()
#print record
handle.close()


########## download sequences from NCBI, output to a new file 
def get_sequences(EntrezRecord,FileName):
    ids = EntrezRecord["IdList"]
    f = open('%s.fasta' % FileName, 'a')
    for seq_id in ids:
         handle = Entrez.efetch(db="nucleotide", id=seq_id, rettype="fasta", retmode="text")
         f.write(handle.read().strip())#need to create the file before you can do anything to it
         f.write("\n")
         handle.close()
         
    f.close() 
        
get_sequences(record, r'als1.fasta')


###tutorials: http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc120


####################################BLAST
from Bio.Blast import NCBIWWW

fasta_string = open("ALS1.fasta.txt").read()
result_handle = NCBIWWW.qblast("blastn", "nt", fasta_string) ##qbblast function

save_file = open("my_blast.xml", "w")
save_file.write(result_handle.read())
save_file.close()
result_handle.close()


from Bio.Blast import NCBIXML
def getblast(filename):
    handle = open(filename)
    blast_record = NCBIXML.read(handle)
    E_VALUE_THRESH = 0.0001
    for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                if hsp.expect <  E_VALUE_THRESH : 
                    print '**************Alignment***********'
                    print 'sequence:', alignment.title
                    print 'length:', alignment.length
                    print 'score:', hsp.score
                    print 'gaps:', hsp.gaps
                    print 'e-value:', hsp.expect
                    print hsp.query[0:100] +'...'
                    print hsp.match[0:100] +'...'
                    print hsp.sbjct[0:100] +'...'
getblast(r'my_blast.xml')
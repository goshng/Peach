# Use the following script to create
# smutans.gene2go and smutans.go2ngene
# gene_association.goa_uniprot
# gene_ontology.1_2.obo
# uniref90.fasta.gz
# ITAG2.3_proteins.fasta
# ITAG2.3_gene_models.gff3.txt.no-negative

export GOAUNIPROT=$(basename /Volumes/Elements/Documents/Projects/RNASeq-Analysis/data/goseq/gene_association.goa_uniprot.gz)
export GOBO=$(basename /Volumes/Elements/Documents/Projects/RNASeq-Analysis/data/goseq/gene_ontology.1_2.obo)
export UNIREFFASTA=$(basename /Volumes/Elements/Documents/Projects/RNASeq-Analysis/data/goseq/uniref90.fasta.gz)
export GOSEQFASTA=$(basename /Volumes/Elements/Documents/Projects/RNASeq-Analysis/data/analysis/ITAG2.3_proteins.fasta)
export REF1GFF=$(basename /Volumes/Elements/Documents/Projects/RNASeq-Analysis/data/analysis/ITAG2.3_gene_models.gff3.txt.no-negative)

# cp /v4scratch/sc2265/rnaseq/output/tomato/1/data/uniref90.fasta.gz output/tomato/1/data
gunzip output/tomato/1/data/uniref90.fasta.gz
UNIREFFASTA=${UNIREFFASTA%.gz}

cp /v4scratch/sc2265/rnaseq/output/tomato/1/data/$GOAUNIPROT output/tomato/1/data
gunzip output/tomato/1/data/$GOAUNIPROT

# 3. Make BLASTDB fo uniref90.fasta
./makeblastdb -in output/tomato/1/data/$UNIREFFASTA \
  -dbtype prot -title uniref90 -input_type fasta \
  -out output/tomato/1/data/uniref90

# 4. BLAST protein sequences
cp /v4scratch/sc2265/rnaseq/output/tomato/1/data/$GOSEQFASTA output/tomato/1/data
./blastp -db output/tomato/1/data/uniref90 -query output/tomato/1/data/$GOSEQFASTA \
  -task blastp \
  -outfmt 6 \
  -num_threads 8 \
  -evalue 1e-3 \
  -out output/tomato/1/run-analysis/goseq.blast
cp output/tomato/1/run-analysis/goseq.blast /v4scratch/sc2265/rnaseq/output/tomato/1/run-analysis

# 5. Create a file that associates genes with gene ontology terms
cp /v4scratch/sc2265/rnaseq/output/tomato/1/data/$REF1GFF output/tomato/1/data
perl pl/geneontology.pl gene2go \
  -gff output/tomato/1/data/$REF1GFF \
  -blast output/tomato/1/run-analysis/geneseq.blast \
  -goa output/tomato/1/run-analysis/gene_association.goa_uniprot \
  -out output/tomato/1/run-analysis/smutans.gene2go
cp output/tomato/1/run-analysis/smutans.gene2go /v4scratch/sc2265/rnaseq/output/tomato/1/run-analysis

# 6. Create a file that associates gene ontology terms with descriptions
perl pl/geneontology.pl go2ngene \
  -gene2go output/tomato/1/run-analysis/smutans.gene2go \
  -obo output/tomato/1/run-analysis/gene_ontology.1_2.obo \
  -out output/tomato/1/run-analysis/smutans.go2ngene
cp output/tomato/1/run-analysis/smutans.go2ngene /v4scratch/sc2265/rnaseq/output/tomato/1/run-analysis

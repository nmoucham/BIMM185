from Bio import SeqIO
import gzip
import os

num_genomes = 0

genbank_file = gzip.open('GCF_000005845.2_ASM584v2_genomic.gbff.gz',mode='rt')

num_genomes = num_genomes + 1
num_replicons = 0
global_CDS = 0
local_CDS = 0
num_genes_in_rep = 0
replicon_structure = ''
exon = 0

genes_table = open('genes_table.txt','w')
genomes_table = open('genomes_table.txt','w')
replicons_table = open('replicons_table.txt','w')
func_table = open('functions_table.txt','w')
syn_table = open('synonyms_table.txt','w')
extref_table = open('extref_table.txt','w')
exons_table = open('exons_table.txt','w')

for record in SeqIO.parse(genbank_file,'genbank'):
	genome_name = ' '.join(str(record.description).split()[0:2])
	size_bp = len(record.seq)
	domain = record.annotations['taxonomy'][0]
	date = record.annotations['date']
	accession = record.annotations['accessions'][0]
	shape = record.annotations['topology']
	if 'complete genome' in record.description:
		type = 'chromosome'

	assembly = ''.join(filter(lambda x:'Assembly' in x, record.dbxrefs)).split(':')[1]
	num_replicons = num_replicons + 1
	for feature in record.features:
		if feature.type == 'source':
			taxon = str(feature.qualifiers.get('db_xref'))
			tax_id = taxon.split("']")[0].split(':')[1]
			
		if feature.type == 'CDS':
			num_genes_in_rep = num_genes_in_rep + 1
			global_CDS = global_CDS + 1
			local_CDS = local_CDS + 1
			exon = 0
			
			locus_tag = feature.qualifiers.get('locus_tag')
			if locus_tag ==None:
				locus_tag = '-'
			
			protein_id = feature.qualifiers.get('protein_id')
			if protein_id == None:
				protein_id = '-'

			gene_name = feature.qualifiers.get('gene')
			if gene_name == None:
				gene_name = '-'
			
			coordinates = str(feature.location)
			direction = coordinates.split(']')[1]
			if direction == '(+)':
				strand = 'F'
			if direction == '(-)':
				strand = 'R' 

			start = str(feature.location).split('[')[1].split(':')[0]
			stop = str(feature.location).split('[')[1].split(':')[1].split(']')[0]
			start = int(start)
			stop = int(stop)
			length = -start + stop
			
			product = feature.qualifiers.get('product')
			if product == None:
				product = '-'
			
			function = feature.qualifiers.get('function')
			if function == None:
				function = '-'
			
			synonyms = feature.qualifiers.get('gene_synonym')
			if synonyms == None:
				synonyms = '-'

			genes_table.write(str(global_CDS))
			genes_table.write('\t')
			genes_table.write(str(num_genomes))
			genes_table.write('\t')
			genes_table.write(str(num_replicons))
			genes_table.write('\t')
			genes_table.write(str(''.join(locus_tag)))
			genes_table.write('\t')
			genes_table.write(str(''.join(protein_id)).split('.')[0])
			genes_table.write('\t')
			genes_table.write(str(''.join(gene_name)))
			genes_table.write('\t')
			genes_table.write(str(strand))
			genes_table.write('\t')
			genes_table.write(str(length))
			genes_table.write('\t')
			genes_table.write(str(''.join(product)))
			genes_table.write('\n')
			
			func_table.write(str(global_CDS))
			func_table.write('\t')
			func_table.write(str(''.join(function)))
			func_table.write('\n')
		
			syn_table.write(str(global_CDS))
			syn_table.write('\t')
			syn_table.write(str(''.join(synonyms)))
			syn_table.write('\n')

			extref_info = feature.qualifiers.get('db_xref')
			if extref_info == None:
				xdb = '-'
				xid = '-'
				
				extref_table.write(str(global_CDS))
				extref_table.write('\t')
				extref_table.write(str(xdb))
				extref_table.write('\t')
				extref_table.write(str(xid))
				extref_table.write('\n')
			else:
				for ref in range(0,len(extref_info)):
					xdb = extref_info[ref].split(':')[0]
					xid = extref_info[ref].split(':')[1].split('.')[0]
				
					extref_table.write(str(global_CDS))
					extref_table.write('\t')
					extref_table.write(str(xdb))
					extref_table.write('\t')
					extref_table.write(str(xid))
					extref_table.write('\n')
			
			if 'join' in str(feature.location):
				for i in feature.location.parts:
					exon = exon + 1
					
					left_pos = str(feature.location).split('[')[1].split(':')[0]
					right_pos = str(feature.location).split(':')[1].split(']')[0]
					if '<' in left_pos:
						continue
					if '>' in right_pos:
						continue
					
					left_pos = int(left_pos)
					right_pos = int(right_pos)
					exon_len = -left_pos + right_pos
					print(exon)
					exons_table.write(str(global_CDS))
					exons_table.write('\t')
					exons_table.write(str(exon))
					exons_table.write('\t')
					exons_table.write(str(left_pos))
					exons_table.write('\t')
					exons_table.write(str(right_pos))
					exons_table.write('\t')
					exons_table.write(str(exon_len))
					exons_table.write('\n')
			else:
				exon = 1
				left_pos = str(feature.location).split('[')[1].split(':')[0]
				right_pos = str(feature.location).split(':')[1].split(']')[0]
				if '<' in left_pos:
					continue
				if '>' in right_pos:
					continue
				
				left_pos = int(left_pos)
				right_pos = int(right_pos)
				exon_len = -left_pos + right_pos
				
				exons_table.write(str(global_CDS))
				exons_table.write('\t')
				exons_table.write(str(exon))
				exons_table.write('\t')
				exons_table.write(str(left_pos))
				exons_table.write('\t')
				exons_table.write(str(right_pos))
				exons_table.write('\t')
				exons_table.write(str(exon_len))
				exons_table.write('\n')

	
	replicons_table.write(str(num_replicons))
	replicons_table.write('\t')
	replicons_table.write(str(num_genomes))
	replicons_table.write('\t')
	replicons_table.write(str(genome_name))
	replicons_table.write('\t')
	replicons_table.write(str(type))
	replicons_table.write('\t')
	replicons_table.write(str(shape))
	replicons_table.write('\t')
	replicons_table.write(str(num_genes_in_rep))
	replicons_table.write('\t')
	replicons_table.write(str(size_bp))
	replicons_table.write('\t')
	replicons_table.write(str(accession))
	replicons_table.write('\t')
	replicons_table.write(str(date))
	replicons_table.write('\n')
	num_genes_in_rep = 0

genomes_table.write(str(num_genomes))
genomes_table.write('\t')
genomes_table.write(str(genome_name))
genomes_table.write('\t')
genomes_table.write(str(tax_id))
genomes_table.write('\t')
genomes_table.write(str(domain))
genomes_table.write('\t')
genomes_table.write(str(num_replicons))
genomes_table.write('\t')
genomes_table.write(str(local_CDS))
genomes_table.write('\t')
genomes_table.write(str(size_bp))
genomes_table.write('\t')
genomes_table.write(str(assembly))
genomes_table.write('\n')



os.chdir('/Users/nicole/bimm185/atumafaciens')

genbank_file2 = gzip.open('GCF_000576515.1_ASM57651v1_genomic.gbff.gz',mode = 'rt')
num_genomes = num_genomes + 1
genome_bp = 0
local_CDS = 0
num_genes_in_rep = 0

for record in SeqIO.parse(genbank_file2,'genbank'):
	genome_name = ' '.join(str(record.description).split()[0:2])
	replicon_size = len(record.seq)
	genome_bp = genome_bp + len(record.seq)
	domain = record.annotations['taxonomy'][0]
	date = record.annotations['date']
	accession = record.annotations['accessions'][0]
	if 'plasmid' in record.description:
		type = 'plasmid'
	if 'chromosome' in record.description:
		type = 'chromosome'
	shape = record.annotations['topology']
	assembly = ''.join(filter(lambda x: 'Assembly' in x, record.dbxrefs)).split(':')[1]
	num_replicons = num_replicons + 1

	for feature in record.features:
		if feature.type == 'source':
			taxon = str(feature.qualifiers.get('db_xref'))
			tax_id = taxon.split("']")[0].split(':')[1]
		
		if feature.type == 'CDS':
			local_CDS = local_CDS + 1
			global_CDS = global_CDS + 1
			num_genes_in_rep = num_genes_in_rep + 1
			exon = 0
				
			locus_tag = feature.qualifiers.get('locus_tag')
			if locus_tag == None:
				locus_tag = '-'

			protein_id = feature.qualifiers.get('protein_id')
			if protein_id == None:
				protein_id = '-'

			gene_name = feature.qualifiers.get('gene')
			if gene_name == None:
				gene_name = '-'
			
			direction = str(feature.location).split(']')[1]
			if direction == '(+)':
				strand = 'F'
			if direction == '(-)':
				strand = 'R'

			start = str(feature.location).split('[')[1].split(':')[0]
			stop = str(feature.location).split('[')[1].split(':')[1].split(']')[0]
			if '<' in start:
				start = start.replace('<','')
			if '>' in stop:
				stop = stop.replace('>','')
			start = int(start)
			stop = int(stop)
			length = -start + stop

			product = feature.qualifiers.get('product')
			if product == None:
				product = '-'

			function = feature.qualifiers.get('function')
			if function == None:
				function = '-'
			
			synonyms = feature.qualifiers.get('gene_synonym')
			if synonyms == None:
				synonyms = '-'
		
			genes_table.write(str(global_CDS))
			genes_table.write('\t')
			genes_table.write(str(num_genomes))
			genes_table.write('\t')
			genes_table.write(str(num_replicons))
			genes_table.write('\t')
			genes_table.write(str(''.join(locus_tag)))
			genes_table.write('\t')
			genes_table.write(str(''.join(protein_id)).split('.')[0])
			genes_table.write('\t')
			genes_table.write(str(''.join(gene_name)))
			genes_table.write('\t')
			genes_table.write(str(strand))
			genes_table.write('\t')
			genes_table.write(str(length))
			genes_table.write('\t')
			genes_table.write(str(''.join(product)))
			genes_table.write('\n')

			func_table.write(str(global_CDS))
			func_table.write('\t')
			func_table.write(str(''.join(function)))
			func_table.write('\n')
			
			syn_table.write(str(global_CDS))
			syn_table.write('\t')
			syn_table.write(str(''.join(synonyms)))
			syn_table.write('\n')

			extref_info = feature.qualifiers.get('db_xref')
			if extref_info == None:
				xdb = '-'
				xid = '-'
				
				extref_table.write(str(global_CDS))
				extref_table.write('\t')
				extref_table.write(str(xdb))
				extref_table.write('\t')
				extref_table.write(str(xid))
				extref_table.write('\n')
			else:
				for ref in range(0,len(extref_info)):
					xdb = extref_info[ref].split(':')[0]
					xid = extref_info[ref].split(':')[1].split('.')[0]
			

					extref_table.write(str(global_CDS))
					extref_table.write('\t')
					extref_table.write(str(xdb))
					extref_table.write('\t')
					extref_table.write(str(xid))
					extref_table.write('\n')

			if 'join' in str(feature.location):
				for i in feature.location.parts:
					exon = exon + 1
					left_pos = str(feature.location).split('[')[1].split(':')[0]
					right_pos = str(feature.location).split(':')[1].split(']')[0]
					if '<' in left_pos:
						continue
					if '>' in right_pos:
						continue
					
					left_pos = int(left_pos)
					right_pos = int(right_pos)
					exon_len = -left_pos + right_pos
					print(exon)		
					exons_table.write(str(global_CDS))
					exons_table.write('\t')
					exons_table.write(str(exon))
					exons_table.write('\t')
					exons_table.write(str(left_pos))
					exons_table.write('\t')
					exons_table.write(str(right_pos))
					exons_table.write('\t')
					exons_table.write(str(exon_len))
					exons_table.write('\n')
			else:
				exon = 1
				left_pos = str(feature.location).split('[')[1].split(':')[0]
				right_pos = str(feature.location).split(':')[1].split(']')[0]
				if '<' in left_pos:
					continue
				if '>' in right_pos:
					continue
				
				left_pos = int(left_pos)
				right_pos = int(right_pos)
				exon_len = -left_pos + right_pos
				
				exons_table.write(str(global_CDS))
				exons_table.write('\t')
				exons_table.write(str(exon))
				exons_table.write('\t')
				exons_table.write(str(left_pos))
				exons_table.write('\t')
				exons_table.write(str(right_pos))
				exons_table.write('\t')
				exons_table.write(str(exon_len))
				exons_table.write('\n')


	replicons_table.write(str(num_replicons))
	replicons_table.write('\t')
	replicons_table.write(str(num_genomes))
	replicons_table.write('\t')
	replicons_table.write(str(genome_name))
	replicons_table.write('\t')
	replicons_table.write(str(type))
	replicons_table.write('\t')
	replicons_table.write(str(shape))
	replicons_table.write('\t')
	replicons_table.write(str(num_genes_in_rep))
	replicons_table.write('\t')
	replicons_table.write(str(replicon_size))
	replicons_table.write('\t')
	replicons_table.write(str(accession))
	replicons_table.write('\t')
	replicons_table.write(str(date))
	replicons_table.write('\n')
	num_genes_in_rep = 0

genomes_table.write(str(num_genomes))
genomes_table.write('\t')
genomes_table.write(str(genome_name))
genomes_table.write('\t')
genomes_table.write(str(tax_id))
genomes_table.write('\t')
genomes_table.write(str(domain))
genomes_table.write('\t')
genomes_table.write(str(num_replicons))
genomes_table.write('\t')
genomes_table.write(str(local_CDS))
genomes_table.write('\t')
genomes_table.write(str(genome_bp))
genomes_table.write('\t')
genomes_table.write(str(assembly))
genomes_table.write('\n')

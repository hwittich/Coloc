import os
import os.path
import argparse

'''
Set up argparser
'''
parser = argparse.ArgumentParser(description = 'Run GWAS-eQTL colocalization pipeline.')
#Add arguments
parser.add_argument('--gwas', required=True, help = 'GWAS summary statistics directory path')
parser.add_argument('--eqtl', required=True, help='Exact eQTL txt file path')
parser.add_argument('--vcf', required=True, help = 'Directory containing the vcf files from the eQTL population')
parser.add_argument('--snp_annot', required=True, help = 'SNP annotation files directory path')
parser.add_argument('--gene_annot', required=True, help = 'Exact gene annotation file path')
parser.add_argument('--geno', required=True, help = 'Genotype files directory path')
parser.add_argument('--frq', required = True, help = 'Exact frq file path')
parser.add_argument('--out', required=True, help = 'Specify main output directory for all output file')
parser.add_argument('--pop1', required=True, help = 'Population used for vcf files')
parser.add_argument('--pop4', required=True, help = 'Population used for frq and eQTL files')
parser.add_argument('--pop_size', type = int, required=True, help = 'Size of population')
parser.add_argument('--phenotypes', nargs = '*', required=True, help = 'Phenotypes to test')
parser.add_argument('--chrs', type = int, nargs= '*', required=True, help = 'Indicate what chromosomes to run')
#Flag to specify whether to run all genes or genes specified by user
parser.add_argument('-gene_id', type = str, nargs='*', default= False, help = 'Input specific genes ids to run colocalization on, default is set to False')
args = parser.parse_args()

'''
Store arguments
'''
gwas = args.gwas
eqtl = args.eqtl
vcf = args.vcf
if vcf[-1] == "/":
    vcf = vcf[:-1]
snp_annot = args.snp_annot
gene_annot = args.gene_annot
geno = args.geno
frq = args.frq
out = args.out
pop1 = args.pop1
pop4 = args.pop4
pop_size = args.pop_size
phenos = args.phenotypes
chrs = args.chrs
#Convert lists from int to str
chrs = list(map(str, chrs))
pop_size = str(pop_size)
#Get Coloc repo working directory
wd = os.getcwd()

'''
Check/create necessary output directories
''' 
if os.path.isdir(out) == False:
    os.mkdir(out)
path = args.out + '/'
#Make main subdirectories
dirs = ['LD_matrix', 'gene_lists', 'GWAS_TOPMED', 'eQTL', 'Coloc_output']
for dir in dirs:
    temp_path = path + dir
    if os.path.isdir(temp_path) == False:
        os.mkdir(temp_path)
#Make population directories in subdirectories 
temp_path = path + 'LD_matrix/' + pop1
if os.path.isdir(temp_path) == False:
    os.mkdir(temp_path)
temp_path = temp_path + '/' + pop1 + '_1Mb_coords_LDMatrix'
if os.path.isdir(temp_path) == False:
    os.mkdir(temp_path)
temp_path = path + 'gene_lists/' + pop1
if os.path.isdir(temp_path) == False:
    os.mkdir(temp_path)
#Make population directories in subdirectories 
temp_path = path + 'GWAS_TOPMED/' + pop4
if os.path.isdir(temp_path) == False:
    os.mkdir(temp_path)
temp_path = path + 'eQTL/' + pop4
if os.path.isdir(temp_path) == False:
    os.mkdir(temp_path)
#Save LD ouput dir
ld = path + 'LD_matrix/'

print('Directories created.')

'''
Run Scripts
'''
print('Run ' + pop1 + '.', flush=True)

## Script 2
for chr in chrs:
    #Get vcf file specific to chr
    vcf_files = os.listdir(vcf)
    for file in vcf_files:
        if chr+"." in file:
            chr_vcf = vcf+"/"+file
    #maybe implement subprocess for parallelization
    os.system("chmod u+x 02_make_bed.sh") #Make the script executable
    cmd = "./02_make_bed.sh "+pop1+" "+chr+" "+chr_vcf+" "+out
    os.system(cmd)
print('Bfiles made.')

## Scripts 1
#Get files in dirs
geno_files = os.listdir(geno)
snp_annot_files = os.listdir(snp_annot)

chr_snp=""
chr_geno=""

#Get file specific to chr
for chr in chrs:
    print(chr,flush=True)
    for file in geno_files:
        print(file,flush=True)
        if "chr"+chr+"." in file:
            print("Found genotype file",flush=True)
            print(file,flush=True)
            chr_geno = geno+"/"+file
    for file in snp_annot_files:
        print(file,flush=True)
        if "chr"+chr+"." in file:
            print("Found SNP annotation file",flush=True)
            print(file,flush=True)
            chr_snp = snp_annot+"/"+file
    #Run command
    script1cmd = 'Rscript 01b_run_pull_snps_driving.R ' + chr + ' ' + chr_snp + ' ' + gene_annot + ' ' + chr_geno + ' ' + out + ' ' + pop1
    print(script1cmd)
    os.system(script1cmd)

print('Pulling SNPs completed.')

## Script 3
for chr in chrs:
    os.system("chmod u+x 03_make_LD_matrix.sh")
    cmd = "./03_make_LD_matrix.sh "+pop1+" "+chr+" "+out
    os.system("./03_make_LD_matrix.sh "+pop1+" "+chr+" "+out)

print('LD matrices created.')

## Run scripts w/specified gene list
if args.gene_id != False:
    
    #Get genes
    gene_ids = args.gene_id
    
    for pheno in phenos:
        
        ## Script 4

        #Get list of files in gwas directory
        gwas_files = os.listdir(gwas)
        #Get gwas file specific to phenotype
        for file in gwas_files:
            if pheno in file:
                pheno_gwas_file = file
        #Run script 4 command
        chrs_unlist = ' '.join(chrs)
        cmd = 'Rscript 04_prep_files_coloc.R '+ gwas + '/' + pheno_gwas_file + ' ' + eqtl + ' ' + frq + ' ' + out \
                    + ' ' + pop4 + ' ' + pop_size + ' ' + pheno + ' ' + chrs_unlist
        os.system(cmd)

        print('Input files formatted.')


        ## Scripts 5
        
        #Get formatted gwas file specific to pop and phenotype
        out_gwas = out + '/GWAS_TOPMED/' + pop4 
        gwas_files = os.listdir(out_gwas)
        for file in gwas_files:
            if pheno in file:
                if pop4 in file: 
                    pheno_pop_gwas = file
        #Get formatted eqtl file specific to pop & phenotype
        out_eqtl = out + '/eQTL/' + pop4
        eqtl_files = os.listdir(out_eqtl)
        for file in eqtl_files:
            if pheno in file:
                if pop4 in file:
                    pheno_pop_eqtl = file
        #Get LD dir for this pop
        pop_ld = ld + pop1 + '/' + pop1 + '_1Mb_coords_LDMatrix'
        #Convert list of genes to string
        genes_unlist = ' '.join(gene_ids)
        #Run script 5 command
        cmd = 'Rscript 05b_run_coloc.R ' + pheno_pop_gwas + ' ' + pheno_pop_eqtl + ' ' + pop_ld + ' ' + out \
                    + ' ' + pop1 + ' ' + pop_size + ' ' + pheno + ' ' + genes_unlist
        os.system(cmd)

        print('Coloc analysis finished')

## Run all genes in chromosomes
else:
    for pheno in phenos:

        ## Script 4
       
        #Get list of files in gwas directory
        gwas_files = os.listdir(gwas)
        #Get gwas file specific to phenotype
        for file in gwas_files:
            if pheno in file:
                pheno_gwas_file = file
        #Run script 4 command
        chrs_unlist = ' '.join(chrs)
        cmd = 'Rscript 04_prep_files_coloc.R '+ gwas + '/' + pheno_gwas_file + ' ' + eqtl + ' ' + frq + ' ' + out \
                    + ' ' + pop4 + ' ' + pop_size + ' ' + pheno + ' ' + chrs_unlist
        os.system(cmd)

        print('Input files formatted.')


        ## Scripts 5
        
        #Get formatted gwas file specific to pop and phenotype
        out_gwas = out + '/GWAS_TOPMED/' + pop4 
        gwas_files = os.listdir(out_gwas)
        for file in gwas_files:
            print(file,flush=True)
            if pheno in file:
                print("Found phenotype: "+pheno,flush=True)
                if pop4 in file: 
                    print("Found population: "+pop4,flush=True)
                    print("Found file",flush=True)
                    print(file,flush=True)
                    pheno_pop_gwas = file
        #Get formatted eqtl file specific to pop & phenotype
        out_eqtl = out + '/eQTL/' + pop4
        eqtl_files = os.listdir(out_eqtl)
        for file in eqtl_files:
            if pheno in file:
                if pop4 in file:
                    pheno_pop_eqtl = file
        #Get LD dir for this pop
        pop_ld = ld + pop1 + '/' + pop1 + '_1Mb_coords_LDMatrix'
        #Run script 5 command
        cmd = 'Rscript 05b_run_coloc.R ' + pheno_pop_gwas + ' ' + pheno_pop_eqtl + ' ' + pop_ld + ' ' + out \
                    + ' ' + pop1 + ' ' + pop_size + ' ' + pheno
        os.system(cmd)

        print('Coloc analysis finished')


print('Pipeline completed running.')

import pandas as pd
import numpy as np
import os
from Bio import Phylo
from os.path import dirname
#from ete3 import Tree
import json
from orthofinder_utils import proteomes_for_orthofinder as pfo
from jw_utils import file_utils as fu
from jw_utils import parse_gff as pgf
import shutil

#########################################################
# !!!!CHANGE WORKING DIRECTORY TO dash_app main folder

# make a directory insite the main /dash_app folder called data (dash_app/data). Place the ncbi_dataset/data
# containing all the data downloaded from ncbi into the ./data folder () dash_app/data/ncbi_datasets/data
#
# you will have run orthofinder_utils/proteomes_for_orthofinder.py to generate the Proteomes directory
#
#You will have to write a summary.json file using

# after running orthofinder algorithm, orthofinder will have deposited the orthofinder results into the 
# Proteomes folder .../Proteomes/OrthoFinder/Results_date. Make sure this folder is in .../dash_app/data 

# if pushing app to heroku:
#       the orthofinder folder should be added to .gitignore in order to not push all that data
#       to their servers
#
############################################################################
#
# =============================================================================
# Enter variables here:

####Run this function to set up files and directories
 



def run_all(path_to_results):
    path_to_app_folder = '/'.join(path_to_results.split('/')[:-4])
    check_dir_structure(path_to_app_folder, path_to_results)
    write_AssemblyAccession_to_SpeciesName_from_summaries(path_to_results)
    copy_folders_from_orthoresults(path_to_results)
    os.chdir(path_to_app_folder)
    initial_check(path_to_results)
    #make_genome_annotation_files(path_to_results)
    
#
# =============================================================================


def check_dir_structure(path_to_parent, path_to_results):
    if not os.path.exists('./data'):
        raise Exception('There is no "./data" dir in the parent folder')
    if not os.path.exists('./data/ncbi_dataset/data'):
        raise Exception('There is no "./data/ncbi_dataset/data" dir in the parent folder')
    if not os.path.exists('./data/Proteomes/OrthoFinder'):
        raise Exception('There is no "./data/Proteomes/OrthoFinder" dir in the parent folder')
    if not os.path.exists('./data/summary_data'):
        os.makedirs('./data/summary_data')
    if not os.path.exists('./data/summary_data/summaries.json'):
        write_accessions_to_file(path_to_results)
        m = 'There is no "./data/summary_data/summaries.json" dir in the parent folder\n'
        m = m+f'in the terminal:\n'
        m = m+ f'$ cd .../dash_app\n'
        m = m+'$ datasets summary genome accession --inputfile ./data/summary_data/accessions.txt > ./data/summary_data/summaries.json'
        raise Exception(m)

    
def write_accessions_to_file(path_to_results):
    accessions = accessions_from_orthofinder_files(path_to_results)
    with open("./data/summary_data/accessions.txt", 'w') as f:
        for acc in accessions:
            f.write(acc+'\n')  



def accessions_from_proteome_dir(proteomes_dir):
    accessions = os.listdir(proteomes_dir)
    return [acc[:-4] for acc in accessions if acc.endswith('.faa')]


def accessions_from_orthofinder_files(path_to_results):
    return list(pd.read_csv(os.path.join(path_to_results, 
                    'Orthogroups/Orthogroups.tsv'), sep = '\t').columns[1:])


def initial_check(path_to_results):
    """
    check that assemblies in orthogroups.tsv, proteome folder, and ncbi/data folder are equivalent
    """
    orth_assemblies = list(pd.read_csv(os.path.join(path_to_results, 
                        'Orthogroups/Orthogroups.tsv'), sep = '\t').columns[1:]) 
    proteome_assemblies = accessions_from_proteome_dir('./data/Proteomes')
    
    ncbidata_assemblies = [folder for folder in os.listdir('./data/ncbi_dataset/data')
                   if folder.startswith('GC')]
    orth_assemblies.sort()
    proteome_assemblies.sort()
    ncbidata_assemblies.sort()
    if not ncbidata_assemblies == proteome_assemblies == orth_assemblies:
        print(
              'There are discrpencies in the assembly data!!\n'
              f'assemblies present in ./data/ncbi_dataset/data/: {len(ncbidata_assemblies)}\n'
              f'assemblies present in ./data/Proteomes: {len(proteome_assemblies)}\n'
              f'assemblies present in\n'
              f'{os.path.join(path_to_results, "Orthogroups/Orthogroups.tsv")}: {len(orth_assemblies)}'
              )
    else:
        print('all checks produce equivalent assembly identifiers')
        return orth_assemblies
        
        
def make_dict_fromNCBI_summaries(path_to_json):
    with open(path_to_json, 'r') as f:
        wierd_nested_d = json.load(f)['assemblies']
    summaries_d = {}
    for ass in wierd_nested_d:
        summaries_d[ass['assembly']['assembly_accession']] = ass['assembly']
    return summaries_d
        

def write_AssemblyAccession_to_SpeciesName_from_summaries(path_to_results, output_file_path=None):
    path_to_json = './data/summary_data/summaries.json'
    summary_dict = make_dict_fromNCBI_summaries(path_to_json)
    acc2name = {}
    for acc, summ, in summary_dict.items():
        if len(summ['org']['sci_name'].split(' '))>2:
            acc2name[acc] = summ['org']['sci_name']
        elif summ['org'].get('strain'):
            acc2name[acc] = summ['org']['sci_name'] +' '+ summ['org'].get('strain')
        elif summ['org'].get('isolate'):
            acc2name[acc] = summ['org']['sci_name'] +' '+ summ['org'].get('isolate')
        else:
            acc2name[acc] = summ['org']['sci_name']
    with open('./data/summary_data/AssemblyAccession_to_SpeciesName.json', 'w') as f:
        json.dump(acc2name, f)  
        
        
def write_AssemblyAccession_to_SpeciesName(path_to_results, output_file_path=None):
    'writes json file with {assemb_accession:name}'
    
    acc_2_name_dict, assembly_dict = pfo.get_summary_loop(
            accessions_from_orthofinder_files(path_to_results), write_file=False)
    if output_file_path:
        output_file_path = output_file_path
    else:
        output_file_path = './data/summary_data/AssemblyAccession_to_SpeciesName.json'
        if not os.path.isdir('./data/summary_data'):
            os.makedirs('./data/summary_data')
        if os.path.isfile(output_file_path):
            raise Exception('File already exists')
    with open(output_file_path, 'w') as f:
        json.dump(acc_2_name_dict, f)    
    return acc_2_name_dict



def make_genome_annotation_files(path_to_results):
    op = OrthogroupParser(path_to_results)
    files = os.listdir('./data/ncbi_dataset/data')
    assembly_accessions = [file for file in files if file.startswith('GC')]
    ortho_accessions = accessions_from_orthofinder_files(path_to_results)
    for assembly in ortho_accessions:
        if assembly not in assembly_accessions:
            raise Exception('Orthofinder proteome no "in ncbi_datasets/data"')
    #make genome annotation dir if it doesnt exist already
    genome_annotations = './data/genome_annotations'
    if not os.path.isdir(genome_annotations):
        os.makedirs(genome_annotations)
    else:
        ans = input(
                f'the directory {genome_annotations} already exists,\n'
                f'do you still want to write files to it?  y or n?'
                )
        if ans != "y":
            raise Exception('Cannot write into existing directory')
            
        
    for assembly_accession in ortho_accessions:
        path_to_gff3 = os.path.join('./data/ncbi_dataset/data/',assembly_accession, 'genomic.gff')
        seq_dict = pgf.make_seq_object_dict(path_to_gff3, feature_type='CDS')
        #make dict to turn in
        df_dict = {
                'IDs':[seq_dict[key].ID for key in seq_dict.keys()],
                'starts':[int(seq_dict[key].start) for key in seq_dict.keys()],
                'ends':[int(seq_dict[key].end) for key in seq_dict.keys()],
                'strand': [seq_dict[key].strand for key in seq_dict.keys()],
                'products':[seq_dict[key].product for key in seq_dict.keys()],
                'Parents':[seq_dict[key].Parent for key in seq_dict.keys()]
                }
        df = pd.DataFrame(df_dict)
        df = df.set_index('IDs')
        #add orthogroup and HOG column to df

        gene_to_orthogroup_dict = op.gene_to_orthogroup_dict(assembly_accession)
        gene_to_HOG_dict = op.gene_HOG_dict(assembly_accession)
        ortho_list = []
        HOG_list = []
        for protein in df.index:
            if protein.startswith('cds-'):
                protein = protein[4:]
            if gene_to_orthogroup_dict.get(protein):
                ortho_list.append(gene_to_orthogroup_dict[protein])
            else:
                ortho_list.append('no orthogroup')
            if gene_to_HOG_dict.get(protein):
                HOG_list.append(gene_to_HOG_dict.get(protein))
            else:
                HOG_list.append('no HOG')
        df['Orthogroups'] = ortho_list
        df['HOGs'] = HOG_list
        filename = assembly_accession +  '_annotations.csv'
        df.to_csv(f'./data/genome_annotations/{filename}')



def copy_folders_from_orthoresults(path_to_results):
    """Copy data from orthofinder/results (.gitignored) folder into ./Data folder """
    import shutil
    #src_ortho = os.path.join(path_to_results,'Orthogroups')
    #shutil.copytree(src_ortho, './data/Orthogroups')
    src_HOGs = os.path.join(path_to_results,'Phylogenetic_Hierarchical_Orthogroups')
    shutil.copytree(src_HOGs, './data/Phylogenetic_Hierarchical_Orthogroups')
    src_gene_tree = os.path.join(path_to_results,'Resolved_Gene_Trees')
    shutil.copytree(src_gene_tree, './data/Resolved_Gene_Trees')
    src_spec_tree = os.path.join(path_to_results,'Species_Tree')
    shutil.copytree(src_spec_tree, './data/Species_Tree')
    OrthogroupParser(path_to_results).N0_HOG_counts
    src_N0_counts = os.path.join(path_to_results,'N0_HOG_counts.tsv')
    shutil.copyfile(src_N0_counts, './Data/N0_HOG_counts.tsv')
        
    



class OrthogroupParser():
    def __init__(self,path_to_results, *args, **kwargs):
        
        self.path_to_results_folder = path_to_results
        self.path_to_dir = self.get_path_to_dir()
        self.path_to_orthogroups = os.path.join(self.path_to_results_folder, 'Orthogroups/Orthogroups.tsv')                          
        self.path_to_gene_counts = os.path.join(self.path_to_results_folder,'Orthogroups/Orthogroups.GeneCount.tsv')       
        self.accession_to_name = self.accession_to_name()
        self.name_to_accession = self.name_to_accession()
        self.genomes = list(self.accession_to_name.keys())
        self.gene_counts_per_orthogroup = self.get_gene_counts()
        self.orthogroups = list(self.gene_counts_per_orthogroup.index)
        self.N0_HOG_path = os.path.join(self.path_to_results_folder, 'Phylogenetic_Hierarchical_Orthogroups/N0.tsv') 
        self.HOG_dir_path = os.path.join(self.path_to_results_folder,'Phylogenetic_Hierarchical_Orthogroups')
        self.N0_HOG_counts = self.N0_HOG_counts_check()
        self.HOGs = list(self.N0_HOG_counts['HOG'])
        
    
    def all_prots_in_HOG(self, HOG):
        'return a list of all proteins in a given HOG e.g. GCF_000332095_2_WP_008307711.1'
        HOG_df = pd.read_csv(self.N0_HOG_path, sep = '\t', low_memory=False).set_index('HOG')
        all_HOG_prots = HOG_df.loc[HOG,HOG_df.columns[2:]]
        #append accession and gene name, e.g GCF_000332095_2_WP_008307711.1
        all_HOG_prots = all_HOG_prots[~all_HOG_prots.isnull()]
        HOG_prot_names = []
        for accession in all_HOG_prots.index:
            prots = all_HOG_prots[accession].split(', ')
            for prot in prots:
                HOG_prot_names.append(accession.replace('.','_') + '_' + prot.strip())
        return HOG_prot_names
                
    def all_prots_in_orthogroup(self, orthogroup):
        'return a list of all proteins in a given orthogroup e.g. GCF_000332095_2_WP_008307711.1'
        orth_df = pd.read_csv(self.path_to_orthogroups, sep = '\t').set_index('Orthogroup')
                
        all_orth_prots = orth_df.loc[orthogroup,:]
        all_orth_prots = all_orth_prots[~all_orth_prots.isnull()]
        orth_prot_names = []
        for accession in all_orth_prots.index:
            prots = all_orth_prots[accession].split(', ')
            for prot in prots:
                orth_prot_names.append(accession.replace('.','_') + '_' + prot.strip())
        return orth_prot_names
    
    def HOG_OG_dict(self):
        df = pd.read_csv(self.HOG_dir_path + '/N0.tsv', sep = '\t', usecols=['HOG', 'OG']).set_index('HOG')
        HOG_ortho_dict = {}
        for HOG in self.HOGs:
            HOG_ortho_dict[HOG] = df.loc[HOG,'OG']
        return HOG_ortho_dict
    
    
    def OG_HOG_dict(self):
        OG_HOG_dict = {}
        HOG_OG_dict = self.HOG_OG_dict()
        for HOG, OG in HOG_OG_dict.items():
            if OG_HOG_dict.get(OG):
                OG_HOG_dict[OG].append(HOG)
            else:
                OG_HOG_dict[OG] = [HOG]
        return OG_HOG_dict

    
    
    def N0_HOG_counts_check(self):
        path_to_counts = os.path.join(self.path_to_results_folder, 'N0_HOG_counts.tsv')
        try:
            df = pd.read_csv(path_to_counts)
            accessions = list(self.accession_to_name.keys())
            diff = set(accessions).difference(set(df.columns[1:]))
            if len(diff) >0:
                raise Exception('The species in "N0_HOG_counts.tsv" do not match the species in acccession to name dict')
        except: 
            df = self.N0_HOG_counts()
            df.to_csv(path_to_counts)
        return df
    
    
    def get_gene_counts(self):
        return pd.read_csv(self.path_to_gene_counts, sep='\t').set_index('Orthogroup')
        
    
    
    def get_path_to_dir(self):
        parent = dirname(self.path_to_results_folder)
        parent2 = dirname(parent)
        return dirname(parent2)
        
        
    def get_path_components(self):
        normalized_path = os.path.normpath(self.path_to_results_folder) 
        path_components = normalized_path.split(os.sep)
        return path_components

    
    def accession_to_name(self):
        'return a dictionary with accessions as keys and names as values'
    
        '''
        path_to_dir (str): path to directory containing 'ncbi_dataset' and 'summary' directories
    
        '''
                       
        path_to_json = './data/summary_data/AssemblyAccession_to_SpeciesName.json'
        with open(path_to_json) as f:
            return json.load(f)
    
    
    def name_to_accession(self):
    
        name2assembly = {}
        for key in self.accession_to_name.keys():
            name2assembly[self.accession_to_name[key]] = key 
    
        return name2assembly    
    
      
    def list_name_to_accession(self, list_of_names):
        
        return {name:self.name_to_accession[name] for name in list_of_names}
    
    
    
    def orthogroup_gene_dict(self, species_accession): 
       #make df with orhtogroup and all genes in a given organism
        df = pd.read_csv(self.path_to_orthogroups, sep = '\t', usecols = ["Orthogroup", species_accession], low_memory=False)
        df = df.astype(str)
        df = df.set_index('Orthogroup')
        orthogroup_gene_dict = {}
        for orthogroup in df.index:
            orthogroup_gene_dict[orthogroup] = df.loc[orthogroup,species_accession].split(',')
        return orthogroup_gene_dict
    
    
    def HOG_gene_dict(self, species_accession): 
       #make df with orhtogroup and all genes in a given organism
        df = pd.read_csv(self.N0_HOG_path, sep = '\t', usecols = ["HOG", species_accession], low_memory=False)
        df = df.astype(str)
        df = df.set_index('HOG')
        HOG_gene_dict = {}
        for HOG in df.index:
            HOG_gene_dict[HOG] = df.loc[HOG,species_accession].split(', ')
        return HOG_gene_dict
    
    def gene_HOG_dict(self, species_accession):
        gene_HOG_dict = {}
        HOG_gene_dict = self.HOG_gene_dict(species_accession)
        for HOG, proteins in HOG_gene_dict.items():
            proteins =  HOG_gene_dict[HOG]
            if proteins[0] != 'nan':
                for protein in proteins:
                    gene_HOG_dict[protein.strip()] = HOG
        return gene_HOG_dict
    
    
    
    def gene_to_orthogroup_dict(self, species_accession):
        'return a dictionary of all genes/proteins (key) to orthogroups (value)'
        g2o_dict = {}
        d = self.orthogroup_gene_dict(species_accession)
        for orthogroup, protein_lst in d.items(): 
            if len(protein_lst) >1:
                for protein in protein_lst:
                    protein = protein.strip()
                    g2o_dict[protein]=orthogroup
            else:
                g2o_dict[protein_lst[0].strip()]= orthogroup
        return g2o_dict
    
       
    def get_og_saturation(self, group, copies):
        'return % of genomes in group that contain each orthogroup'
        '''
        group (list):
        return (dict): orthogroup:%present (string:float)
        '''
        gene_counts_df = pd.read_csv(self.path_to_gene_counts, sep='\t').set_index('Orthogroup')
        df_grp = gene_counts_df[group]
        filt =  (df_grp >= copies)
        return ((filt.sum(axis=1)/len(group)*100)).to_dict()
    
    
    
    def get_HOG_saturation(self, group, copies):
        gene_counts_df = self.N0_HOG_counts.set_index('HOG')
        df_grp = gene_counts_df[group]
        filt =  (df_grp >= copies)
        return ((filt.sum(axis=1)/len(group)*100)).to_dict()
    
    
    def genomes_missing_HOGs(self,ing, outg, HOG, copies):
        'get ingroup and outgroup genome accessions with below threshold paralogs in orthogroup'
        gene_counts_df = self.N0_HOG_counts.set_index('HOG')
        df_ingrp = gene_counts_df.loc[HOG,ing]
        df_outgrp = gene_counts_df.loc[HOG,outg]
        if copies >= 2:
            outgroup_copies = 2
        else:
            outgroup_copies = copies
        filt_in =  (df_ingrp >= copies)
        filt_out =  (df_outgrp >= outgroup_copies)
        missing_dict = {}
        missing_dict['ingroup'] = (list(df_ingrp[~filt_in].index))
        missing_dict['outgroup'] = (list(df_outgrp[~filt_out].index))
        return missing_dict
        

        
    
    
    def og_enrichment(self, ingroup, outgroup, copies):
        'return % dif in presence of orthogroups in ingroup vs outgroup genomes'
        '''
        ingroup (list):
        outgroup (list):
        return (dict): orthogroup:%difference  (string:float)
        '''
        in_dict = self.get_og_saturation(ingroup, copies)
        
        if copies >= 2:
            outgroup_copies = 2
        else:
            outgroup_copies = copies
        out_dict = self.get_og_saturation(outgroup, outgroup_copies)
        orthogroups = self.orthogroups
        enrich_dict =  {orthgrp:(in_dict[orthgrp] - out_dict[orthgrp]) for orthgrp in orthogroups}
        return enrich_dict
    

    def HOG_enrichment(self, ingroup, outgroup, copies):
        'return % dif in presence of HOGs in ingroup vs outgroup genomes'
        '''
        ingroup (list):
        outgroup (list):
        return (dict): orthogroup:%difference  (string:float)
        '''
        
        in_dict = self.get_HOG_saturation(ingroup, copies)
        
        
        if copies >= 2:
            outgroup_copies = 2
        else:
            outgroup_copies = copies
        out_dict = self.get_HOG_saturation(outgroup, outgroup_copies)
        
        
        df1 = pd.DataFrame().from_dict(in_dict, orient='index')
        df2 = pd.DataFrame().from_dict(out_dict, orient='index')
        df1.columns = ['ingroup_saturation']
        df2.columns = ['Outgroup_saturation']
        df = df1.merge(df2, left_index = True, right_index= True)
        filt = (df['ingroup_saturation'] != 0) & (df['Outgroup_saturation'] !=0)
        df_filt = df.loc[filt,:]
        enrichment_df = df_filt['ingroup_saturation'] - df_filt['Outgroup_saturation']
        enrich_dict = enrichment_df.to_dict()
        
        return enrich_dict
    
    
    def get_common_orthogroups(self, ingroup, threshold = 100, copies=1):
        
        'get orthogroups common to ingroup at a user defined threshold'
        '''
        ingroup (list): 
        outgroup (list):
        threshold (int of float):
        copies (int):
        return (list):
        '''
        gene_counts_df = pd.read_csv(self.path_to_gene_counts, sep='\t').set_index('Orthogroup')
        ingroup_counts_df = gene_counts_df.loc[:,ingroup]
        OG_percent_ingroup = {}
        for orthogroup in ingroup_counts_df.index:
            og_count_ingroup=0
            for genome in ingroup:
                if ingroup_counts_df.loc[orthogroup,genome] >= copies:
                    og_count_ingroup +=1
            og_count_ingroup = (og_count_ingroup/len(ingroup))*100
            if og_count_ingroup >= threshold:
                OG_percent_ingroup[orthogroup] = og_count_ingroup
        return list(OG_percent_ingroup.keys())
    
    #Need to fix the path to orthogroups_df
  
    def N0_HOGs_proteins(self, tax_level = None):
        if not tax_level:
            tax_level = 'N0'
        path = os.path.join(self.HOG_dir_path, tax_level + '.tsv')
        HOGs_df = pd.read_csv(path, sep = '\t')
        HOGs_df = HOGs_df.set_index('HOG')
        HOG_proteins = {}
        
        for assembly_identifier in HOGs_df.columns[2:]:
            filt = HOGs_df[assembly_identifier].notnull()
            populated_OGs = HOGs_df.loc[filt, assembly_identifier]
            for HOG in populated_OGs.index:
                if HOG_proteins.get(assembly_identifier):
                    HOG_proteins[assembly_identifier].update({HOG:populated_OGs.loc[HOG].split(',')})
                
                else:
                    HOG_proteins[assembly_identifier] = {HOG:populated_OGs.loc[HOG].split(',')}
        return HOG_proteins
    
    
    
    def N0_HOG_counts(self, tax_level = None):
        HOG_proteins = self.N0_HOGs_proteins()
        lst_of_dfs = []
        for accession in HOG_proteins.keys():
            df = pd.DataFrame()
            counts = []
            HOGs = []
            for HOG, proteins in HOG_proteins[accession].items():
                counts.append(len(proteins))
                HOGs.append(HOG)
            df['HOG'] = HOGs
            df[accession] = counts
            #df = df.set_index('HOG')
            lst_of_dfs.append(df)
        HOG_counts_df = lst_of_dfs[0]
        for df in lst_of_dfs[1:]:
            HOG_counts_df = pd.merge(HOG_counts_df, df, left_on='HOG', right_on='HOG', how  = 'outer')
        return HOG_counts_df
        
 
    
    def get_unique_orthogroup_dups(self, input_accessions, comparision_accessions):
        #get total set of duplicated orthogroups present in any of the comparison genomes
        all_duped_orthos_in_outgroup = self.get_all_dups(comparision_accessions)
    
        #get set of orthogroups that are duplicated in all input genomes
        input_common_dups = self.dups_common_to_group(input_accessions)
    
        return input_common_dups.difference(all_duped_orthos_in_outgroup) 
    
    
    
    
 
    def get_HOG_annot_for_genome(self, HOG, assemb_accession):
        'return df with annotations for each protein in the given HOG for input genome'       
        annot_df = pd.read_csv(f'./data/genome_annotations/{assemb_accession}_annotations.csv').set_index('HOGs')
        if HOG in annot_df.index:
            annots = annot_df.loc[HOG, :]
            
        else:
            data = np.full(annot_df.shape[1], 'NA')
            annots = pd.Series(data = data, index= annot_df.columns)
            
            print(f'the HOG {HOG} does not exist in {self.accession_to_name[assemb_accession]}')
        return annots   
    
    

        
        


import json
import os
import pandas as pd
from os.path import dirname
import numpy as np
from jw_utils import parse_gff as pgf

class Orthogroup_Parser():
    
    '''***note: all the folders associated with orthofinder must be in a parent folder 
    called "data", and the current working dir must be the parent of this data folder'''
    
    def __init__(self, path_to_results_folder):
        
        
        self.path_to_results_folder = path_to_results_folder
        self.path_to_main_dir = self.get_path_to_dir()
          
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
        
        
    def get_gff_path(self, accession):
        return os.path.join(self.path_to_main_dir, f'./ncbi_dataset/data/{accession}/genomic.gff')
        
    def make_prot2gene_orthologues(self, accession1, accession2):
        '''
        return a df with the protein names in the orthologues df changed to gene names
        for doing comparisons of two genomesd
        '''
        path = os.path.join(
                self.path_to_results_folder,
                f'Orthologues/Orthologues_{accession1}/{accession1}__v__{accession2}.tsv'
                )
        prot2gene_acc1 = pgf.make_prot2gene_dict(self.get_gff_path(accession1))
        prot2gene_acc2 = pgf.make_prot2gene_dict(self.get_gff_path(accession2))
        df = pd.read_csv(path, sep = '\t')
        def get_gene_col(df, accession, prot2gene_dict):
            org = []
            for ind in df.index:
                prots = df.loc[ind,accession].split(', ')
                t = []
                for prot in prots:
                    prot = 'cds-' + prot.strip()
                    t.append(prot2gene_dict[prot])
                org.append(', '.join(t))
            return org
        
        org1 = get_gene_col(df, accession1, prot2gene_acc1)
        org2 = get_gene_col(df, accession2, prot2gene_acc2)
        df2 = pd.DataFrame()
        df2['Orthogroup'] = df.loc[:,'Orthogroup']
        df2[accession1] = org1
        df2[accession2] = org2
        return df2
        
        
        
        
        
        
    def get_protseqs_from_genetree(self, orthogroup, out_filepath):
        '''
        get all protein sequences present in orhtofinder gene tree
        
        parameters
            orthogroup (str):
            out_filepath (str): 
        
            output: fasta file with protein sequence on one line
        '''
        path_to_gene_tree = os.path.join(self.path_to_results_folder, 
                                            f'Resolved_Gene_Trees/{orthogroup}_tree.txt')
        t = Tree(path_to_gene_tree, format=1)
        proteome_lst = []
        for name in t.get_leaf_names():
            acc = name[:15][::-1].replace('_','.',1)[::-1]
            prot = name[16:]
            #get proteome file
            path_to_proteome = os.path.join(self.path_to_main_dir, f'Proteomes/{acc}.faa')
            seq_dict = self.get_fasta_seq_dict(path_to_proteome)
            prot_seq = seq_dict[prot]
            #proteome_lst.append(f'{name} {prot_seq}\n')
            proteome_lst.append(f'>{name}\n {prot_seq}\n')
        
        with open(out_filepath, 'w') as f:
            for line in proteome_lst:
                f.write(line) 
            
            

    def get_fasta_seq_dict(self, path_to_proteome):
        '''
        return a dict with seqID as key and sequence as value
    
        parameters:
            path_to_proteome (str): path to a fasta proteome
        '''
        seq = ''
        seq_dict = {}
        with open(path_to_proteome, 'r') as f:
            for line in f:
                if line[0] == '>':
                    seq = ''
                    prot_id = line.split(' ')[0][1:].strip()
                    seq_dict[prot_id] = seq
                else:
                    seq = seq + line.strip()
                    seq_dict[prot_id] = seq
        return seq_dict


    
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
            df = pd.read(path_to_counts, sep = '\t')
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
                       
        #path_to_json = './Orthofinder_data/summary_data/AssemblyAccession_to_SpeciesName.json'
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
   
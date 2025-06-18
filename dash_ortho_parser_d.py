#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 22 09:44:51 2022

@author: jonwinkelman
"""
import json
import os
import pandas as pd
from os.path import dirname
import numpy as np
from jw_utils import parse_gff as pgf
from jw_utils import parse_fasta as pfa
from collections import Counter

class DashOrthoParser():
    """
    a class to process orthofinder results for a dash application
    
    object attributes:
    path_to_data_folder (str): path to folder holding all data for app
        

    """
    def __init__(self, path_to_data_folder, tax_level='N0', *args, **kwargs):

        
        #self.path_to_data_folder = path_to_data_folder  # 
        self.path_to_data = path_to_data_folder
        self.path_to_orthogroups = os.path.join(self.path_to_data, 'Orthogroups/Orthogroups.tsv')                          
        self.path_to_gene_counts = os.path.join(self.path_to_data,'Orthogroups/Orthogroups.GeneCount.tsv')        
        self.accession_to_name = self.accession_to_name()
        self.name_to_accession = self.name_to_accession()
        self.accessions = list(self.accession_to_name.keys())
        self.HOG_node = tax_level
        self.path_to_HOG_counts = os.path.join(self.path_to_data,
                                               f'{self.HOG_node}_HOG_counts.tsv')
        self.N_HOG_path = os.path.join(self.path_to_data, f'Phylogenetic_Hierarchical_Orthogroups/{self.HOG_node}.tsv') 
        self.HOG_dir_path = os.path.join(self.path_to_data,'Phylogenetic_Hierarchical_Orthogroups')     
        self.HOGs = list(pd.read_csv(self.N_HOG_path, sep='\t', usecols=['HOG'])['HOG'])

           
    def all_prots_in_HOG(self, HOG):
        'return a list of all proteins in a given HOG e.g. GCF_000332095_2_WP_008307711.1'
        HOG_df = pd.read_csv(self.N_HOG_path, sep = '\t', low_memory=False).set_index('HOG')
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
        df = pd.read_csv(self.N_HOG_path, sep = '\t', usecols=['HOG', 'OG']).set_index('HOG')
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
      
    
    def prot_counts_in_all_HOGs(self):
        "Return a df containing, for each accession, the number of proteins in each HOG"
        
        tax_level = self.HOG_node
        HOG_proteins = self.N_HOGs_proteins(tax_level)
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

    
    def check_full_id_format(self, s):
        "checks that s is in the format GC#_########_#_<prot_id>"
        
        return s.startswith("GC") and len(s) > 14 and s[13] == "_"
    
    def prot_counts_in_HOG(self, HOG):
        """Return protein counts within a given HOG accross proteomes, fast but assumes NCBI accessions!!
        
        ***This assumes that protein name has NCBI assembly accession in beginning of ID
        """ 
        all_prots_in_HOG = self.all_prots_in_HOG(HOG)
        #check that beginning of prot ID is a NCBI assembly accession
        if not self.check_full_id_format(all_prots_in_HOG[0]):
            raise ValueError(f'{all_prots_in_HOG[0]} does not appear to be formatted properly')
                
        accessions_present = [full_id[:15][::-1].replace('_', '.', 1)[::-1] for full_id in all_prots_in_HOG]
        counts = Counter(accessions_present)
    
        for acc in self.accessions:
            if acc not in accessions_present:
                counts[acc] = 0
        return counts
 

    def get_gene_counts(self):
        return pd.read_csv(self.path_to_gene_counts, sep='\t').set_index('Orthogroup')

    
    
    def accession_to_name(self):
        'return a dictionary with accessions as keys and names as values'
        path_to_json = os.path.join(self.path_to_data, 'summary_data/AssemblyAccession_to_SpeciesName.json')             
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
    

    
    def HOG_protnames_in_genome(self, HOG, species_accession):
  
        ser=pd.read_csv(self.N_HOG_path, sep='\t', usecols=['HOG', species_accession], low_memory=False).set_index('HOG').loc[HOG,:]
 
        if ser.shape[0] > 0: 
            if type(ser.values[0]) == str:
                return ser.values[0].split(', ')
            else:
                return []
                
        else:
            return []

    def HOG_genenames_in_genome(self, HOG, species_accession):
        prot_list = self.HOG_protnames_in_genome(HOG, species_accession)
        genes = []
        if prot_list:
            path_to_gff3 = os.path.join(self.path_to_data, f'ncbi_dataset/data/{species_accession}/genomic.gff')
            prot2gene_dict = pgf.make_prot2gene_dict(path_to_gff3)
            for prot in prot_list:
                prot = 'cds-'+prot
                genes.append(prot2gene_dict[prot].replace('gene-', ''))
        return genes
    
    
    
    def HOG_proteins_in_genome(self, HOG, species_accessions):
        """Return series with  numer of proteins in a given hog in each genome
        
        HOG (str): name of a single HOG
        species_accessions (str or list):
        """
        #
        
        if type(species_accessions) == list:
            col_list = species_accessions
        if type(species_accessions) ==str:
            col_list = [species_accessions]
        else:
            col_list = list(species_accessions)
        col_list.append('HOG')
        df=pd.read_csv(self.N_HOG_path, sep='\t', usecols=col_list, low_memory=False)
        series = df.set_index('HOG').loc[HOG, :]#.fillna('')
        copies = {}
        for acc in series.index:
            if type(series[acc]) == str:
                num_copies = len(series[acc].split(','))
                copies[acc] = num_copies
            else:
                copies[acc] = 0
        
        return pd.Series(copies) 
    
 

 

    
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

    def protein_to_HOG_df(self, accession):
        "return the HOG that a protein belongs to within a genome"
    
        df = pd.read_csv(self.N_HOG_path, sep ='\t', usecols=['HOG',accession],low_memory=False).dropna()
        df["protein_id"] = df[accession].str.split(", ")
        df = df.explode("protein_id", ignore_index=True).drop(accession, axis=1).set_index('protein_id')
        return df

    
  
    def N_HOGs_proteins(self, tax_level=None, HOG_list=None, accession_list=None):
        """return nested dict {accession:{HOG:proteins}}

        Parameters:
        HOG_list (list): List of HOGs to include, if Null, all HOGs are included
        accession_list (list): List of accessions to include, if Null, all accessions are included
        """
        if not tax_level:
            tax_level = self.HOG_node
        if not tax_level:
            tax_level = 'N0'
        path = os.path.join(self.HOG_dir_path, tax_level + '.tsv')
        HOGs_df = pd.read_csv(path, sep = '\t', low_memory=False)
        HOGs_df = HOGs_df.set_index('HOG')
        if HOG_list:
            HOGs_df = HOGs_df.loc[HOG_list,:]
        if accession_list:
            HOGs_df = HOGs_df.loc[:,accession_list]
        else: accession_list=self.accessions
        HOG_proteins = {acc:{} for acc in accession_list}

        for assembly_identifier in accession_list: #for each accession column
            null_filt = HOGs_df[assembly_identifier].notnull()
            populated_OGs = HOGs_df.loc[null_filt, assembly_identifier]
            populated_OGs = populated_OGs.to_dict()
            for HOG in populated_OGs:
                HOG_proteins[assembly_identifier][HOG] = populated_OGs[HOG].split(',')
        return HOG_proteins
    
       
    def get_unique_orthogroup_dups(self, input_accessions, comparision_accessions):
        #get total set of duplicated orthogroups present in any of the comparison genomes
        all_duped_orthos_in_outgroup = self.get_all_dups(comparision_accessions)
    
        #get set of orthogroups that are duplicated in all input genomes
        input_common_dups = self.dups_common_to_group(input_accessions)
    
        return input_common_dups.difference(all_duped_orthos_in_outgroup) 
    
    
    
    def get_HOG_annot_for_genome(self, HOG, assemb_accession):
        'return df with annotations for each protein in the given HOG for input genome'
        path_to_data = os.path.join(self.path_to_data,f'genome_annotations/{assemb_accession}_annotations.csv' )
        annot_df = pd.read_csv(path_to_data).set_index('HOGs')
        if HOG in annot_df.index:
            annots = annot_df.loc[HOG, :]
            
        else:
            data = np.full(annot_df.shape[1], 'NA')
            annots = pd.Series(data = data, index= annot_df.columns)
            
            print(f'the HOG {HOG} does not exist in {self.accession_to_name[assemb_accession]}')
        return annots

    def _get_common_gene_names(self, prot_dict, gene_dict):
        Common_names = []
        for prot in prot_dict:
            gene_obj = gene_dict.get(prot_dict[prot].Parent)
            if gene_obj:
                Common_names.append(gene_obj.Name)
            else:
                Common_names.append(None)
        return Common_names

    def make_genome_annotation_df(self, assembly_accession, get_common_names=False):
        path_to_gff3 = os.path.join(self.path_to_data, f'ncbi_dataset/data/{assembly_accession}/genomic.gff')
        gene_dict = pgf.make_seq_object_dict(path_to_gff3, feature_type='gene')
        seq_dict = pgf.make_seq_object_dict(path_to_gff3, feature_type='CDS')
        df_dict = {
                'IDs':[seq_dict[key].ID for key in seq_dict.keys()],
                'starts':[int(seq_dict[key].start) for key in seq_dict.keys()],
                'ends':[int(seq_dict[key].end) for key in seq_dict.keys()],
                'strand': [seq_dict[key].strand for key in seq_dict.keys()],
                'products':[seq_dict[key].product for key in seq_dict.keys()],
                'Parents':[seq_dict[key].Parent for key in seq_dict.keys()],
                }
        if get_common_names:
            df_dict['Common_names'] = self._get_common_gene_names(seq_dict, gene_dict)
        df = pd.DataFrame(df_dict)
        df = df.set_index('IDs')
        #add orthogroup and HOG column to df
        gene_to_orthogroup_dict = self.gene_to_orthogroup_dict(assembly_accession)
        gene_to_HOG_dict = self.gene_HOG_dict(assembly_accession)
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
        return df   
    

    def HOG_gene_dict(self, species_accession): 
       #make df with orhtogroup and all genes in a given organism
        df = pd.read_csv(self.N_HOG_path, sep = '\t', usecols = ["HOG", species_accession], low_memory=False)
        df = df.astype(str)
        df = df.set_index('HOG')
        HOG_gene_dict = {}
        for HOG in df.index:
            HOG_gene_dict[HOG] = df.loc[HOG,species_accession].split(', ')
        return HOG_gene_dict
 


    def gene_HOG_dict(self, species_accession, feature_type='protein'):
        """
        Right now this *stupidly* returns protein to HOG by default. 
        
        feature_type (str): can be 'protein' or 'gene'
        """
        gene_HOG_dict = {}
        HOG_gene_dict = self.HOG_gene_dict(species_accession)
        for HOG, proteins in HOG_gene_dict.items():
            proteins =  HOG_gene_dict[HOG]
            if proteins[0] != 'nan':
                for protein in proteins:
                    gene_HOG_dict[protein.strip()] = HOG
        if feature_type=='gene':
            path_to_gff3 = os.path.join(self.path_to_data, f'ncbi_dataset/data/{species_accession}/genomic.gff')
            prot2gene_dict = pgf.make_prot2gene_dict(path_to_gff3)
            prot2gene_dict = {p.replace('cds-',''):g.replace('gene-','') for p,g in prot2gene_dict.items()}
            gene2HOG_dict = {prot2gene_dict[prot]:HOG for prot, HOG in gene_HOG_dict.items()}
            return gene2HOG_dict
            
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
    
 

    def get_HOG_protein_seqs(self, HOG, accession_length=15):
        """Return df with proteinID as index, accession, HOG and protein sequence """
        prot_names = self.all_prots_in_HOG(HOG)
        d = {}
        for prot_name in prot_names:
            acc = prot_name[:accession_length][::-1].replace('_', '.',1)[::-1]
            prot_id = prot_name[(accession_length+1):]
            path_to_fasta = os.path.join(self.path_to_data, f'Proteomes/{acc}.faa')
            seq = pfa.get_seq_dict(path_to_fasta).get(prot_id)
            d[prot_name] = [acc, prot_id, HOG, seq]
        df = pd.DataFrame.from_dict(d, orient='index')
        df.columns = ['Accession', 'Prot_id', 'HOG', 'Protein_seq']
        return df
    
    
    

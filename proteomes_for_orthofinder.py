#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan  1 14:36:07 2022

@author: jonwinkelman
workflow:
    1. get assembly accessions for ncbi genomes and save them to the directory you are working in 
    2. as a default, the directory structure will be created within the directory you are currently in
    3. things are easiest if all genomes are in GCF_#########.version format. Otherwise several things have 
        to be done manually
        a. a proteome file (filename.faa) has to be manually added to the ./Proteomes directory
        b. a directory containing the file: ./ncbi_datasets/data/filename.version/genomic.gff
        c. the species name and accession have to be manually added to ./summary_data/AssemblyAccession_to_SpeciesName.json
            i. to manually add: run: get_genome_summmary(accession_list, manual_proteomes=True) and follow prompts
    running Orthofinder:
        I am currently running orthofinder from a docker container due to problems in the conda environment.
        running the newist version of orthofinder is supposed to fix this
        running docker container:
        1. there is a google meet recording with Bradford going through the process 2022-02-07 08:44 GMT-8 in email and google drive:
            https://drive.google.com/file/d/1Id38aTNVMva9okclRxjoEavfSuaxarqY/view
        2. in breif: run the below command in terminal
            docker run --ulimit nofile=1000000:1000000 -it --rm -v "/path/to/data/directory":/data davidemms/orthofinder:2.5.4 bash


for downloading large datasets from NCBI, download dehydrated data packages using the cli as such:
    eg in the termal:
    datasets download genome accession --inputfile accessions.txt --dehydrated --filename my_dehydrated_dateset.zip --exclude-gff3  --exclude-rna --exclude-protein
    unzip my_dehydrated_dateset.zip -d my_rehydrated_dateset
    datasets rehydrate --directory my_rehydrated_dateset/

"""
import os
from zipfile import ZipFile
import json
try:
    import ncbi.datasets
except ImportError:
    print('ncbi.datasets module not found. To install, run `pip install ncbi-datasets-pylib`.')

def check_for_duplicates(acc_list):
    if len(acc_list) != len(set(acc_list)):
        raise Exception('There are duplicate values in the input list, exiting...') 
    else:
        pass

def remove_duplicates(acc_list):
    old_len = len(acc_list)
    acc_list_deduped = list(set(acc_list))
    new_len = len(acc_list_deduped)
    print(f'{old_len-new_len} accessions were removed because they were duplicates')
    return acc_list_deduped
        
    
def get_accessions(path_to_accessions = None):
    if path_to_accessions:
        path_to_accessions = path_to_accessions
    else:
        path_to_accessions = 'accessions.txt'
    with open(path_to_accessions, 'r') as f:
        accessions = [line.strip() for line in f]
        return list(set(accessions))


def get_genome_data(acc_list = None,  include_annotation_type = ['PROT_FASTA', 'GENOME_GFF' ]): 
    if acc_list:
        acc_list = acc_list
    else:
        acc_list = get_accessions()
    check_for_duplicates(acc_list)
    #['DEFAULT', 'GENOME_GFF', 'GENOME_GBFF', 'GENOME_GB', 'RNA_FASTA', 'PROT_FASTA', 'GENOME_GTF', 'CDS_FASTA']
    #include_annotation_type = ['PROT_FASTA', 'GENOME_GFF' ]
    api_instance = ncbi.datasets.GenomeApi(ncbi.datasets.ApiClient())
    api_response = api_instance.download_assembly_package(
            acc_list,
            include_annotation_type=include_annotation_type,
            # Because we are streaming back the results to disk, 
            # we should defer reading/decoding the response
            _preload_content=False
        )
    with open('assemblies.zip', 'wb') as f:
        f.write(api_response.data)
        
    

def unzip_ncbi_files(path_to_file=None):
    if path_to_file:
        pass
    else:
        path_to_file = 'assemblies.zip'
    with ZipFile(path_to_file, 'r') as zipObj:
        listOfFileNames = zipObj.namelist()
        for file in listOfFileNames:
            if file.endswith('genomic.gff'):
                zipObj.extract(file)

            elif file.endswith('protein.faa'):
                zipObj.extract(file) 
                
            elif file.endswith('genomic.fna'):
                zipObj.extract(file)
                
def grab_proteomes(path_to_ncbi_data=None):
    import os
    import shutil
    if not path_to_ncbi_data:
        path = './ncbi_dataset/data'
    else:
        path=path_to_ncbi_data
    os.makedirs('./Proteomes')
    folders =  os.listdir(path)
    folders = [fold for fold in folders if fold.startswith('GC')]
    print(f'{len(folders)} assemblies found in {path}')
    failed = []
    for folder in folders:
        proteome_file = [file for file in os.listdir(os.path.join(path,folder)) if file.endswith('.faa')]
        if len(proteome_file) == 1:
            path_to_proteome = os.path.join(path,folder,proteome_file[0])
            shutil.copy(path_to_proteome, './Proteomes')
            os.rename('./Proteomes/protein.faa', './Proteomes/' + folder + '.faa')
        else:
            print(f'{folder} did not have 1 and only 1 .faa file')
            failed.append(folder)
    proteomes = os.listdir('./Proteomes')
    accessions = [acc.replace('.faa', '') for acc in proteomes]
    print(f'{len(folders)} folders were processed and there are {len(proteomes)} proteomes  in "./Proteomes" ')
    return accessions
       
            
def write_AssemblyAccession_to_SpeciesName(acc_2_name_dict, output_file_path=None):
    'writes json file with {assemb_accession:name}'
    if output_file_path:
        output_file_path = output_file_path
    else:
        output_file_path = './summary_data/AssemblyAccession_to_SpeciesName.json'
        if not os.path.isdir('./summary_data'):
            os.makedirs('./summary_data')
        if os.path.isfile(output_file_path):
            raise Exception('File already exists')
    with open(output_file_path, 'w') as f:
        json.dump(acc_2_name_dict, f)
    
           


def get_summary_loop(accession_list, write_file = True, output_file_path=None):
    '''
    pulls summaries from NCBI datasets and builds python dictionary

    returns (tuple): assembly_to_name dict and an much more complex assembly dict
    
    '''
    acc_2_name = {}
    assembly_dict = {}
    fail_list = []
    for acc in accession_list:
        acc_lst = [acc]
        try:
            temp_ass_2_name, temp_assembly_dict = build_summary_data(acc_lst)
            acc_2_name.update(temp_ass_2_name)
            assembly_dict.update(temp_assembly_dict)
        except:
            fail_list.append(acc)
            print(f'the accession {acc} failed to download from ncbi')
            
        
    #checks
    success = len(acc_2_name.keys())
    attempted = len(accession_list)
    if success != attempted: 
        print(f'{success} assembly summaries were successfully downloaded, {attempted - success} failed')
        
   
    #get other file path if entered
    if write_file:
        write_AssemblyAccession_to_SpeciesName(acc_2_name, output_file_path=output_file_path)
    
    return acc_2_name, assembly_dict
        
       

def build_summary_data(acc_list, manual_proteomes = False):
    
    """only run this with 'get_summary_loop', Do not call directly or it will only return partial results
    acc_list (lst): list of ncbi assembly accessions. 
    return (tuple):
        tuple[0] (dict): assembly accession:organism name
        tuple[1] (dict): assembly_dict[acc_accession]['org']['strain']
    """
    
    api_instance = ncbi.datasets.GenomeApi(ncbi.datasets.ApiClient()) 
    assembly_dict={}  
    acc_2_name = {}
    genome_summary = api_instance.assembly_descriptors_by_accessions(acc_list)
    for i, assembly in enumerate(genome_summary['assemblies']):
        assembly_dict[acc_list[i]] = assembly['assembly']
        org_name = assembly_dict[acc_list[i]]['org']['sci_name'].split(' ')[:2] #keep only genus and species
        strain = assembly_dict[acc_list[i]]['org'].get('strain') #get strain info if available
        name = f'{org_name[0]} {org_name[1]} {strain}'
        acc_2_name[acc_list[i]] = name
    return acc_2_name, assembly_dict





def run_all(acc_list, manual_proteomes=False):
    get_genome_data(acc_list)
    unzip_ncbi_files() 
    grab_proteomes() 
    return get_summary_loop(acc_list)          
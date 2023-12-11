def rgb_to_hex():
    return {'rgba(132,137,145,1)':'#848991',
                'rgb(0,208,132)':'#00d084',
                'rgb(0,122,255)':'#007aff',
                'rgb(171,184,195)':'#abb8c3',
                'rgb(255,105,0)':'#ff6900',
                'rgb(252,185,0)':'#fcb900',
                'rgb(123,220,181)':'#7bdcb5',
                'rgb(142,209,252)':'#8ed1fc',
                'rgb(6,147,227)':'#0693e3',
                'rgb(155,81,224)':'#9b51e0',
                'rgb(207,46,46)':'#cf2e2e',
                'rgb(247,141,167)':'#f78da7'}




def _copynumber_bargraph_data_dict(HOGs, dop_obj, id_type='accession'):
    """Return a nested dict 'data_dict' {HOG{leaf_id:counts}} for input into create_tree_w_bargraphs function
    
    Parameters:
    HOGs (list): List of one or more HOGs
    dop_obj (dash_orthoparser oject): must be generated atthe same Hierarchical level as the input HOG
    id_type (str): Can be 'accession' or 'name'. Determines whether the dict keys will be accessions or sci. names.
    """
    
    data_dict =  {}
    data_dict_names = {}
    for HOG in HOGs:
        data_dict[HOG] = dop_obj.HOG_proteins_in_genome(HOG, dop_obj.accessions).to_dict()
    if id_type =='accession':
        return data_dict

    if id_type == 'name':
        for HOG in HOGs:
            data_dict_names[HOG] = {}
            for  acc, count in data_dict[HOG].items():
                name = dop_obj.accession_to_name[acc]
                data_dict_names[HOG][name] = data_dict[HOG][acc]
        return data_dict_names

def make_count_dict(tree, HOGs, dop_obj):
    """Return tuple with [HOG1, HOG2, ...] and a dict {genome_accession:[counts, for, HOG1, HOG2,...]}

    parameters:
    tree (Bio.Phylo.Newick.Tree): Species tree with leaf names that match keys in the second level of the nested count dict
    HOGs (dict): HOG names you wish to 
    dop_obj (orthofinder_utils.dash_ortho_parser.DashOrthoParser): DashOrthoParser object for calculating counts from orthofinder results
    
    """
    HOG_dict = _copynumber_bargraph_data_dict(HOGs, dop_obj, id_type='accession')
    count_d = {leaf_name:[] for leaf_name in [cl.name for cl in tree.get_terminals()]}
    HOG_list = []
    for HOG, d in HOG_dict.items():
        HOG_list.append(HOG)
        for acc in d:
            count_d[acc].append(d[acc])
    return HOG_list, count_d


def make_binary_count_dict(tree, HOGs, dop_obj):
    """Return tuple with [HOG1, HOG2, ...] and a dict {genome_accession:[binary counts, for, HOG1, HOG2,...]}
    *binary counts is just 0 of 1 for presence or absence of ortholog(s) in genome
    
    parameters:
    tree (Bio.Phylo.Newick.Tree): Species tree with leaf names that match keys in the second level of the nested count dict
    HOGs (dict): HOG names you wish to 
    dop_obj (orthofinder_utils.dash_ortho_parser.DashOrthoParser): DashOrthoParser object for calculating counts from orthofinder results
    
    """

    HOG_list, count_dict = make_count_dict(tree, HOGs, dop_obj)
    binary_count_dict = {}  
    for HOG, count_list in count_dict.items():
        binary_count_dict[HOG] = [0 if count == 0 else 1 for count in count_list]      
    return HOG_list, binary_count_dict



def _make_field_labels(name_list):
    """return itol line with field labels pulled from from name list"""
    
    field_labels='FIELD_LABELS,'
    for name in name_list:
        field_labels =f'{field_labels}{name},'
    field_labels=field_labels.strip(',')
    return field_labels


def _make_field_colors(name_list,hexcolors=None):
    """Return itol line with hexcolors for each field"""
    if not hexcolors:
        hexcolors = ['#58D68D','#F4D03F','#F5B041','#AAB7B8','#566573','#A93226','#EC7063', '#A569BD', '#5DADE2','#48C9B0']
    field_colors = 'FIELD_COLORS'
    for i, _ in enumerate(name_list):
        field_colors = f'{field_colors},{hexcolors[i]}'
    field_colors=field_colors.strip(',')
    return field_colors


def _make_legend_colors(name_list,hexcolors=None):
    """Return itol line with hexcolors for each field"""
    if not hexcolors:
        hexcolors = ['#58D68D','#F4D03F','#F5B041','#AAB7B8','#566573','#A93226','#EC7063', '#A569BD', '#5DADE2','#48C9B0']
    legend_colors = 'LEGEND_COLORS'
    for i, _ in enumerate(name_list):
        legend_colors = f'{legend_colors},{hexcolors[i]}'
    legend_colors=legend_colors.strip(',')
    return legend_colors

def _make_legend_labels(name_list):
    """return itol line with legend labels pulled from name list"""
    
    legend_labels='LEGEND_LABELS,'
    for name in name_list:
        legend_labels =f'{legend_labels}{name},'
    legend_labels=legend_labels.strip(',')
    return legend_labels


def make_multi_itol_bargraph_dataset(outfile_path, count_dict, name_list, dataset_label='dataset_label',
                                     color='#848991', legend_title='Dataset legend', hexcolors=None,):
    """
    Add data to simple bargraph template for itol dataset adn write to file.

    Parameters
    ----------
    template_path : str
        path to template file from itol.
    file_path : str
        path for new file that is to be created
    data_dict :dict
        {name:{terminal node name: count}}

    """
    if not hexcolors:
        hexcolors = ['#58D68D','#F4D03F','#F5B041','#AAB7B8','#566573','#A93226','#EC7063', '#A569BD', '#5DADE2','#48C9B0']

    with open(outfile_path, 'w') as f:
        f.write('DATASET_MULTIBAR\n')
        f.write('SEPARATOR COMMA\n')
        f.write(f'DATASET_LABEL,{dataset_label}\n')
        f.write(f'COLOR,{color}\n')
        f.write(f'{_make_field_colors(name_list, hexcolors=hexcolors)}\n')
        f.write(f'{_make_field_labels(name_list)}\n')
        f.write(f'LEGEND_TITLE,{legend_title}\n')
        f.write(f"LEGEND_SHAPES,{','.join(['1' for _ in name_list])}\n")
        f.write(f'{_make_legend_colors(name_list, hexcolors=hexcolors)}\n')
        f.write(f'{_make_legend_labels(name_list)}\n')
        f.write('ALIGN_FIELDS,1\n')
        f.write('DATA\n')

        for name, counts in count_dict.items():
            line = f'{name},'
            for count in counts:
                line = line + str(count)+','
            line=line.strip(',') 
            f.write(line+'\n')

def relabel_itol_treeleafs(tree, relabel_dict, outfile_path):
    """Write itol annotation file to relabel terminal leaves in tree
    
    parameters:
    tree (Bio.Phylo.Newick.Tree): tree that is to be relabeled
    relabel_dict (dict): dict {old_leaf_name:new_leaf_name}
    outfile_path (str): path for new itol annotation file
    """
    import warnings
    tree_leafnames = [cl.name for cl in tree.get_terminals()]
    if len(tree_leafnames) != len(relabel_dict.keys()):
        warnings.warn('The number of tree leafs and the number of dict key names to be replaced are not equal')
    with open(outfile_path, 'w') as f:
        f.write('LABELS\n')
        f.write('SEPARATOR COMMA\n')
        f.write('DATA\n')
        for old_name, new_name in relabel_dict.items():
            f.write(f'{old_name},{new_name}\n')
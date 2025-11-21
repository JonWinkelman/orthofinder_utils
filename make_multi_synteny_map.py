from pathlib import Path
import pandas as pd 
from jw_utils import parse_gff as pgf
from jw_utils import plotly_preferences as pprefs
import plotly.graph_objs as go
from plotly.subplots import make_subplots
from jw_utils import jw_draw_tree as jtw
from tqdm import tqdm
from jw_utils import jw_draw_tree as jdt
from Bio import Phylo

def get_acc_feature_list(tree): 
    ""
    acc_feature_lst = []
    leaf_names = [cl.name for cl in tree.get_terminals()]
    for i, leaf in  enumerate(leaf_names, start=1):
        feature = leaf[16:]
        acc = leaf[:15][::-1].replace('_', '.', 1)[::-1]
        acc_feature_lst.append([acc, feature])
    return acc_feature_lst


def get_tree_from_HOG(dop_obj, HOG):
    """Returns a resolved orthofinder gene tree pruned to the HOG """
    OG = dop_obj.HOG_OG_dict()[HOG]
    OG_tree_fp = Path(dop_obj.path_to_data) / f'Resolved_Gene_Trees/{OG}_tree.txt'
    tree = Phylo.read(str(OG_tree_fp), format='newick')
    prune_OG_tree_to_HOGtree(tree, dop_obj, HOG)
    return tree



def prune_OG_tree_to_HOGtree(OG_tree, dop_obj, HOG):

    """Prune away leaves from tree (in place) that are in OG but not in corresponding HOG"""
    OG = dop_obj.HOG_OG_dict()[HOG]
    HOG_prots = dop_obj.all_prots_in_HOG(HOG)
    OG_prots = dop_obj.all_prots_in_orthogroup(OG)
    to_prune=set(OG_prots).difference(HOG_prots)
    for leaf in to_prune: 
        OG_tree.prune(leaf)


def genetreeID_2_orgname(tree, dop_obj):
    """Return dict with organism name
    
    *assumes standard OrthoFInder leaf name format: 'GCF_964211495_1_WP_368892597.1'
    """
    tree_leaves = [cl.name for cl in tree.get_terminals()]
    relable_d = {}
    for leafname in tree_leaves: 
        raw_acc = leafname[:15]
        acc     = raw_acc[::-1].replace('_', '.', 1)[::-1]
        orgname = dop_obj.accession_to_name[acc].split(' ')
        orgname = '_'.join([f for f in orgname])
        newname = f"{leafname[16:]}_{orgname}"
        relable_d[leafname] = newname
    return relable_d


def add_annotations(fig, xcoord, ycoord, text, fontsize=3):
    fig.add_annotation(
        x=xcoord,
        y=ycoord,
        text=text,
        showarrow=False,
        xanchor="left",    # anchor text on the left
        yanchor="middle",  # vertically centered
        xref="x",
        yref="y", 
        font=dict(size=fontsize) 
    )
    return fig



def get_connector_traces(tree): 
    """Returns list of go.Scatter traces from each node to  edge of graph (length of longest leaf)"""

    traces = []
    term_ycoords = pd.Series({cl.name:v for cl,v in jdt.get_y_coordinates(tree).items() if cl.is_terminal()})
    term_xcoords = pd.Series({cl.name:v for cl,v in jdt.get_x_coordinates(tree).items() if cl.is_terminal()})
    xmax = term_xcoords.max()
    for name, leaf_ycoord in term_ycoords.items(): 
        connector_trace = go.Scatter(y=[leaf_ycoord, leaf_ycoord], 
                                     x = [xmax, term_xcoords[name]], 
                                     mode='lines', 
                                     line={'dash':'dash',
                                          'width':0.2, 
                                          'color':'rgb(100,100,100)'}) 
        traces.append(connector_trace)
    return traces







def make_tree_synteny_fig(dop_obj, tree, acc_feature_lst, num_genes = 20, x_spread =8000, height=5000,
                           HOGs_to_highlight=None, annotate_leaves=False, rename=False,
                            feature_line_widths=0.1, **plot_opts):
    """
    dop_obj (DashOrthoParser object):
    acc_feature_lst (list of lists or tuples): internal list contains [[genome accessio, feature ID],[],,,]  feature ID is the protein that the figure will 
    be centered around. Presumably they are all orthologs to each other. e.g. [['GCF_003024515.2', 'WP_107010098.1'],
                     ['GCF_964211495.1', 'WP_368892597.1'],...]
    num_genes (int): number of genes to process on each side of feature ID in which each fig is centered
    x_spread (int): number of nts on each side of feature ID. 
    height (int): Figure height
    HOGs_to_highlight (dict): e.g.{HOG1:'rgb(100,100,100)', HOG2:'rgb(200,100,100)',...}
    annotate_leaves : bool
            default False. If True adds permament leaf name at end of leaf node

    feature_line_widths : float
            width of lines in synteny figure
    **kwargs :
        Additional keyword arguments passed to go.Scatter.
        Supported kwargs:
            intern_node_label : str
                label to appear in hover data. "name", "confidence", ...fields set by newick file and Bio.Phylo
            connector_lines : str
                lines connectin leaf nodes to synteny plot
            in_node_size : float
                Internal node size on pylo tree.
            nt_node_sizeame : float
                Terminal node size on pylo tree.
            opacity : float
                Opacity of the trace.
            hover_text : iterable of strings
                text to appear when hovering over nodes on tree
            cl_to_highlight : clade object or iterable of clade objects,
                should be internal node
            node_color_dict : dict
                for coloring nodes, if node name not in dict will be given
                defaulg color. e.g. {node_name:'rgb(x,x,x)'}

            


    """


    fig = make_subplots(
        rows =1, 
        cols =2, 
        shared_xaxes=False, 
        shared_yaxes=True,
        horizontal_spacing=0
        )
    
 
    tree_trace  = jtw.create_tree(tree, **plot_opts)
    ycoords = {cl.name:v for cl,v in jdt.get_y_coordinates(tree).items()}
    ## update tree x-range to remove whitespace
    all_xcoords = pd.Series({cl.name:v for cl,v in tree.depths().items()})
    tree_xspan=all_xcoords.max()-all_xcoords.min()
    pad = tree_xspan*0.01
    fig.update_xaxes(range=[all_xcoords.min()-pad, all_xcoords.max()+pad], row=1, col=1)  # subplot (1,1)
    # connector traces
    connector_traces = get_connector_traces(tree)
    fig.add_traces(connector_traces, rows=[1 for r in connector_traces], cols=[1 for c in connector_traces])
    # add permanant leaf names
    if annotate_leaves:
        pass
    fig.add_trace(tree_trace['data'][0], row=1, col=1)
    fig.update_layout(tree_trace['layout'])
    for acc_feature in  tqdm(acc_feature_lst):
        
        acc, feature = acc_feature
        fts = int(num_genes/2)
        leafname = f"{acc_feature[0].replace('.', '_')}_{feature}"
        trimmed_df = get_df_for_feature(dop_obj, acc, feature=feature, fts=fts,)
        trim = trimmed_df.loc[feature, 'starts'] - x_spread
        trimmed_df['starts'] = trimmed_df['starts'] - trim
        trimmed_df['ends'] = trimmed_df['ends'] - trim
        tfig, xrange = make_map(trimmed_df, feature_name=feature,  x_spread=x_spread)
        for t in tfig.data:
            if "line" in t:
                t.line.width = feature_line_widths
            if t['text']:
                HOG=t['text'].split('<br>')[-1].split('HOG: ')[-1]
                if HOGs_to_highlight:
                    if HOGs_to_highlight.get(HOG): 
                        t.line.color = HOGs_to_highlight.get(HOG)
            ycoord = ycoords[leafname]
            t['y'] = [y+ycoord for y in t['y']]
            fig.add_trace(t, row=1, col=2)
    fig.update_xaxes(range=xrange, row=1, col=2)
    fig = blank_layout(fig)
    if annotate_leaves:
        relable_d = genetreeID_2_orgname(tree, dop_obj)
        for cl in tree.get_terminals():
            xcoord, ycoord = all_xcoords[cl.name], ycoords[cl.name]
            add_annotations(fig, xcoord, ycoord, relable_d[cl.name] )
    return fig



def blank_layout(fig):
    fig.update_xaxes(
        showgrid=False,    # hide gridlines
        zeroline=False,    # hide x=0 line
        showline=False,    # hide axis line
        showticklabels=False,  # hide tick labels
        ticks=""
    )

    fig.update_yaxes(
        showgrid=False,
        zeroline=False,
        showline=False,
        showticklabels=False,
        ticks=""
    )
    fig.update_layout(
        plot_bgcolor='rgba(0,0,0,0)',
        paper_bgcolor='rgba(0,0,0,0)',
        margin={'t':0, 'b':0, 'l':0, 'r':0},
        showlegend = False,
    )

    # fig.update_xaxes(showticklabels=False)
    # fig.update_yaxes(showticklabels=False)
    return fig

def make_multi_synteny_map(dop_obj, acc_feature_lst, num_genes = 20, x_spread=8000, height=5000, HOGs_to_highlight=None):
    """
    dop_obj (DashOrthoParser object):
    acc_feature_lst (list of lists or tuples): internal list contains [[genome accessio, feature ID],[],,,]  feature ID is the protein that the figure will 
    be centered around. Presumably they are all orthologs to each other. e.g. [['GCF_003024515.2', 'WP_107010098.1'],
                     ['GCF_964211495.1', 'WP_368892597.1'],...]
    num_genes (int): number of genes to process on each side of feature ID in which each fig is centered
    x_spread (int): number of nts on each side of feature ID. 
    height (int): Figure height
    HOGs_to_highlight (dict): e.g.{HOG1:'rgb(100,100,100)', HOG2:'rgb(200,100,100)',...}
    """
    fig = make_subplots(rows = len(acc_feature_lst), cols=1, shared_xaxes=False )
    for i, acc_feature in  tqdm(enumerate(acc_feature_lst, start=1)):
        acc, feature = acc_feature
        fts = int(num_genes/2)
        trimmed_df = get_df_for_feature(dop_obj, acc, feature=feature, fts=fts,)
        tfig, xrange = make_map(trimmed_df, feature_name=feature,  x_spread=x_spread)
        for t in tfig.data:
            if t['text']:
                HOG=t['text'].split('<br>')[-1].split('HOG: ')[-1]
                if HOGs_to_highlight:
                    if HOGs_to_highlight.get(HOG): 
                        t.line.color = HOGs_to_highlight.get(HOG)
            fig.add_trace(t, row=i, col=1)
            fig.update_xaxes(range=xrange, row=i, col=1)

    fig = blank_layout(fig)
    fig.update_layout(height=height)
    return fig




def flip_block(df):
    """
    Flip a genomic block (reverse coordinates and strand)
    for a dataframe containing columns:
    ['ID', 'starts', 'ends', 'strand', 'product', 'gene_ID', 'HOG']
    """
    
    # Determine boundaries of the block
    block_min = df[['starts', 'ends']].min().min()
    block_max = df[['starts', 'ends']].max().max()

    df = df.copy()  # to avoid modifying original dataframe

    # Flip coordinates
    df['new_starts'] = block_max - df['ends']
    df['new_ends']   = block_max - df['starts']

    # Flip strand
    df['new_strand'] = df['strand'].map({'+': '-', '-': '+'}).fillna(df['strand'])

    # Replace old columns with new ones
    df['starts'] = df['new_starts']
    df['ends']   = df['new_ends']
    df['strand'] = df['new_strand']

    # Drop temp columns
    df = df.drop(columns=['new_starts', 'new_ends', 'new_strand'])

    # Sort by the new starts coordinate
    df = df.sort_values('starts')

    return df



def get_df_for_feature(dop_obj, assembly_accession, feature, fts=15):
    """Return a local df containing features near the passed feature
    
    feature (str): gene locus tag (no "gene-" prefix btw...)
    fts (int): num features to show in map on each side of of feature of interest 
    """

    df = dop_obj.make_genome_annotation_df(assembly_accession)
    df.index=df.index.str.replace('cds-', '')
    df['Parents']=df['Parents'].str.replace('gene-', '')
    if feature in df.index:  

        f_index = df.index.get_loc(feature)
        if f_index>=fts and f_index<= (len(df.index)-fts):
            trimmed_df = df.iloc[(f_index-fts):(f_index+fts),:]
        
        elif f_index-fts<0:
            trimmed_df = df.iloc[:(f_index+fts),:]
        elif f_index + fts > df.shape[0]:
            trimmed_df = df.iloc[(f_index-fts):,:]
    else:
        print(f'{feature} not in in the dataframe index, e.g. index = {df.index[0]}')
        trimmed_df = df.iloc[:10,:]
    if trimmed_df.loc[feature, 'strand'] == '-': 
        trimmed_df=flip_block(trimmed_df)
    return trimmed_df


def make_arrow_trace(protein, par, strand, prod, start, end, br, HOG, **kwargs): #color=None,
    if kwargs['color']:
        color = kwargs['color']
    else:
        color = 'rgb(100,100,100)'

    'return a trace with an arrow drawn for the input feature'
    
    trace = go.Scatter(
            x = [start, start, br, br, end,br,br,start,start],
            y = [0, 0.25,0.25,0.5,0,-0.5,-0.25,-0.25,0],
            mode = 'lines',
            fill='toself',
            line = {'width': 1.5,
                    'color':color},
            name = protein,
            text = f'{protein}<br><br>{par}<br>strand: {strand}<br>{prod}<br>start: {start}<br>end: {end}<br>HOG: {HOG}',
            hoverinfo = 'text')
    return trace


def make_map(trimmed_df, feature_name = None, height = 150, yrange = [-1,1], x_spread = None, bgcolor='rgb(250,250,250)'):
    traces = []
    xrange = [0,0]
    #make traces for each feature


    for i,protein in enumerate(trimmed_df.index):
        #hoverinfo variables
        par =  trimmed_df.loc[protein,'Parents']
        prod = trimmed_df.loc[protein,'products']
        start = trimmed_df.loc[protein, 'starts']
        end = trimmed_df.loc[protein, 'ends']
        strand = trimmed_df.loc[protein,'strand']
        #orthogroup = trimmed_df.loc[protein,'Orthogroups']
        HOG = trimmed_df.loc[protein,'HOGs']
        #make backbone trace between features
        if i< (len(trimmed_df.index)-1):
            next_orf_start = trimmed_df.iloc[(i+1),0]
            traces.append(go.Scatter(x = [end, next_orf_start],
                          y = [0,0],
                          mode = 'lines',
                          line = {'width': 3,
                                  'color':'black'},
                          showlegend = False,
                          hoverinfo = None))
        if strand == '-':
            start =  end
            end = trimmed_df.loc[protein, 'starts']   
        l = (end-start)     #lenght of the arrow 
        br = start+(0.65*l) #defines where head of arrow is placed
         #make feature traces, highlighting the feature of interest
    
        if protein == feature_name:
            traces.append(make_arrow_trace(protein, par, strand, prod, start, end, br, HOG, color='red'))
            xrange  =[start-x_spread, end+x_spread]
        else:
            traces.append(make_arrow_trace(protein, par, strand, prod, start, end, br, HOG, color='grey'))
        
    dl = {'data':traces,
            'layout':go.Layout(#title = f'local genome map around {feature_name}',
                               paper_bgcolor='rgb(255,255,255)',
                               plot_bgcolor=bgcolor,
                               width = 700,
                               height = height,
                               margin={'t':0, 'b':0, 'l':0, 'r':0},
                               showlegend = False,
            )} 
 
    fig = go.Figure(dl)  
    fig.update_yaxes(showgrid=False, showticklabels=False)
    fig.update_yaxes(range=yrange)
    fig.update_xaxes(range=xrange)
    return fig, xrange
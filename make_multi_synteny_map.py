from pathlib import Path
import pandas as pd 
from jw_utils import parse_gff as pgf
from jw_utils import plotly_preferences as pprefs
import plotly.graph_objs as go
from plotly.subplots import make_subplots
from tqdm import tqdm



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
        tfig, x_spread = make_map(trimmed_df, feature_name=feature,  x_spread=x_spread)
        for t in tfig.data:
            if t['text']:
                HOG=t['text'].split('<br>')[-1].split('HOG: ')[-1]
                if HOGs_to_highlight:
                    if HOGs_to_highlight.get(HOG): 
                        t.line.color = HOGs_to_highlight.get(HOG)
            fig.add_trace(t, row=i, col=1)
            fig.update_xaxes(range=x_spread, row=i, col=1)

    fig.update_layout(height=height, 
                     paper_bgcolor='rgb(255,255,255)',
                     plot_bgcolor='rgb(250,250,250)',
                     margin={'t':0, 'b':0, 'l':0, 'r':0}, 
                     showlegend = False,
                     font={'size':5}
                      
    )
    fig.update_xaxes(showticklabels=False)
    fig.update_yaxes(showticklabels=False)
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
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
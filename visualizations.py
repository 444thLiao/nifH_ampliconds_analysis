"""
Parse the amplicons visualization 
"""
from plotly.subplots import make_subplots
from scipy.stats import mannwhitneyu
from os.path import *
from Bio import SeqIO
import pandas as pd
from subprocess import check_call
import json
from glob import glob
from tqdm import tqdm
import os
import warnings
import plotly.express as px
from plotly import graph_objects as go

from for_software.for_EPA.parse_jplace import *

warnings.filterwarnings("ignore")
os.chdir('/mnt/ivy/thliao/project/nif_jjtao/amplicons_study')


pd.set_option('display.max_rows', 1001)

t2color = {"Sym": "#f4a7a6",
           "Sym_Basal": "#f4a7a6",
           "FL": "#71dab8",
           "PB": "#60b0fd",
           "Cluster FL": "#71dab8",
           "Cluster PB": "#60b0fd",
           }


## parse functions

def get_num(tree, nodes, node_name2edge_num):
    t = tree.get_common_ancestor(nodes)
    return [node_name2edge_num[_.name] for _ in t.traverse()]


def classification_criteria(tree,node_name2edge_num):
    " for wang2020"
    cluster_PB = get_num(tree,['Bradyrhizobium_oligotrophicum_S58','Bradyrhizobium_sp_BTAi1'],node_name2edge_num)
    cluster_PB.append('Bradyrhizobium_sp_STM_3843')
    
    
    FL = get_num(tree,['Bradyrhizobium_iriomotense_SZCCT0346','Bradyrhizobium_sp_BK707'],node_name2edge_num)
    _tmp = get_num(tree,['Bradyrhizobium_sp_CCBAU_53426','Bradyrhizobium_sp_CCBAU_53424'],node_name2edge_num)
    for _ in _tmp:
        FL.remove(_)
#     FL.append(node_name2edge_num['Bradyrhizobium_AUGA_SZCCT0283'])  # FL near PB
    FL.extend(get_num(tree,['Bradyrhizobium_sp_S23321','Bradyrhizobium_sp_W'],node_name2edge_num))
    
    FL_in_Sym = []
    FL_in_Sym.append(node_name2edge_num['Bradyrhizobium_mercantei_SEMIA_6399'])
    FL_in_Sym.append(node_name2edge_num['Bradyrhizobium_yuanmingense_P10_130'])
    FL_in_Sym.append(node_name2edge_num['Bradyrhizobium_liaoningense_CCBAU_83689'])
    
    sym_basal = get_num(tree,['Bradyrhizobium_sp_Cp5_3','Bradyrhizobium_stylosanthis_BR_446'],node_name2edge_num)
    sym_basal.extend(_tmp)
    
    Sym_all = []
    Sym_all.extend(get_num(tree,['Bradyrhizobium_sp_WSM4349','Bradyrhizobium_canariense_UBMAN05'],node_name2edge_num))
    Sym_all.extend(get_num(tree,['Bradyrhizobium_sp_WSM2254','Bradyrhizobium_pachyrhizi_BR3262'],node_name2edge_num))
    
    for _ in FL_in_Sym:
        Sym_all.remove(_)
    
    criteria = {"PB":cluster_PB,
                "FL":FL, #+FL_in_Sym,
                "Sym_Basal":sym_basal,
                "Sym":Sym_all}
    return criteria


def get_classification(tree, node_name2edge_num, df):
    criteria = classification_criteria(tree, node_name2edge_num)

    seq2type = {}
    for type_name, edge_nums in criteria.items():
        sub_df = df.loc[df['edge_num'].isin(edge_nums)]
        names = list(sub_df['name'])
        seq2type.update({n: type_name for n in names})
    return seq2type


def get_closest_matches(tree, node_name2edge_num, df):
    edge_num2node_name = {v: k for k, v in node_name2edge_num.items()}
    seq2edge_num = dict(zip(df['name'],
                            [edge_num2node_name[_] for _ in df['edge_num']]))
    if not tree.name:
        tree.name = 'Root'
    name2node = {n.name: n for n in tree.traverse()}
    leafs = tree.get_leaf_names()
    zero_seq = [k for k, v in seq2edge_num.items() if v == '0']
    seq2edge_num = {k: v for k, v in seq2edge_num.items() if v != '0'}
    seq2closest_match = {k: v
                         if v in leafs else name2node[v].get_closest_leaf()
                         for k, v in seq2edge_num.items()}
    return seq2closest_match


# parse jplace results
guppy_exe = "/home-user/thliao/download/pplacer-Linux-v1.1.alpha19/guppy "
classified_dict = {}
for jplace in tqdm(glob(f"./phylogenetic_placement/wang2020/epa_o/*/epa_result.jplace")):
    srr_id = basename(dirname(jplace))
    cmd = f"{guppy_exe} to_csv {jplace} > {dirname(jplace)}/jplace.csv"
    if not exists(f"{dirname(jplace)}/jplace.csv"):
        check_call(cmd,shell=1)
    df = pd.read_csv(f"{dirname(jplace)}/jplace.csv")
    df = df.sort_values('like_weight_ratio',ascending=False).groupby('name').head(1)
    obj = json.load(open(jplace))
    used_tree = obj["tree"]
    tree,node_name2edge_num = parse_tree_with_edges(used_tree)
    seq2type_epa = get_classification(tree,node_name2edge_num,df)
    # use blast to filter
    #seq2type_epa = {k:v for k,v in seq2type_epa.items() if k in normal_seqs[srr_id]}
    
    classified_dict[srr_id] = seq2type_epa
    # read in dada2 profilling table
    infa = f"./{srr_id}/{srr_id}_repr_dada2.fa"
    name2len = {_.name: len(_.seq)
                for _ in SeqIO.parse(infa, 'fasta')}
    count_tab = f"./{srr_id}/{srr_id}_table_dada2.csv"
    count_df = pd.read_csv(count_tab, sep='\t', index_col=0)
    if count_df.index[0] in [_.split(';')[0] for _ in name2len]:
        count_df = count_df
    else:
        count_df = count_df.T
    srr2total = count_df.sum(0)
    ratio_df = count_df/count_df.sum(0) * 100
    sub_ratio_df = ratio_df.loc[ratio_df.index.isin(seq2type_epa), :]
    sub_ratio_df.loc[:, 'type'] = [seq2type_epa[s] for s in sub_ratio_df.index]
    
    type2ratio = sub_ratio_df.groupby('type').sum().T
    tmp_df = pd.DataFrame(columns= ['PB', 'FL', 'Sym_Basal', 'Sym'])
    tmp_df = tmp_df.reindex(type2ratio.index)
    for c in type2ratio.columns:
        tmp_df.loc[:, c] = list(type2ratio.loc[:, c])
    tmp_df = tmp_df.fillna(0)
    tmp_df.loc[:, 'total reads'] = srr2total
    
    _odir = f'./wang2020/epa_out/{srr_id}'
    if not exists(_odir):
        os.system(f'mkdir -p {_odir}')
    tmp_df.to_csv(
       f'{_odir}/type_count.csv', sep='\t', index=1)

# parse profilling tables

df2 = pd.read_excel(
    '/mnt/ivy/thliao/project/nif_jjtao/amplicons_study/SraRunTable.xlsx')
df = pd.read_excel(
    '/mnt/ivy/thliao/project/nif_jjtao/amplicons_study/nifH metagenome_filtered_16S.xlsx')

df2.loc[df2['Experiment'].isin(set(df['Experiment Accession']))]

removed_srr = {'SRP135966': ['SRR6854367'],
               "ERP108575": ["ERR2560196"]}

run2info_df = df2.set_index('Run')
useful_columns = ['BioProject', 'BioSample', "Experiment",
                  'Library Name', "Organism", "Sample Name", "SRA Study",
                  "AvgSpotLen", "Bases", 'Collection_Date',
                  'BioSampleModel',
                  'geo_loc_name_country',
                  'geo_loc_name_country_continent',
                  'geo_loc_name', "Isolation_source", 'env_biome',
                  'env_feature',
                  'env_material',
                  ]

all_dfs = []
for tab in glob(f'./wang2020/epa_out/*/type_count.csv'):
    _df = pd.read_csv(tab, sep='\t', index_col=0)
    _df.index = [_.split('_')[0] for _ in _df.index]
    all_dfs.append(_df)
_type_count_df = pd.concat(all_dfs, 0)
_type_count_df.columns = [_ + ' (%)'
                          for _ in _type_count_df.columns[:-1]] + ['total number of remained reads']
no_zero_sra = _type_count_df.index[_type_count_df.iloc[:, :-1].sum(1) != 0]
sub_run2info_df = run2info_df.reindex(columns=list(
    _type_count_df.columns) + list(run2info_df.columns))
sub_run2info_df.loc[_type_count_df.index,
                    _type_count_df.columns] = _type_count_df
sub_srainfo_df = sub_run2info_df.reindex(no_zero_sra)


classified_critera = {"soil": ['rhizosphere metagenome',
                               'soil metagenome',
                               'mine tailings metagenome',
                               'mine metagenome',
                               'sediment metagenome',
                               'soil crust metagenome',
                               'compost metagenome',
                               'peat metagenome',
                               'oil sands metagenome'
                               ],
                      "marine": ['marine metagenome',
                                 'coral metagenome',
                                 'seawater metagenome',
                                 "marine sediment metagenome",
                                 'marine plankton metagenome',
                                 'beach sand metagenome',
                                 'coral reef metagenome',
                                 'mangrove metagenome'
                                 ],
                      "plant": ['plant metagenome',
                                'root metagenome',
                                'phyllosphere metagenome',
                                'Matricaria chamomilla',
                                'Solanum sp. 0048',
                                'Calendula officinalis'
                                ],
                      'others': ['bioreactor metagenome',
                                 'wastewater metagenome',
                                 'Bacteria',
                                 'metagenome',
                                 'synthetic metagenome',
                                 'biofilm metagenome',
                                 'sludge metagenome',
                                 'freshwater metagenome',
                                 'mollusc metagenome'
                                 ]}

classified_critera_rev = {v: k for k,
                          _l in classified_critera.items() for v in _l}
srr2env = {}
sra2env = {}
for srr, row in sub_srainfo_df.iterrows():
    srr2env[srr] = classified_critera_rev.get(row['Organism'])
    sra2env[row['SRA Study']] = classified_critera_rev.get(row['Organism'])

remained_columns = ['Relative ratio PB (%)', 'Relative ratio FL (%)',  'Relative ratio Sym (%)', 'Relative ratio SymBasal (%)',
                    'PB (%)', 'FL (%)',  'Sym (%)', 'Sym_Basal (%)',
                    "Total Number (read)",
                    "Environmental type", "BioProject", "BioSample", "SRA Study", 'Organism', 'Isolation_source', "lat_lon", "Center Name", ]

ssub_srainfo_df = sub_srainfo_df.reindex(
    columns=['PB (%)', 'FL (%)', 'Sym (%)', 'Sym_Basal (%)'])
ssub_srainfo_df = ssub_srainfo_df.div(ssub_srainfo_df.sum(1), 0) * 100

final_df = sub_srainfo_df.reindex(columns=remained_columns)
final_df.loc[:, "Environmental type"] = [srr2env[_] for _ in final_df.index]
final_df.loc[:,
             "Total Number (read)"] = sub_srainfo_df['total number of remained reads']

final_df.loc[:, "Relative ratio PB (%)"] = ssub_srainfo_df['PB (%)'].round(2)
final_df.loc[:, "Relative ratio FL (%)"] = ssub_srainfo_df['FL (%)'].round(2)
final_df.loc[:, "Relative ratio Sym (%)"] = ssub_srainfo_df['Sym (%)'].round(2)
final_df.loc[:,
             "Relative ratio SymBasal (%)"] = ssub_srainfo_df['Sym_Basal (%)'].round(2)
final_df.to_excel('./wang2020/epa_out/SraRun_metadata.xlsx',
                  index=True, index_label='SRR ID')


ssub_srainfo_df = sub_srainfo_df.reindex(
    columns=['PB (%)', 'FL (%)', 'Sym (%)', 'Sym_Basal (%)'])
ssub_srainfo_df = ssub_srainfo_df.div(ssub_srainfo_df.sum(1), 0) * 100
new_df = []
for c in ['PB (%)', 'FL (%)', 'Sym (%)', 'Sym_Basal (%)']:
    _sub_df = pd.DataFrame(ssub_srainfo_df[c])
    _sub_df.columns = ['ratio']
    _sub_df.loc[:, 'nifH type'] = c.replace(' (%)', '')
    new_df.append(_sub_df)
new_df = pd.concat(new_df, axis=0)
new_df.loc[:, 'env type'] = [srr2env[_] for _ in new_df.index]
new_df.to_csv('./wang2020/fig5_data.csv', index=1, sep='\t')


# visualization  (violin plot)
def get_stars(p):
    if p <= 0.05 and p > 0.01:
        return '*'
    elif p <= 0.01 and p > 0.001:
        return '**'
    elif p <= 0.001:
        return '***'
    return 'ns'


draw_df = new_df.replace('FL', "Cluster FL")
draw_df = draw_df.replace('PB', "Cluster PB")
draw_df = draw_df.fillna(0)

fig = make_subplots(rows=1, cols=4, shared_yaxes=True, horizontal_spacing=0.01,
                    subplot_titles=['marine (n)', 'soil', 'others', 'plant']
                    )
for i, env in enumerate(['marine', 'soil', 'others', 'plant']):
    traces = []
    sub_df = draw_df.loc[draw_df['env type'] == env, :]
    for nifH_t, idxs in sub_df.groupby('nifH type').groups.items():
        _sub_df = sub_df.loc[sub_df['nifH type'] == nifH_t, :]
        traces.append(go.Violin(x=[nifH_t]*len(idxs),
                                y=_sub_df.loc[idxs, 'ratio'],
                                legendgroup=nifH_t,
                                name=nifH_t,
                                line_color=t2color[nifH_t],
                                spanmode='hard',
                                width=0.7,
                                points=False,
                                meanline=dict(visible=True),
                                showlegend=False
                                ))
    fig.add_traces(traces, 1, i+1)

for idx, env in enumerate(['marine', 'soil', 'others', 'plant']):
    pos = [1.5, 0.5, 1, 2, 3, 4, 5]
    base_height = 110
    per_height = 5
    _count = 0
    for f1, f2 in [('Cluster PB', 'Sym'),
                   ('Cluster PB', 'Cluster FL'),
                   ('Cluster FL', 'Sym'),
                   ]:
        r1 = draw_df.loc[(draw_df['env type'] == env) &
                         (draw_df['nifH type'] == f1), 'ratio']
        r2 = draw_df.loc[(draw_df['env type'] == env) &
                         (draw_df['nifH type'] == f2), 'ratio']
        t = mannwhitneyu(r1, r2,
                         )
        p = t.pvalue

        fig.add_traces(go.Scatter(x=[f1, f1, f2, f2],
                                  y=[base_height,
                                     base_height+per_height,
                                     base_height+per_height,
                                     base_height],
                                  showlegend=False,
                                  mode='lines', line=dict(color='#000000')), 1, idx+1)
    #     fig.add_annotation(
    #         x=pos[_count],
    #         y=base_height+12,
    #         xanchor='center',
    #         xref=f'x{idx+1}',
    #         font=dict(size=15),
    #         showarrow=False,

    #         text="<Br> {:.2e} ".format(p),
    #     )
        fig.add_annotation(
            x=pos[_count],
            y=base_height+10,
            xanchor='center',
            xref=f'x{idx+1}',
            font=dict(size=15),
            showarrow=False,
            text=get_stars(p),
        )
        #print(get_stars(p),p,env,f1,f2)
        base_height += 10
        _count += 1

fig.layout.template = "simple_white"
fig.layout.width = 1400
fig.layout.height = 500
fig.layout.font.size = 20
#fig.layout.yaxis.title = "Relative abundance of <Br> different types of nif clusters"
for _ in fig.layout.annotations[:4]:
    _['font']['size'] = 23
fig.layout.yaxis.tickvals = [0, 20, 40, 60, 80, 100]
fig.layout.yaxis.title = "Relative abundance (%)"
# fig.show()


# visualization  (geological distributions)

# draw global distribution
def convert_lat_lon(lat_lon):
    infos = lat_lon.split(' ')

    lat = float(infos[0]) if infos[1] == 'N' else -float(infos[0])
    lon = float(infos[2]) if infos[3] == 'E' else -float(infos[2])
    return lat, lon


printable_lat_lon_df = sub_srainfo_df[(~sub_srainfo_df['lat_lon'].isna()) & (
    ~sub_srainfo_df['lat_lon'].isin(['Not Applicable', 'missing', 'not applicable', 'not collected']))]

printable_lat_lon_df.loc[:, 'lat'] = [convert_lat_lon(
    _)[0] for _ in printable_lat_lon_df['lat_lon']]
printable_lat_lon_df.loc[:, 'lon'] = [convert_lat_lon(
    _)[1] for _ in printable_lat_lon_df['lat_lon']]

printable_lat_lon_df.loc[:, 'env_type'] = [srr2env[_]
                                           for _ in printable_lat_lon_df.index]


for t_nifH in ['PB (%)', 'FL (%)', 'Sym (%)', 'Sym_Basal (%)']:
    fig = px.scatter_geo(printable_lat_lon_df.loc[printable_lat_lon_df[t_nifH] != 0],
                         lat='lat',
                         lon='lon',
                         hover_data=['PB (%)',
                                     'FL (%)',
                                     'Sym_Basal (%)',
                                     'Sym (%)',
                                     'env_type'],
                         color='env_type',
                         color_discrete_map={'marine': '#636efa',
                                             'soil': '#ef583e',
                                             'others': '#ac65fa',
                                             'plant': '#01cc96'},
                         #size="FLnif (%)",
                         projection="natural earth")

    fig.layout.width = 1200
    fig.layout.height = 800
    fig.layout.font.size = 30
    fig.layout.title = f"{t_nifH.split(' ')[0]} nifH"
    # fig.show()

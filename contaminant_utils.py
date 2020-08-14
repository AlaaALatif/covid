import re  # regex
import numpy as np  # linear algebra
import pandas as pd  # data processing
from path import Path  # file I/O
import plotly  # viz
import plotly.figure_factory as ff  # viz
import plotly.express as px  # viz
import plotly.graph_objects as go  #viz
from jinja2 import Environment, FileSystemLoader  # html template engine




def generate_html(general_hmap, cont_hmap, cont_table, num_samples, num_conts, expt_name):
    # dir containing our template
    file_loader = FileSystemLoader('templates')
    # load the environment
    env = Environment(loader=file_loader)
    # load the template
    template = env.get_template('contamination.html')
    # render data in our template format
    html_output = template.render(general_hmap=general_hmap, cont_hmap=cont_hmap, 
                                  cont_table=cont_table, num_samples=num_samples,
                                  num_conts=num_conts, expt_name=expt_name)
    return html_output


def save_html(html_output: str, filename: str):
    with open(filename, 'w') as f:
        f.write(html_output)
        
        
def generate_heatmap(data, x, y):
    heatmap = go.Heatmap(z=data.T, x=x, y=y)
    plot = [heatmap]
    fig = go.Figure(data = plot)
    return plotly.offline.plot(fig, include_plotlyjs=False, output_type='div')



def generate_table(ans, contaminated_samples):
    return (ans.loc[(ans['sample'].isin(contaminated_samples)) & (ans['paired_read_count'] > 0)]
                 .set_index(['sample', 'paired_read'])
                 .sort_index()[['paired_read_count']])


def load_all_data(input_pths: list):
    """Load barcode read data from the given list of sample file paths"""
    ans = pd.DataFrame()
    for sample_pth in input_pths:
        try:
            df = prepare_data(*load_data(sample_pth))
        except: continue
        df = (df.reset_index()#[['forward_barcode', 'reverse_barcode']]
                .drop_duplicates())
        df['sample'] = sample_pth.basename().split('_')[0]
        ans = pd.concat([ans, df], axis=0)
    return ans


def load_data(input_file: Path):
    with open(input_file, 'r') as f:
        input_data = f.readlines()
    header = re.split('\t|\n', input_data[0])[:-1]
    data = []
    for row in input_data[1:]:
        data.append(re.split('\t|\n', row)[:-1])
    df = pd.DataFrame(data=data, columns=header)
    return df, header


def prepare_data(df: pd.DataFrame, header: list):
    # convert read counts to integer
    df[header[-1]] = df[header[-1]].astype(int)
    # consolidate barcodes and their reverse complements 
    df = df.apply(_consolidate_reverse_complements, axis=1)
    # merge barcodes and their rcs and sum their read counts
    df = (df.groupby(['forward_barcode', 'reverse_barcode'])
            .agg({'paired_read_count': 'sum'})
            .reset_index())
    # create df with all possible paired combos of forward and reverse barcodes 
    forward_bcodes = df[['forward_barcode']].drop_duplicates()
    reverse_bcodes = df[['reverse_barcode']].drop_duplicates()
    forward_bcodes['key'] = 0
    reverse_bcodes['key'] = 0
    all_pairs = forward_bcodes.merge(reverse_bcodes, how='outer', on='key').drop(columns='key')
    all_pairs = (all_pairs.merge(df, on=['forward_barcode', 'reverse_barcode'], how='left')
                          .fillna(0)
                          .set_index(['forward_barcode', 'reverse_barcode']))
    return all_pairs

def _consolidate_reverse_complements(x: str):
    x['forward_barcode'] = x['forward_barcode'].split('_')[0]
    x['reverse_barcode'] = x['reverse_barcode'].split('_')[0]
    return x


def get_unique_barcodes(x):
    x = set(x.unique())
    x.discard('unknown')
    return x


def is_contaminant(x):
    if len(x['uniq_forward_bcodes']) > 1 or len(x['uniq_reverse_bcodes']) > 1:
        return True
    return False
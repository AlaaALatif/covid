import sys
import numpy as np
from contaminant_utils import *
from path import Path



if __name__ == "__main__":
    # grab user args
    out_pth = Path(snakemake.output)
    expt_pth = out_pth.parent
    input_pths = [Path(p) for p in snakemake.input.barcodes]
#     print(input_pths)
    # print(input_pths)
    num_samples = len(input_pths)
    # grab data from each sample
    ans = load_all_data(input_pths)
    # generate paired reads
    ans['paired_read'] = ans.apply(lambda x: x['forward_barcode'] + '-' + x['reverse_barcode'], axis=1)
    # compute log of read counts 
    ans['log_count'] = ans['paired_read_count'].apply(lambda x: np.log(x+1))
    # generate heatmap matrix of (logged: optional) read counts per sample per paired read
    hmap, data, x, y = get_heatmap_data(ans)
    general_hmap = generate_heatmap(data, x, y)
    # generate heatmap matrix of only samples suspected of contamination
    cont_data, cont_x, cont_y = get_contaminated_data(data, hmap)
    cont_hmap = generate_heatmap(cont_data, cont_x, cont_y)
    # table of barcode read counts for contaminated samples only
    cont_table = generate_table(ans, cont_x)
    # number of contaminated samples
    num_conts = len(cont_x)
    # generate scatter plot of coverage vs mapped reads ratio
    mapped_pth = Path(snakemake.input.mapping)
    coverage_pth = Path(snakemake.input.coverage)
    mapped_df = pd.read_csv(mapped_pth, delimiter='\t')
    mapped_df['mapped ratio'] = mapped_df['mapped'] / (mapped_df['unmapped']+mapped_df['mapped'])
    mapped_df['SAMPLE'] = mapped_df['SAMPLE'].apply(lambda x: x.split('/')[-1].split('_')[0])
    coverage_df = pd.read_csv(coverage_pth, delimiter='\t')
    coverage_df['SAMPLE'] = coverage_df['SAMPLE'].apply(lambda x: x.split('_')[0])
    ans = pd.merge(mapped_df, coverage_df, on='SAMPLE', how='inner')[['SAMPLE', 'COVERAGE', 'mapped ratio']]
    scattr_plot = go.Figure(data=go.Scatter(x=ans['mapped ratio'],
                                y=ans['COVERAGE'],
                                mode='markers',
                                text=ans['SAMPLE'])) # hover text goes here

    scattr_plot.update_layout(title='Coverage versus Mapped Reads Ratio',
                  xaxis_title="Ratio of Mapped Reads",
                  yaxis_title="Coverage",
                  template='plotly',
                  height=800)
    # generate html string
    html_output = generate_html(general_hmap, cont_hmap, cont_table, scattr_plot,
                                num_samples, num_conts, out_pth)
    # save report to file
    save_html(html_output, out_pth)
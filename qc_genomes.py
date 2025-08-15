"""Get QC information from various outputs."""
import pandas as pd
from tqdm import tqdm
from os import path as op
from functools import partial
import matplotlib.pyplot as plt
from collections import defaultdict, Counter
from matplotlib.backends.backend_pdf import PdfPages
pbar = partial(tqdm, bar_format='{l_bar}{bar:15}{r_bar}')

perc_keys = {
    # BUSCO and compleasm
    'C' : 'Complete',
    'S' : 'Complete, single-copy',
    'D' : 'Complete, duplicated',
    'F' : 'Fragmented',
    'M' : 'Missing',
    'n' : 'total BUSCOs',
    'E' : 'Essentially missing',
    # compleasm
    'I' : 'Intermediate',
    'N' : 'total'
}

def read(file: str, lines=True, ignore_blank=False):
    """Read lines from a file.

    Return a list of lines, or one large string
    """
    with open(file, "r") as o:
        text = o.read()

    if lines is True:
        text = text.split("\n")
        if ignore_blank is True:
            text = [line for line in text if line != ""]

    return text


def busco(short_summary, busco_data=None, name=None):
    """Parse the short_summary file output from BUSCO.

    Parameters
    ----------
    short_summary : Path
    name : str | int
        name for series
    busco_data : pd.Series
        series to add busco information; created if None

    Returns
    -------
    busco_data : pd.Series
        series with BUSCO information 
            index = ['Complete', 'Complete, single-copy', 'Complete, duplicated', 'Fragmented', 'Missing', 'total BUSCOs']
    """
    if busco_data is None:
        busco_data = pd.Series()

    busco_data.name = name

    text = read(short_summary, lines=True, ignore_blank=True)

    for i, line in enumerate(text):
        if 'Results' in line:
            results_data = text[i+1 : i+8]

            clean_percs = results_data[0].split()[0].replace(']', '').replace('[',',').split(',')

            perc_results = {}
            for key_val in clean_percs:
                key, val = key_val.split(":")
                perc_results[perc_keys[key]] = val

            for line in results_data[1:]:
                num, key, *_ = line.split('\t')[1:]

                if key.startswith('Total BUSCO'):
                    key_letter = 'n'
                    entry = f'{num} (100%)'
                else:
                    key_letter = key[-2:-1]
                    entry = f'{int(num):,} (%s)' % perc_results[perc_keys[key_letter]]

                busco_data.loc[perc_keys[key_letter]] = entry

    return busco_data


def compleasm(summary_file, compleasm_data=None, name=None):
    """Parse the summary_file output from compleasm.

    Parameters
    ----------
    summary_file : Path
    name : str | int
        name for series
    compleasm_data : pd.Series
        series to add compleasm information; created if None

    Returns
    -------
    compleasm_data : pd.Series
        series with compleasm information 
            index = ['Complete', 'Complete, single-copy', 'Complete, duplicated', 'Fragmented', 'Missing', 'total BUSCOs']
    """
    if compleasm_data is None:
        compleasm_data = pd.Series()

    compleasm_data.name = name

    text = read(summary_file, lines=True, ignore_blank=True)

    total_perc = 0
    total_count = 0
    for line in text:
        # print(line)
        if line.startswith("## lineage"):
            continue
        if line.startswith('N:'):
            key, count = line.split(":")
            perc = '100%'
        else:
            key, perc, count = line.replace(',', '').replace(':', '\t').split()

        if key in ['S', 'D']:
            total_perc = round(total_perc + float(perc.removesuffix("%")), 2)
            total_count += int(count)

        # print(count, perc)
        compleasm_data.loc[perc_keys[key]] = f'{int(count):,} ({perc})'

    idx = [perc_keys['C']] + compleasm_data.index.tolist()

    compleasm_data.loc[perc_keys['C']] = f'{int(total_count):,} ({total_perc}%)'

    return compleasm_data.loc[idx]


def quast(report, quast_data=None, name=None):
    """Parse the report.tsv output from quast.

    Parameters
    ----------
    report : Path
    quast_data : pd.Series
        dataframe to add quast information; created if None
    name : str | int
        name for series

    Returns
    -------
    quast_df : pd.DataFrame
        dataframe with quast information 
            index = ['# contigs', 'Largest contig', 'Total length', 'GC (%)', 'N50', 'N90', 'auN', 'L50', 'L90', "# N's per 100 kbp"]
    """
    if quast_data is None:
        quast_data = pd.Series()

    quast_data.name = name

    text = read(report, lines=True, ignore_blank=True)

    begin = False
    for line in text:
        if line.startswith('# contigs\t'):
            # print('beginning', line)
            begin = True

        if begin is True:
            key, val = line.split('\t')
            if '.' in val:
                val = f"{float(val):,}"
            else:
                val = f"{int(val):,}"
            quast_data.loc[key] = val

    return quast_data


def nanoplot_stats(statfiles, genome_size=None, verbose=True):
    """Parse NanoStats file from NanoPlot output.
    
    Parameters
    ----------
    statsfiles : str | list
        the NanoPlot output files that end with `_NanoStats.txt`. Note underscore is assumed.
    genome_size : int | float
        the size of the genome in base-pairs; to be used to calculate coverage - see Notes

    Returns
    -------
    pd.DataFrame
        index = nano plot metrics, columns = names inferred from basenames of `statfiles`

    Notes
    -----
    coverage is calculated across all files provide (assumed to go towards a single assembly)
    """
    total_reads = 0
    total_bases = 0
    read_info = defaultdict(Counter)
    read_qual = Counter()
    read_qual_size = Counter()
    
    if isinstance(statfiles, str):
        statfiles = [statfiles]

    dfs = []
    for statfile in statfiles:
        name = op.basename(statfile).split('_NanoStats.txt')[0]
        
        df = pd.read_table(statfile)
        data = df.set_index('Metrics').to_dict()['dataset']
    
        total_reads += float(data['number_of_reads'])
        total_bases += float(data['number_of_bases'])

        try:
            for qual in [10, 15, 20, 25, 30]:
                num, perc, mb = data[f'Reads >Q{qual}:'].split()
                assert mb.endswith('Mb')
                read_info[f'Reads >Q{qual}:']['total_read_length'] += float(num)
                read_info[f'Reads >Q{qual}:']['total_read_size (Mb)'] += float(mb.rstrip('Mb'))
            print_read_info = True
        except KeyError as e:
            print_read_info = False
            pass
    
        df.reset_index(drop=False)
        df.dataset = df.dataset.apply(lambda entry: "%s %s" % (f'{float(entry.split()[0]):,}', ' '.join(entry.split()[1:])))
        if verbose is True:
            print(f'\n**{name} (n50 = %s)**\n' % data['n50'])
            print(df.to_markdown())
        df.columns = ['Metrics', name]
        dfs.append(df.set_index('Metrics'))

    if print_read_info is True and verbose is True:
        print('\n**read info**\n')
        print(pd.DataFrame(read_info).to_markdown())

    if verbose is True:
        print(f'\n{total_reads = :,}')
        print(f'{total_bases = :,}')

    if genome_size is not None and verbose is True:
        print(f'coverage = ', round(total_bases / genome_size, 2))

    return pd.concat(dfs, axis=1)


def nanoplot_traceplots(statfiles, statnames=None, dataset_name=None, genome_size=None, pdf_path=None, verbose=True):
    """Plot data output from NanoPlot output, assuming sequential QC order specified in `statfiles`.

    Parameters
    ----------
    statfiles : list
        a list of the paths to NanoStats.txt files. The order of the list is the assumed to be the
        order the dataset went through sequential QC steps - eg first lenght-filtering, then 
        rebasecalling, centrifuge, etc
    statnames : list
        names for each file in statfiles - used to label x-axis labels in figure
    dataset_name : str
        the name of the dataset pertaining to the statfiles - for labeling the figure
    genome_size : int
        if provided, `nanoplot_traceplots` will calculage and plot coverage changes
    pdf_path : str
        if provide, the path where the figure will be saved
    verbose : bool
        if True, markdown-formatted tables from each of the `statfiles` will be printed along with
        other information

    Returns
    -------
    stats : pandas.DataFrame
        index = metrics, cols = statnames (or if statnames not provided, Step numbers are assigned)
    """
    
    metrics = ['number_of_reads', 'number_of_bases', 'median_read_length', 'read_length_stdev', 'n50']

    if statnames is None:
        statnames = [f'Step {x}' for x in list(range(len(statfiles)))]

    # create a dataframe with `metrics` as index and `statnames` as columns
    stats = pd.DataFrame(index=metrics)
    for name, statfile in zip(pbar(statnames), statfiles):
        if verbose is True:
            print(f'\033[4m\033[1m{name}\033[0m')  # bold underline
            print(' ', statfile)

        df = nanoplot_stats(statfile, verbose=verbose)

        # add each metric into dataframe where col is one of the `statnames` and index is the `metric`
        for metric in metrics:
            stats.loc[metric, name] = df.loc[metric, df.columns[0]].replace(',', '')

    stats = stats.astype(float).astype(int)

    # add in percentages of the original (first step or stepfile)
    for i, perc_metric in enumerate(['number_of_reads', 'number_of_bases', 'median_read_length']):
        stats.loc[f'percent of original\n{perc_metric}'] = stats.loc[perc_metric] / stats.loc[perc_metric].iloc[0]
        metrics.insert(i+3, f'percent of original\n{perc_metric}')

    if genome_size is not None:
        stats.loc['coverage'] = stats.loc['number_of_bases'] / genome_size
        metrics.append('coverage')

    # create lineplot figure with values of metrics for y-axis and `statnames` for x-axis labels
    fig, axes = plt.subplots(ncols=3, nrows=3, sharex=False, sharey=False, figsize=(15, 12))
    for metric, ax in zip(metrics, axes.flat):
        ax.plot(stats.loc[metric], marker='o', linestyle='-')
        ax.set_title(metric)

    fig.suptitle(dataset_name, fontsize=18)

    if len(stats.columns) > 4:
        for ax in axes.flat:
            plt.setp(ax.get_xticklabels(), rotation=45, ha='center')

    # if genome_size not provided, remove empty figure
    if genome_size is None:
        last_ax = axes.flat[-1]
        for spine in last_ax.spines.values():
            spine.set_visible(False)
        
        last_ax.set_xticks([])
        last_ax.set_yticks([])
        last_ax.set_xlabel('')
        last_ax.set_ylabel('')
        last_ax.set_title('')
    
    plt.tight_layout()

    # save
    if pdf_path is not None:
        with PdfPages(pdf_path) as pdf:
            pdf.savefig(bbox_inches="tight")
        print(f'\033[1mSaved to\033[0m: ', pdf_path)  # bold
        
    plt.show()
    
    return stats

if __name__ == '__main__':
    pass

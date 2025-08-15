help documentation as of 

commit 5b9e456e4d99f6b1e126b4bd53d940335dc0e795  
Author: Brandon Lind <lind.brandon.m@gmail.com>  
Date:   Fri Aug 15 10:45:32 2025 -0400

----
### Python Library Documentation: module qc_genomes
```

NAME
    qc_genomes - Get QC information from various outputs.

FUNCTIONS
    busco(short_summary, busco_data=None, name=None)
        Parse the short_summary file output from BUSCO.

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

    compleasm(summary_file, compleasm_data=None, name=None)
        Parse the summary_file output from compleasm.

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

    nanoplot_stats(statfiles, genome_size=None, verbose=True)
        Parse NanoStats file from NanoPlot output.

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

    nanoplot_traceplots(statfiles, statnames=None, dataset_name=None, genome_size=None, pdf_path=None, verbose=True)
        Plot data output from NanoPlot output, assuming sequential QC order specified in `statfiles`.

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

    quast(report, quast_data=None, name=None)
        Parse the report.tsv output from quast.

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

    read(file: str, lines=True, ignore_blank=False)
        Read lines from a file.

        Return a list of lines, or one large string

DATA
    pbar = functools.partial(<class 'tqdm.std.tqdm'>, bar_format='{l_bar}{...
    perc_keys = {'C': 'Complete', 'D': 'Complete, duplicated', 'E': 'Essen...

```


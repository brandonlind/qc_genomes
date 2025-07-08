help documentation as of 

commit f2afd8ccfb70d3b698599e7a6782dacb1e1833d7  
Author: Brandon Lind <lind.brandon.m@gmail.com>  
Date:   Tue Jul 8 18:09:52 2025 -0400

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

    nanoplot_stats(statfiles, genome_size=None)
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
    perc_keys = {'C': 'Complete', 'D': 'Complete, duplicated', 'E': 'Essen...

```


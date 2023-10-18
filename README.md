# ptmanalysis
PTM analysis scripts

## Dependencies

```
# Download plasmodb fasta and gffs into data dir
wget -r -np -nH --cut-dirs=3 --mirror --convert-links --execute="robots = off" -R "index.html*" -P data https://plasmodb.org/common/downloads/release-64/Pfalciparum3D7/
```

## Usage
```
# writes coverage plots to coverage_plots, and m6a_genes_modified.csv to data dir
Rscript PlotSitesOnGenome.R 

# writes volcano plots to volcano_plots
Rscript PlotModifiedGenesVolcano.R
```
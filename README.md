# ptmanalysis
PTM analysis scripts

Dependencies:

```
# Download plasmodb fasta and gffs into data dir
wget -r -np -nH --cut-dirs=3 --mirror --convert-links --execute="robots = off" -R "index.html*" -P data https://plasmodb.org/common/downloads/release-64/Pfalciparum3D7/
```

Usage:
```
Rscript PlotSitesOnGenome.R # writes coverage plots to coverage_plots
Rscript PlotModifiedGenesVolcano.R # writes volcano plots to volcano_plots
```
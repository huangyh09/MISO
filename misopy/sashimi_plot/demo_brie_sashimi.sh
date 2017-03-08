#!/bin/sh


### generate gff annotations
index_gff --index SE.gold.gff3 event-data/

plot_dir=plots/
gff_dir=event-data/

sashimi=sashimi_plot.py

python $sashimi --plot-event "ENSMUSG00000030689" $gff_dir settings/sashimi_settings.txt --output-dir $plot_dir  --plot-label "test.pdf" #--plot-title "DNMT3B-exon2" --no-posteriors

#!/bin/sh

### generate gff annotations
# index_gff --index test-data/events.gff test-data/event-data/

### create sashimi plot
python sashimi_plot.py --plot-event "chr17:45816186:45816265:-@chr17:45815912:45815950:-@chr17:45814875:45814965:-" test-data/event-data/ settings/sashimi_plot_settings.txt --output-dir test-plot

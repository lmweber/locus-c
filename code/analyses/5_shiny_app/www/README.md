# Locus coeruleus

**Shiny app to explore locus coeruleus (LC) dataset**

For interactive annotation:
- use either `Clusters (interactive)` or `Gene (interactive)` tab
- select sample and genes (`Continuous variable to plot`) in dropdown menus on left
- select alpha transparency (`Spot transparency level`), e.g. value 0.25, 0.33, or 0.5
- select spot size to approximately match the Visium spot resolution (55um spot diameter, 100um spot center to center) (see below for suggested spot sizes for each zoom level)
- set `Minimum count value` to -1 to show zero-expression spots if needed
- select `Gene color scale`, e.g. viridis or magma
- to annotate spots, enter a name (e.g. `LC`, `manual_region`, `manual_spots`, etc) in the `Manual annotation label` field. Then either use the lasso tool or click `Enable spot manual annotations` and click individual spots, then click `Label selected points (from lasso)` to store the annotations in the running shiny instance.
- export annotations using the `Download manual annotations` button to export as .csv file. If NA / zero-expression spots need to be included (e.g. within larger lasso region), unselect the `Drop NA layer entries in the CSV file?` option.


Approximate spot sizes for each zoom level:
- zoom level 0: spot size 0.8
- zoom level 1: spot size 1.7
- zoom level 2: spot size 3.3
- zoom level 3: spot size 6.2
- zoom level 4: spot size 12.0


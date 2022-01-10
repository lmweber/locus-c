# Locus coeruleus

**Shiny app to explore locus coeruleus (LC) dataset**

For interactive annotation:
- use either `Clusters (interactive)` or `Gene (interactive)` tab
- select sample and genes (`Continuous variable to plot`) in dropdown menus on left
- select alpha transparency (`Spot transparency level`), e.g. value 0.25, 0.33, or 0.5
- select spot size to approximately match the Visium spot resolution (55um spot diameter, 100um spot center to center) (see below for suggested spot sizes for each zoom level)
- set `Minimum count value` to -1 to show zero-expression spots if needed, e.g. when using lasso selection tool
- select `Gene color scale`, e.g. viridis or magma
- to annotate spots, enter a name (e.g. `manual_lasso` or `manual_region`) in the `Manual annotation label` field, and select `Enable spot manual annotations`. Then click `Label selected points (from lasso)` to store annotations in the running shiny instance
- export annotations using `Download manual annotations` button to export as .csv file. Unselect `Drop NA layer entries in the CSV file?` to make sure NA / zero-expression spots are included (e.g. zero-expression spots within lasso region)


Approximate spot sizes for each zoom level:
- zoom level 0: spot size 0.8
- zoom level 1: spot size 1.7
- zoom level 2: spot size 3.3
- zoom level 3: spot size 6.2
- zoom level 4: spot size 12.0


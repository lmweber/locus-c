Spot-level data
===============

**Standard Operating Procedure (SOP) for Annotating Regions and Single Spots in Shiny app**

1. Go to 'spot-level data' tab
2. On the Documentation tab: select and set all the parameters (note setting these parameters on the documentation prevents reloading and is faster than setting them on 'clusters or gene(interactive)'):
    - Select the sample you want to work on using the 'Sample to Plot' drop down box
    - Set 'Discrete variable to plot' to all
    - Select the desired gene of interest or leave it on 'sum_umi for annotating the region and spots using the 'Continuous variable to plot' drop down box. For defining a region, using marker genes for that spatial domain can be helpful.
    - Set 'Gene scale' to logcounts
    - Set 'Image name' to hires
    - Set the spot transparency level to 0.25 (It can also be set to 0.33 or 0.5 based on the purpose and visibility). Sometimes if marker genes are not visible, try changing the spot transparency back to 1.
    - Approximate spot sizes for each zoom level:
          - zoom level 0: spot size 0.8
          - zoom level 1: spot size 1.7
          - zoom level 2: spot size 3.3
          - zoom level 3: spot size 6.2
          - zoom level 4: spot size 12.0
    - Set 'Minimum Count Value' to -1 to show all spots including those with zero expression of a given gene. This is useful for the lasso selection, where we want the whole region including zeros within the region
    - For 'Gene color scale' viridis and magma work quite well. Check or uncheck the box below the Gene color scale option to change the color order from darkest to lightest and vice versa
    - Save often, as the site can reload or disconnect automatically. To export the annotations, use the ‘Download manual annotations’ button to export to a .csv file. Unselect the Drop NA layer entries in the CSV file? button to make sure the NA / zero expression spots are included, since we want these included within the lasso region.
    - Name your .csv file with sample name and any specifics relevant to your annotation.
3. Go to 'Clusters(interactive)' or 'Gene(interactive)' tab and enter a name for your annotation under 'Manual annotation label'. E.g  Br6522_LC_round1_lasso
4. Hover over your image on the right hand side there are small icons that allow you to zoom in and out, lasso tool, pan etc.
5. Use the '+' and '-' icons to zoom in and out on the image. Use the '4 way arrow' icon to pan on the image from one point to another
6. To lasso a region, make sure to check the name in the Manual annotation label field and then click on 'lasso' icon (looks something like this: ථ) and lasso the desired region. Then click 'Label selected points (from lasso) with a manual annotation...' to store the labels annotation in the running shiny app. Now go to 'Discrete variable to plot' and select  Manual annotations and you will be able to see your lasso-ed region. At the top of the image, different clusters such as 'Lasso' and 'NA' are visible. You may have to uncheck 'NA' to see and confirm your Lasso region. Note, when lasso-ing a region the spots disappear and you will only see the histology image.
7. To annotate single spots, enter a name in the Manual annotation label field, and then check the 'Enable spot manual annotations by clicking on individual spots' box. Left click on the spot and it registers it automatically. At the top of your image, you will be able to see different clusters.
8. If you are using the 'gene (interactive)' tab for manual spot annotation, after marking a spot it refreshes so you lose the zoom level and the location where you annotated a spot, so it requires you to re navigate to the spot again using the icons. The advantage is that now you see a checkbox at the top with a separate cluster such as 'Spots'. The clusters are the names you entered in step 3 and 7. You can check or uncheck those clusters to see the annotations in each cluster and also note that a annotated spot can only be in one cluster.
9. Save the annotations using the Download manual annotations button.

Continuing from the downloaded annotations:

1. Go to 'Spot-level data' and on the Documentation tab, set all the setting as mentioned above
2. To overwrite or continue where you left at by uploading the .csv annotation file you saved using the 'Browse' button at the end of all the settings.
3. Enter the desired name in the 'Manual annotation label' and continue your annotations.


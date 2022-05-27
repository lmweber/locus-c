initial <- list()

################################################################################
# Settings for Reduced dimension plot 1
################################################################################

initial[["ReducedDimensionPlot1"]] <-
  new(
    "ReducedDimensionPlot",
    Type = "GLMPCA_MNN",
    XAxis = 1L,
    YAxis = 2L,
    FacetRowByColData = "Sample",
    FacetColumnByColData = "Sample",
    ColorByColumnData = "cellType.merged",
    ColorByFeatureNameAssay = "logcounts",
    ColorBySampleNameColor = "#FF0000",
    ShapeByColumnData = "Sample",
    SizeByColumnData = "sum",
    FacetRowBy = "None",
    FacetColumnBy = "None",
    ColorBy = "Column data",
    ColorByDefaultColor = "#000000",
    ColorByFeatureName = "SNAP25",
    ColorByFeatureSource = "---",
    ColorByFeatureDynamicSource = FALSE,
    ColorBySampleName = "3_AAACCCAAGTAGAATC-1",
    ColorBySampleSource = "---",
    ColorBySampleDynamicSource = FALSE,
    ShapeBy = "None",
    SizeBy = "None",
    SelectionAlpha = 0.1,
    ZoomData = numeric(0),
    BrushData = list(),
    VisualBoxOpen = FALSE,
    VisualChoices = "Color",
    ContourAdd = FALSE,
    ContourColor = "#0000FF",
    PointSize = 1,
    PointAlpha = 1,
    Downsample = FALSE,
    DownsampleResolution = 200,
    CustomLabels = FALSE,
    CustomLabelsText = "3_AAACCCAAGTAGAATC-1",
    FontSize = 1,
    LegendPointSize = 1,
    LegendPosition = "Bottom",
    HoverInfo = TRUE,
    LabelCenters = FALSE,
    LabelCentersBy = "Sample",
    LabelCentersColor = "#000000",
    VersionInfo = list(iSEE = structure(
      list(c(2L, 6L, 0L)), class = c("package_version", "numeric_version")
    )),
    PanelId = c(ReducedDimensionPlot = 1L),
    PanelHeight = 800L,
    PanelWidth = 6L,
    SelectionBoxOpen = FALSE,
    RowSelectionSource = "---",
    ColumnSelectionSource = "---",
    DataBoxOpen = FALSE,
    RowSelectionDynamicSource = FALSE,
    ColumnSelectionDynamicSource = FALSE,
    RowSelectionRestrict = FALSE,
    ColumnSelectionRestrict = FALSE,
    SelectionHistory = list()
  )

################################################################################
# Settings for Complex heatmap 1
################################################################################

initial[["ComplexHeatmapPlot1"]] <-
  new(
    "ComplexHeatmapPlot",
    Assay = "logcounts",
    CustomRows = TRUE,
    CustomRowsText = "LINC00682\nMET\nDBH\nCCNA1\nCYB561\nADRA2A\nPDGFD\nTH\nSLC6A2\nTINCR",
    ClusterRows = FALSE,
    ClusterRowsDistance = "spearman",
    ClusterRowsMethod = "ward.D2",
    DataBoxOpen = FALSE,
    VisualChoices = "Annotations",
    ColumnData = c("cellType.collapsed",
                   "donor"),
    RowData = character(0),
    CustomBounds = FALSE,
    LowerBound = 0L,
    UpperBound = 0L,
    AssayCenterRows = FALSE,
    AssayScaleRows = FALSE,
    DivergentColormap = "purple < black < yellow",
    ShowDimNames = "Rows",
    LegendPosition = "Right",
    LegendDirection = "Vertical",
    VisualBoxOpen = FALSE,
    NamesRowFontSize = 10,
    NamesColumnFontSize = 10,
    ShowColumnSelection = TRUE,
    OrderColumnSelection = TRUE,
    VersionInfo = list(iSEE = structure(
      list(c(2L, 6L, 0L)), class = c("package_version", "numeric_version")
    )),
    PanelId = c(ComplexHeatmapPlot = 1L),
    PanelHeight = 800L,
    PanelWidth = 6L,
    SelectionBoxOpen = FALSE,
    RowSelectionSource = "---",
    ColumnSelectionSource = "---",
    RowSelectionDynamicSource = FALSE,
    ColumnSelectionDynamicSource = FALSE,
    RowSelectionRestrict = FALSE,
    ColumnSelectionRestrict = FALSE,
    SelectionHistory = list()
  )

################################################################################
# Settings for Row data table 1
################################################################################

initial[["RowDataTable1"]] <-
  new(
    "RowDataTable",
    Selected = "TH",
    Search = "",
    SearchColumns = c("",
                      "", "", "", "", "", ""),
    HiddenColumns = character(0),
    VersionInfo = list(iSEE = structure(
      list(c(2L, 6L, 0L)), class = c("package_version",
                                     "numeric_version")
    )),
    PanelId = c(RowDataTable = 1L),
    PanelHeight = 600L,
    PanelWidth = 6L,
    SelectionBoxOpen = FALSE,
    RowSelectionSource = "---",
    ColumnSelectionSource = "---",
    DataBoxOpen = FALSE,
    RowSelectionDynamicSource = FALSE,
    ColumnSelectionDynamicSource = FALSE,
    RowSelectionRestrict = FALSE,
    ColumnSelectionRestrict = FALSE,
    SelectionHistory = list()
  )

################################################################################
# Settings for Feature assay plot 1
################################################################################

initial[["FeatureAssayPlot1"]] <-
  new(
    "FeatureAssayPlot",
    Assay = "logcounts",
    XAxis = "Column data",
    XAxisColumnData = "cellType.merged",
    XAxisFeatureName = "SNAP25",
    XAxisFeatureSource = "---",
    XAxisFeatureDynamicSource = FALSE,
    YAxisFeatureName = "TH",
    YAxisFeatureSource = "RowDataTable1",
    YAxisFeatureDynamicSource = TRUE,
    FacetRowByColData = "Sample",
    FacetColumnByColData = "Sample",
    ColorByColumnData = "cellType.merged",
    ColorByFeatureNameAssay = "logcounts",
    ColorBySampleNameColor = "#FF0000",
    ShapeByColumnData = "Sample",
    SizeByColumnData = "sum",
    FacetRowBy = "None",
    FacetColumnBy = "None",
    ColorBy = "Column data",
    ColorByDefaultColor = "#000000",
    ColorByFeatureName = "SNAP25",
    ColorByFeatureSource = "---",
    ColorByFeatureDynamicSource = FALSE,
    ColorBySampleName = "3_AAACCCAAGTAGAATC-1",
    ColorBySampleSource = "---",
    ColorBySampleDynamicSource = FALSE,
    ShapeBy = "None",
    SizeBy = "None",
    SelectionAlpha = 0.1,
    ZoomData = numeric(0),
    BrushData = list(),
    VisualBoxOpen = FALSE,
    VisualChoices = "Color",
    ContourAdd = FALSE,
    ContourColor = "#0000FF",
    PointSize = 1,
    PointAlpha = 1,
    Downsample = FALSE,
    DownsampleResolution = 200,
    CustomLabels = FALSE,
    CustomLabelsText = "3_AAACCCAAGTAGAATC-1",
    FontSize = 1,
    LegendPointSize = 1,
    LegendPosition = "Bottom",
    HoverInfo = TRUE,
    LabelCenters = FALSE,
    LabelCentersBy = "Sample",
    LabelCentersColor = "#000000",
    VersionInfo = list(iSEE = structure(
      list(c(2L, 6L, 0L)), class = c("package_version", "numeric_version")
    )),
    PanelId = c(FeatureAssayPlot = 1L),
    PanelHeight = 600L,
    PanelWidth = 6L,
    SelectionBoxOpen = FALSE,
    RowSelectionSource = "---",
    ColumnSelectionSource = "---",
    DataBoxOpen = FALSE,
    RowSelectionDynamicSource = FALSE,
    ColumnSelectionDynamicSource = FALSE,
    RowSelectionRestrict = FALSE,
    ColumnSelectionRestrict = FALSE,
    SelectionHistory = list()
  )
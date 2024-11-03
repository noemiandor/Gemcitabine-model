# Load DiagramR library
library(DiagrammeR)

# Define the graph
DiagrammeR::grViz("
digraph G {

# Define graph layout and style
graph [layout = dot, rankdir = TB]

# Track Classification section
subgraph cluster_Track_Classification {
    label = 'Track Classification';
    style = filled;
    color = lightgrey;

    x_t [label='Define parental cell at last timepoint (x_t)']
    a_b [label='Identify daughter cells']
    post_division_tracks [label='post-division tracks']
    inter_division_tracks [label='inter-division tracks']

    x_t -> a_b [label='Cell division event']
    a_b -> post_division_tracks [label='Locate tracks after division']
    a_b -> inter_division_tracks [label='Identify tracks spanning divisions']
}

# Classification of live cells to WGD and cell cycle states section
subgraph cluster_Classification_WGD_Cycle {
    label = '';
    style = filled;
    color = lightblue;

    fold_change [label='Calculate fold changes in cell size']
    train_model [label='Train Gaussian model for WGD probability']
    gauss_model [label='Apply Gaussian model']
    wgd_decision [label='Assign WGD']
    svm_classification [label='SVM classification for cell cycle phase (G1/S, G2/M)']

    inter_division_tracks -> fold_change 
    fold_change -> train_model
    train_model -> gauss_model
    gauss_model -> wgd_decision [label='probability > 0.5']
    gauss_model -> svm_classification [label='probability < 0.5']
}

# Connect post-division tracks to Gaussian model
post_division_tracks -> gauss_model
}")

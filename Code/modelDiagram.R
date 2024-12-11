library(DiagrammeR)

diagram <- grViz(
  "digraph cell_cycle_system {
    
    # Define node shapes and styles
    node [shape = rectangle, style = filled, fillcolor = lightblue, fontname = Helvetica]
    G1S [label = \"G1S\"]
    G2M [label = \"G2M\"]
    NP [label = \"N_P\"]
    A [label = \"A (Dead Cells)\"]
    
    
    # Rank nodes to enforce layout
    {rank=same; G1S; G2M; NP}
    {rank=same; A}
    
    # Define edges with labels for transitions
    edge [fontname = Helvetica, fontsize = 10]
    G2M -> G1S [label = \"2 * f_m * k_2\"]
    G1S -> G2M [label = \"f_s * k_1\"]
    G1S -> A [label = \"k_d * Cgem\"]
    G2M -> A [label = \"k_d * Cgem\"]
    G2M -> NP [label = \"k_NP\"]
    NP -> A [label = \"k_d2\"]
    
  }")

# Render the diagram
diagram


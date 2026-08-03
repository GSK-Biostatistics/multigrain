# create_edges() works

    Code
      create_edges(trans_matrix)
    Output
      # A tibble: 12 x 3
         from  to     value
         <chr> <chr>  <dbl>
       1 H1    H2    0.0986
       2 H1    H3    0.439 
       3 H1    H4    0.462 
       4 H2    H1    0.489 
       5 H2    H3    0.465 
       6 H2    H4    0.0457
       7 H3    H1    0.193 
       8 H3    H2    0.165 
       9 H3    H4    0.642 
      10 H4    H1    0.233 
      11 H4    H2    0.466 
      12 H4    H3    0.301 

# create_nodes() works

    Code
      create_nodes(hyp_weight, edges = edges_df)
    Output
      # A tibble: 4 x 3
        hypothesis weight optimised
        <chr>       <dbl> <lgl>    
      1 H1          0.125 TRUE     
      2 H2          0.176 TRUE     
      3 H3          0.270 TRUE     
      4 H4          0.429 TRUE     


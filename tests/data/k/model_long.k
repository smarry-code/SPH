$ model_long: widths strictly per k_cards.hpp
*KEYWORD LONG=Y
*TITLE
Small test model
*PART
Impactor
                   1                   1                   1
*PART
Fiber
                   2
*NODE
                   1                 0.0                 0.0                 0.0                                        
                   2                0.25                 0.0                 0.0                   0                   0
                   3                 0.5                 0.0                 0.0                                        
                   4                0.75                 0.0              2.5E-1                                        
                  10                 0.0                 1.0                 0.0                                        
                  11                 0.5                                     0.0                                        
                  12              1.0D+0                 1.0                -2.0                                        
                  13                 0.0                 0.0                 5.0                                        
                  20                 0.0                 2.0                 0.0                                        
                  21                 1.0                 2.0                 0.0                                        
                  22                 1.0                 3.0                 0.0                                        
                  23                 0.0                 3.0                 0.0                                        
                  24                 0.0                 2.0                 1.0                                        
                  25                 1.0                 2.0                 1.0                                        
                  26                 1.0                 3.0                 1.0                                        
                  27                 0.0                 3.0                 1.0                                        
*ELEMENT_SPH
                   1                   1                 0.5                   3
                   4                   1              -0.125                    
*ELEMENT_BEAM
                   1                   2                  10                  11                  13                                                                                                    
                   2                   2                  11                  12                                                                                                                       3
*ELEMENT_SOLID
                   1                   1                  20                  21                  22                  23                  24                  25                  26                  27
*SET_NODE_LIST
                   5                                                                                                MECH
                  10                  11                  12                  13                  20                  21                  22                  23
                  24                  25
*SET_NODE_LIST_GENERATE
                   6
                  20                  23                  25                  27
*BOUNDARY_SPC_NODE
                  10                                       1                   1                   1
*BOUNDARY_SPC_SET
                   5                                       1                                        
*MAT_ELASTIC
         1     1000.0   1.0E+09       0.3
*END

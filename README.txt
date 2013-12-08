
===================
GENERATION
===================

===================
TYPE: linear in L*
-------------------
SUBTYPE: one direction
...................
clab_blackwhite   -> clab_blackwhite_make
clab_hot          -> clab_hot_make          -> makecmap_pinchedspiral
clab_purple       -> clab_purple_make       -> makecmap_pinchedspiral
-------------------
SUBTYPE: two direction
...................
clab_bluewhitered -> clab_bluewhitered_make -> makecmap_AwpBlin

===================
TYPE: fixed L*
-------------------
clab_greenpurple  -> clab_greenpurple_make  -> makecmap_ABlin
clab_yellowblue   -> clab_yellowblue_make   -> makecmap_ABlin
clab_Larc
clab_Lfix_disting
clab_wheel


===================
DEVELOPMENT
===================

find_AB            for makecmap_ABlin         <- greenpurple, yellowblue
find_AwpB          for makecmap_AwpBlin       <- bluewhitered
find_AwpBcurve     for makecmap_AwpBcurve     <- bluewhitered
find_AwpBtwist     for makecmap_AwpBtwist     <- "hot-cold"
find_pinchedspiral for makecmap_pinchedspiral <- hot

fit_cbrewer        for fitting colormaps to colorbrewer tables

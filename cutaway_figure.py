import pymol
from pymol import cmd

cmd.load("/Users/jonathanborowsky/Documents/grabelab/projects/rutter-lab/draft_paper_segment/figure_5_section/figures/cutaway_view/seg_0035_ucp1_centered.gro", "protein")

cmd.set_view((
     0.999923825,   -0.011909503,   -0.002981398,\
     0.003197649,    0.018239507,    0.999828696,\
    -0.011851273,   -0.999763131,    0.018275371,\
    -0.000038579,   -0.000087558, -313.712249756,\
    48.857982635,   46.785575867,   43.795089722,\
   292.577178955,  334.846099854,  -20.000000000))

cmd.hide("sticks", "elem H")
cmd.hide("lines")
cmd.hide("spheres")
cmd.hide("nb_spheres")

cmd.set_color("c_protein", [57/256, 79/256, 141/256])
cmd.color("c_protein", "poly")

cmd.set_color("c_headgroup", [233/256, 182/256, 78/256])
cmd.set_color("c_tails", [255/256, 126/256, 126/256])


#POPC
cmd.color("c_headgroup", "resn POPC and \
          name N+C11+C12+C13+C14+C15+\
P+O11+O12+O13+O14+\
C1+C2+C3+\
C31+O31+O32+\
C21+O21+O22")
#choline, PO4, glycerol (carbons only), and the two FA esters

cmd.color("c_tails", "resn POPC and \
          not name N+C11+C12+C13+C14+C15+\
P+O11+O12+O13+O14+\
C1+C2+C3+\
C31+O31+O32+\
C21+O21+O22")


#TLCL2
cmd.color("c_headgroup", "resn TLCL2 and \
          name OG12+C1+C2+C3+\
P1+OP11+OP12+OP13+OP14+\
C11+C12+C13+\
CA1+OA1+O12+\
CB1+OB1+O13+\
P3+OP31+OP32+OP33+OP34+\
C31+C32+C33+\
CC1+OC1+O32+\
CD1+OD1+O33")
#central glycerol, 2 x (phosphate, glycerol carbons, and 2 FA esters)

cmd.color("c_tails", "resn TLCL2 and \
          not name OG12+C1+C2+C3+\
P1+OP11+OP12+OP13+OP14+\
C11+C12+C13+\
CA1+OA1+O12+\
CB1+OB1+O13+\
P3+OP31+OP32+OP33+OP34+\
C31+C32+C33+\
CC1+OC1+O32+\
CD1+OD1+O33")


#ARAN
cmd.color("c_headgroup", "resn ARAN and \
          name O1+O2+C1")
cmd.color("c_tails", "resn ARAN and \
          not name O1+O2+C1")


#TLCL2 and POPC
cmd.show("spheres", "elem P or name P3 or name P1")


#TIP3P
cmd.show("spheres", "resn TIP3 and elem O")
cmd.set("sphere_scale", 0.35, "resn TIP3 and elem O")

cmd.set_color("c_water", [68/256, 155/256, 223/256])
cmd.color("c_water", "resn TIP3 and elem O")

cmd.set("ray_shadows", 0)
cmd.png("figures/cutaway_figure_red.png", width=2000, height=2000, dpi=600, ray=1)
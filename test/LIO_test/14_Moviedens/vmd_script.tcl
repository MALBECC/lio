# $Id: au-iso.vmd,v 1.2 2004/05/21 15:50:29 akohlmey Exp $ 
# Display settings 
display projection Orthographic
display nearclip set 0.000000
display farclip set 10.000000
display depthcue off
# store the molecule id for later use 
set updmol [mol new {frame0000.cube} type cube waitfor all]
 mol addfile {frame0000.cube} type cube waitfor all
 mol addfile {frame0001.cube} type cube waitfor all
 mol addfile {frame0002.cube} type cube waitfor all
 mol addfile {frame0003.cube} type cube waitfor all
 mol addfile {frame0004.cube} type cube waitfor all
 mol addfile {frame0005.cube} type cube waitfor all
 mol addfile {frame0006.cube} type cube waitfor all
 mol addfile {frame0007.cube} type cube waitfor all
 mol addfile {frame0008.cube} type cube waitfor all
 mol addfile {frame0009.cube} type cube waitfor all
 mol addfile {frame0010.cube} type cube waitfor all

# Representation 0: The molecule 
mol delrep 0 top
mol representation lines
mol color Name
mol selection {all}
mol addrep top

# Representation 1: Positive isosurface 

mol color ColorID 0
mol representation Isosurface 0.000200 0.000000 0.000000 0.000000 1 1
mol selection {all}
mol addrep top
# Representation 2: Negative isosurface 

mol material Opaque
mol color ColorId 1
mol representation Isosurface -0.0002 0.0 0.0 0.0  1 1
mol selection {all}
mol addrep top

mol material Opaque
mol color ColorId 1
mol representation Isosurface 0.001 0.0 0.0 0.0  1 1
mol selection {all}
mol addrep top

mol material Opaque
mol color ColorId 1
mol representation Isosurface -0.001 0.0 0.0 0.0  1 1
mol selection {all}
mol addrep top



set updrep1 [mol repname top 1]
set updrep2 [mol repname top 2]
set updrep3 [mol repname top 3]
set updrep4 [mol repname top 4]
mol rename top {test}
# use the volumetric data set for the isosurface corresponding to the frame. 
# $updmol contains the id of the molecule and $updrep the (unique) name of 
# the isosurface representation 
proc update_iso {args} {
    global updmol
    global updrep1
    global updrep2
    global updrep3
    global updrep4

    # get representation id and return if invalid 
    set repid1 [mol repindex $updmol $updrep1]
    if {$repid1 < 0} { return }
    set repid2 [mol repindex $updmol $updrep2]
    if {$repid2 < 0} { return }
    set repid3 [mol repindex $updmol $updrep3]
    if {$repid2 < 0} { return }
    set repid4 [mol repindex $updmol $updrep4]
    if {$repid2 < 0} { return }
    # update representation but replace the data set 
    # id with the current frame number. 
    set frame [molinfo $updmol get frame]
    lassign [molinfo $updmol get "{rep $repid1}"] rep
    mol color colorid 0
    mol representation [lreplace $rep 2 2 $frame]
    mol modrep $repid1 $updmol
    lassign [molinfo $updmol get "{rep $repid3}"] rep
    mol color colorid 0
    mol representation [lreplace $rep 2 2 $frame]
    mol modrep $repid3 $updmol

    lassign [molinfo $updmol get "{rep $repid2}"] rep
    mol color colorid 1
    mol representation [lreplace $rep 2 2 $frame]
    mol modrep $repid2 $updmol
    lassign [molinfo $updmol get "{rep $repid4}"] rep
    mol color colorid 1
    mol representation [lreplace $rep 2 2 $frame]
    mol modrep $repid4 $updmol
}
trace variable vmd_frame($updmol) w update_iso
animate goto 0


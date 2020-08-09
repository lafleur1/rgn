#Alyssa La Fleur
#adapting the full atom minimization movers from the trRosetta protocol to try and clean up the PULCHAR outputs a bit, introduce real SS to the backbones

import numpy as np
import random
import json
import tempfile
from pyrosetta import *
from pyrosetta.rosetta.protocols.minimization_packing import MinMover
from pyrosetta.toolbox import cleanATOM
from rosetta.protocols.simple_moves import *
from pyrosetta.teaching import *
from rosetta.protocols.relax import *
#init('-hb_cen_soft -relax:default_repeats 5 -default_max_cycles 200 -out:level 100')
init('-out:level 100')

def remove_clash(scorefxn, mover, pose):
    # moves if there is a clash ????
    for _ in range(0, 5):
        if float(scorefxn(pose)) < 10:
            break
        mover.apply(pose)


def fulltrRosettaPyRosettaRelaxProtocol(pose, outname):
    # score function and movers
    sf = ScoreFunction()
    sf.add_weights_from_file('./trRosettaPyRosetta/data/scorefxn.wts')

    sf1 = ScoreFunction()
    sf1.add_weights_from_file('./trRosettaPyRosetta/data/scorefxn1.wts')

    sf_vdw = ScoreFunction()
    sf_vdw.add_weights_from_file('./trRosettaPyRosetta/data/scorefxn_vdw.wts')

    sf_cart = ScoreFunction()
    sf_cart.add_weights_from_file('./trRosettaPyRosetta/data/scorefxn_cart.wts')

    mmap = MoveMap()
    mmap.set_bb(True)
    mmap.set_chi(False)
    mmap.set_jump(True)

    min_mover = MinMover(mmap, sf, 'lbfgs_armijo_nonmonotone', 0.0001, True)
    min_mover.max_iter(1000)

    min_mover1 = MinMover(mmap, sf1, 'lbfgs_armijo_nonmonotone', 0.0001, True)
    min_mover1.max_iter(1000)

    min_mover_vdw = MinMover(mmap, sf_vdw, 'lbfgs_armijo_nonmonotone', 0.0001, True)
    min_mover_vdw.max_iter(500)

    min_mover_cart = MinMover(mmap, sf_cart, 'lbfgs_armijo_nonmonotone', 0.0001, True)
    min_mover_cart.max_iter(1000)
    min_mover_cart.cartesian(True)

    repeat_mover = RepeatMover(min_mover, 3)

    # mutate GLY to ALA
    for i, a in enumerate(seq):
        if a == 'G':
            mutator = rosetta.protocols.simple_moves.MutateResidue(i + 1, 'ALA')
            mutator.apply(pose)
            print('mutation: G%dA' % (i + 1))

    # remove initial clashes
    remove_clash(sf_vdw, min_mover_vdw, pose)

    print('short + medium + long')
    repeat_mover.apply(pose)
    min_mover_cart.apply(pose)
    remove_clash(sf_vdw, min_mover1, pose)

    # mutate ALA back to GLY
    for i, a in enumerate(seq):
        if a == 'G':
            mutator = rosetta.protocols.simple_moves.MutateResidue(i + 1, 'GLY')
            mutator.apply(pose)
            print('mutation: A%dG' % (i + 1))

    ########################################################
    # full-atom refinement
    ########################################################

    sf_fa = create_score_function('ref2015')
    sf_fa.set_weight(rosetta.core.scoring.atom_pair_constraint, 5)
    sf_fa.set_weight(rosetta.core.scoring.dihedral_constraint, 1)
    sf_fa.set_weight(rosetta.core.scoring.angle_constraint, 1)

    mmap = MoveMap()
    mmap.set_bb(True)
    mmap.set_chi(True)
    mmap.set_jump(True)

    relax = rosetta.protocols.relax.FastRelax()
    relax.set_scorefxn(sf_fa)
    relax.max_iter(200)
    relax.dualspace(True)
    relax.set_movemap(mmap)

    switch = SwitchResidueTypeSetMover("fa_standard")
    switch.apply(pose)

    print('relax...')
    params['PCUT'] = 0.15
    # add_rst(pose, rst, 1, len(seq), params, True) #remove restrains
    relax.apply(pose)

    ########################################################
    # save final model
    ########################################################
    pose.dump_pdb(outname)


def faRefinementFromTrRosetta(cleanPDBName, outname, max_iters=1000):
    print ("fast relax")
    ########################################################
    # full-atom refinement
    ########################################################
    pose = pose_from_pdb(cleanPDBName)

    sf_fa = create_score_function(
        'ref2015')  # optimized energy function or score function called ref2015 for calculating the energy of all atomic interactions in a globular protein made of L-amino acids
    sf_fa.set_weight(rosetta.core.scoring.atom_pair_constraint, 5)
    sf_fa.set_weight(rosetta.core.scoring.dihedral_constraint, 1)
    sf_fa.set_weight(rosetta.core.scoring.angle_constraint, 1)

    mmap = MoveMap()
    mmap.set_bb(True)
    mmap.set_chi(True)
    mmap.set_jump(True)

    relax = rosetta.protocols.relax.FastRelax()
    relax.set_scorefxn(sf_fa)
    relax.max_iter(max_iters)
    relax.dualspace(True)
    relax.set_movemap(mmap)

    switch = SwitchResidueTypeSetMover("fa_standard")
    switch.apply(pose)

    print('relax...')
    # params['PCUT'] = 0.15
    # add_rst(pose, rst, 1, len(seq), params, True)
    relax.apply(pose)


    ########################################################
    # save final model
    ########################################################
    pose.dump_pdb(outname)
    return sf_fa(pose)

def vanillaMinimizationMovers(cleanPDBName, outname, iterations = 100):
    #try using teh default movers for a set number of iterations
    print ("vanilla relax")
    pose = pose_from_pdb(cleanPDBName)
    min_mover = MinMover()
    movemapDefault = MoveMap() #standard move map
    scorefxn = get_fa_scorefxn() #standard fa score function
    min_mover.movemap(movemapDefault)
    min_mover.score_function(scorefxn)
    for i in range(0, iterations):
        min_mover.apply(pose)
    pose.dump_pdb(outname)
    return scorefxn(pose)

def classicRelax(cleanPDBName,outname):
    #try the classic relax protocol from   Bradley,  Misura, &  Baker  2005
    print ('classic relax')
    pose = pose_from_pdb(cleanPDBName)
    relax = ClassicRelax()
    scorefxn = get_fa_scorefxn()
    relax.set_scorefxn(scorefxn)
    relax.apply(pose)
    pose.dump_pdb(outname)
    return scorefxn(pose)

def runSuiteOfMethods(name):
    # try for one, clean fa pdb adn then try relaxing it
    cleanATOM(name)
    cleanPDBName = name.replace(".pdb", ".clean.pdb")
    start = pose_from_pdb(cleanPDBName)
    # fa fast relax refinement process from trRosetta
    faName = "fastRelax_"+name.replace(".pdb","")+"_"
    score1 = faRefinementFromTrRosetta(cleanPDBName, faName+"200.pdb", 200)
    score2 = faRefinementFromTrRosetta(cleanPDBName, faName+"500.pdb", 500)
    score3 = faRefinementFromTrRosetta(cleanPDBName, faName+"1000.pdb", 1000)
    # simple min mover
    vanillaName="vanilla"+name.replace(".pdb","")+"_"
    score4 = vanillaMinimizationMovers(cleanPDBName, vanillaName+"200.pdb", 200)
    score5 = vanillaMinimizationMovers(cleanPDBName, vanillaName+"500.pdb", 500)
    score6 = vanillaMinimizationMovers(cleanPDBName, vanillaName+"1000.pdb", 1000)
    # classic relax (long time)
    score7 = classicRelax(cleanPDBName, "classicRelax_" +name.replace(".pdb","") +".pdb")
    #make a scores report
    #f = open("scoreReport"+name.replace(".pdb","")+".txt", "w")
    print ("original score: ")
    scorefxn = get_fa_scorefxn()
    print (scorefxn(start))
    print ( " ")
    print ("Fast relax scores: ")
    print ("200: ")
    print (score1)
    print("500: ")
    print(score2)
    print("1000: ")
    print(score3)
    print (" ")
    print ("Vanilla min mover scores: ")
    print("200: ")
    print(score4)
    print("500: ")
    print(score5)
    print("1000: ")
    print(score6)
    print(" ")
    print ("Classic relax score: ")
    print (score7)

if __name__ == '__main__':
    print ("al q prediction: ")
    runSuiteOfMethods("negative_279_long.pdb")
    print ("trRosetta prediction: ")
    runSuiteOfMethods("negative_279_trRosetta.pdb")

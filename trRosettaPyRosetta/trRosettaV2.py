#modifying trRosetta to run on multiple predictions

import sys,os,json
import tempfile
import numpy as np
import argparse
from utils_ros import *
from pyrosetta import *
from pyrosetta.rosetta.protocols.minimization_packing import MinMover
import time
os.environ["OPENBLAS_NUM_THREADS"] = "1"
# init PyRosetta (one time, moved from main)
init('-hb_cen_soft -relax:default_repeats 5 -default_max_cycles 200 -out:level 100')

def set_args(params, fastaName, npzName, modelName, maxiter1, maxiter2, maxiter3, maxiter4, maxiter5):
    #manually filling arg parser for new variables
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-NPZ", type=str, help="input distograms and anglegrams (NN predictions)")
    parser.add_argument("-FASTA", type=str, help="input sequence")
    parser.add_argument("-OUT", type=str, help="output model (in PDB format)")

    parser.add_argument('-pd', type=float, dest='pcut', default=params['PCUT'],
                        help='min probability of distance restraints')
    parser.add_argument('-m', type=int, dest='mode', default=2, choices=[0, 1, 2],
                        help='0: sh+m+l, 1: (sh+m)+l, 2: (sh+m+l)')
    parser.add_argument('-w', type=str, dest='wdir', default=params['WDIR'], help='folder to store temp files')
    parser.add_argument('-n', type=int, dest='steps', default=1000, help='number of minimization steps')
    parser.add_argument('--orient', dest='use_orient', action='store_true', help='use orientations')
    parser.add_argument('--no-orient', dest='use_orient', action='store_false')
    parser.add_argument('--fastrelax', dest='fastrelax', action='store_true', help='perform FastRelax')
    parser.add_argument('--no-fastrelax', dest='fastrelax', action='store_false')
    parser.set_defaults(use_orient=True)
    parser.set_defaults(fastrelax=True)

    args = parser.parse_args(["-NPZ", npzName, "-FASTA", fastaName, "-OUT", modelName])

    params['PCUT'] = args.pcut
    params['USE_ORIENT'] = args.use_orient

    return args

def generateOneStructure(fastaName, npzName, modelName, ):
    #generates one pdb structure from a npz, returns score function score
    #Namespace(FASTA='T1008.fasta', NPZ='T1008.npz', OUT='model.pdb', fastrelax=True, mode=2, pcut=0.05, steps=1000, use_orient=True, wdir='/dev/shm')
    timeRunningFile = open(fastaName.replace(".fasta","") + "timeToProcess.txt", "w")
    startTime = time.time()
    ########################################################
    # process inputs
    ########################################################

    # read params
    scriptdir = os.path.dirname(os.path.realpath(__file__))
    with open(scriptdir + '/data/params.json') as jsonfile:
        params = json.load(jsonfile)

    # get command line arguments
    args = set_args(params, fastaName, npzName, modelName)
    print(args)


    # Create temp folder to store all the restraints
    tmpdir = tempfile.TemporaryDirectory(prefix=args.wdir+'/')
    params['TDIR'] = tmpdir.name
    print('temp folder:     ', tmpdir.name)

    # read and process restraints & sequence
    npz = np.load(args.NPZ)
    seq = read_fasta(args.FASTA)
    L = len(seq)
    timeRunningFile.write("LENGTH," + str(L) + "\n")
    params['seq'] = seq
    rst = gen_rst(npz,tmpdir,params)
    seq_polyala = 'A'*len(seq)
    timeRunningFile.write("SETUP," + str(time.time() - startTime) + "\n")
    startTime = time.time()

    ########################################################
    # Scoring functions and movers
    ########################################################

    sf = ScoreFunction()
    sf.add_weights_from_file(scriptdir + '/data/scorefxn.wts')

    sf1 = ScoreFunction()
    sf1.add_weights_from_file(scriptdir + '/data/scorefxn1.wts')

    sf_vdw = ScoreFunction()
    sf_vdw.add_weights_from_file(scriptdir + '/data/scorefxn_vdw.wts')

    sf_cart = ScoreFunction()
    sf_cart.add_weights_from_file(scriptdir + '/data/scorefxn_cart.wts')

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
    timeRunningFile.write("MOVERS," + str(time.time() - startTime) + "\n")
    startTime = time.time()

    ########################################################
    # initialize pose
    ########################################################
    pose = pose_from_sequence(seq, 'centroid' )

    # mutate GLY to ALA
    for i,a in enumerate(seq):
        if a == 'G':
            mutator = rosetta.protocols.simple_moves.MutateResidue(i+1,'ALA')
            mutator.apply(pose)
            print('mutation: G%dA'%(i+1))

    set_random_dihedral(pose)
    remove_clash(sf_vdw, min_mover_vdw, pose)

    timeRunningFile.write("INTIALIZE," + str(time.time() - startTime) + "\n")
    startTime = time.time()
    ########################################################
    # minimization
    ########################################################

    if args.mode == 0:

        # short
        print('short')
        add_rst(pose, rst, 1, 12, params)
        repeat_mover.apply(pose)
        min_mover_cart.apply(pose)
        remove_clash(sf_vdw, min_mover1, pose)

        # medium
        print('medium')
        add_rst(pose, rst, 12, 24, params)
        repeat_mover.apply(pose)
        min_mover_cart.apply(pose)
        remove_clash(sf_vdw, min_mover1, pose)

        # long
        print('long')
        add_rst(pose, rst, 24, len(seq), params)
        repeat_mover.apply(pose)
        min_mover_cart.apply(pose)
        remove_clash(sf_vdw, min_mover1, pose)

    elif args.mode == 1:

        # short + medium
        print('short + medium')
        add_rst(pose, rst, 3, 24, params)
        repeat_mover.apply(pose)
        min_mover_cart.apply(pose)
        remove_clash(sf_vdw, min_mover1, pose)

        # long
        print('long')
        add_rst(pose, rst, 24, len(seq), params)
        repeat_mover.apply(pose)
        min_mover_cart.apply(pose)
        remove_clash(sf_vdw, min_mover1, pose)

    elif args.mode == 2:

        # short + medium + long
        print('short + medium + long')
        add_rst(pose, rst, 1, len(seq), params)
        repeat_mover.apply(pose)
        min_mover_cart.apply(pose)
        remove_clash(sf_vdw, min_mover1, pose)


    # mutate ALA back to GLY
    for i,a in enumerate(seq):
        if a == 'G':
            mutator = rosetta.protocols.simple_moves.MutateResidue(i+1,'GLY')
            mutator.apply(pose)
            print('mutation: A%dG'%(i+1))

    timeRunningFile.write("MINIMIZATION," + str(time.time() - startTime) + "\n")
    startTime = time.time()
    ########################################################
    # full-atom refinement
    ########################################################

    if args.fastrelax == True:

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

        pose.remove_constraints()
        switch = SwitchResidueTypeSetMover("fa_standard")
        switch.apply(pose)

        print('relax...')
        params['PCUT'] = 0.15
        add_rst(pose, rst, 1, len(seq), params, True)
        relax.apply(pose)
    timeRunningFile.write("FASTRELAX," + str(time.time() - startTime) + "\n")
    startTime = time.time()
    ########################################################
    # save final model
    ########################################################

    pose.dump_pdb(args.OUT)
    #return final energy score
    scorefxn = get_fa_scorefxn()
    finalScore = scorefxn(pose)
    timeRunningFile.write("SCORE," + str(finalScore) + "\n")
    timeRunningFile.close()
    return finalScore


if __name__ == '__main__':
    fastaName = "negative_279.fasta"
    npzName = "seq.npz"
    modelName = "demoNeg279_lowQ.pdb"
    score = generateOneStructure(fastaName, npzName, modelName)
    print ("FINAL ENERGY SCORE: ", score)

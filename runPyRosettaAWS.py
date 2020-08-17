
#python trRosetta.py T1008.npz T1008.fasta model.pdb

#adapting trRosetta PyRosetta scripts to run as functions - generating SKEMPI chain predicted PDBs
import numpy as np
import random
import json
import tempfile
from pyrosetta import *
from pyrosetta.rosetta.protocols.minimization_packing import MinMover
import os
import time

def gen_rst(npz, tmpdir, params):
    #generate restraints
    dist,omega,theta,phi = npz['dist'][0],npz['omega'][0],npz['theta'][0],npz['phi'][0]
    print ("dist: ", dist.shape)
    print ("omega: ", omega.shape)
    print ("theta: ", theta.shape)
    print ("phi: ", phi.shape)
    # dictionary to store Rosetta restraints
    rst = {'dist' : [], 'omega' : [], 'theta' : [], 'phi' : [], 'rep' : []}

    ########################################################
    # assign parameters
    ########################################################
    PCUT  = 0.05 #params['PCUT']
    PCUT1 = params['PCUT1']
    EBASE = params['EBASE']
    EREP  = params['EREP']
    DREP  = params['DREP']
    PREP  = params['PREP']
    SIGD  = params['SIGD']
    SIGM  = params['SIGM']
    MEFF  = params['MEFF']
    DCUT  = params['DCUT']
    ALPHA = params['ALPHA']

    DSTEP = params['DSTEP']
    ASTEP = np.deg2rad(params['ASTEP'])

    seq = params['seq']
    print ("number res: ", len(seq))
    ########################################################
    # repultion restraints
    ########################################################
    #cbs = ['CA' if a=='G' else 'CB' for a in params['seq']]
    '''
    prob = np.sum(dist[:,:,5:], axis=-1)
    i,j = np.where(prob<PREP)
    prob = prob[i,j]
    for a,b,p in zip(i,j,prob):
        if b>a:
            name=tmpdir.name+"/%d.%d_rep.txt"%(a+1,b+1)
            rst_line = 'AtomPair %s %d %s %d SCALARWEIGHTEDFUNC %.2f SUMFUNC 2 CONSTANTFUNC 0.5 SIGMOID %.3f %.3f\n'%('CB',a+1,'CB',b+1,-0.5,SIGD,SIGM)
            rst['rep'].append([a,b,p,rst_line])
    print("rep restraints:   %d"%(len(rst['rep'])))
    '''


    ########################################################
    # dist: 0..20A
    ########################################################
    nres = dist.shape[0]
    #print ("number res: ", nres)
    bins = np.array([4.25+DSTEP*i for i in range(32)])
    #print ("bins: ", bins.shape)
    prob = np.sum(dist[:,:,5:], axis=-1)
    #print ("prob: ", prob.shape)
    bkgr = np.array((bins/DCUT)**ALPHA)
    #print ("bkgr: ", bkgr.shape)
    attr = -np.log((dist[:,:,5:]+MEFF)/(dist[:,:,-1][:,:,None]*bkgr[None,None,:]))+EBASE
    repul = np.maximum(attr[:,:,0],np.zeros((nres,nres)))[:,:,None]+np.array(EREP)[None,None,:]
    dist = np.concatenate([repul,attr], axis=-1)
    bins = np.concatenate([DREP,bins])
    i,j = np.where(prob>PCUT)
    prob = prob[i,j]
    nbins = 35
    step = 0.5
    for a,b,p in zip(i,j,prob):
        if b>a:
            name=tmpdir.name+"/%d.%d.txt"%(a+1,b+1)
            with open(name, "w") as f:
                f.write('x_axis'+'\t%.3f'*nbins%tuple(bins)+'\n')
                f.write('y_axis'+'\t%.3f'*nbins%tuple(dist[a,b])+'\n')
                f.close()
            rst_line = 'AtomPair %s %d %s %d SPLINE TAG %s 1.0 %.3f %.5f'%('CB',a+1,'CB',b+1,name,1.0,step)
            rst['dist'].append([a,b,p,rst_line])
    print("dist restraints:  %d"%(len(rst['dist'])))


    ########################################################
    # omega: -pi..pi
    ########################################################
    nbins = omega.shape[2]-1+4
    bins = np.linspace(-np.pi-1.5*ASTEP, np.pi+1.5*ASTEP, nbins)
    prob = np.sum(omega[:,:,1:], axis=-1)
    i,j = np.where(prob>PCUT)
    prob = prob[i,j]
    omega = -np.log((omega+MEFF)/(omega[:,:,-1]+MEFF)[:,:,None])
    omega = np.concatenate([omega[:,:,-2:],omega[:,:,1:],omega[:,:,1:3]],axis=-1)
    for a,b,p in zip(i,j,prob):
        if b>a:
            name=tmpdir.name+"/%d.%d_omega.txt"%(a+1,b+1)
            with open(name, "w") as f:
                f.write('x_axis'+'\t%.5f'*nbins%tuple(bins)+'\n')
                f.write('y_axis'+'\t%.5f'*nbins%tuple(omega[a,b])+'\n')
                f.close()
            rst_line = 'Dihedral CA %d CB %d CB %d CA %d SPLINE TAG %s 1.0 %.3f %.5f'%(a+1,a+1,b+1,b+1,name,1.0,ASTEP)
            rst['omega'].append([a,b,p,rst_line])
    print("omega restraints: %d"%(len(rst['omega'])))


    ########################################################
    # theta: -pi..pi
    ########################################################
    prob = np.sum(theta[:,:,1:], axis=-1)
    i,j = np.where(prob>PCUT)
    prob = prob[i,j]
    theta = -np.log((theta+MEFF)/(theta[:,:,-1]+MEFF)[:,:,None])
    theta = np.concatenate([theta[:,:,-2:],theta[:,:,1:],theta[:,:,1:3]],axis=-1)
    for a,b,p in zip(i,j,prob):
        if b!=a:
            name=tmpdir.name+"/%d.%d_theta.txt"%(a+1,b+1)
            with open(name, "w") as f:
                f.write('x_axis'+'\t%.3f'*nbins%tuple(bins)+'\n')
                f.write('y_axis'+'\t%.3f'*nbins%tuple(theta[a,b])+'\n')
                f.close()
            rst_line = 'Dihedral N %d CA %d CB %d CB %d SPLINE TAG %s 1.0 %.3f %.5f'%(a+1,a+1,a+1,b+1,name,1.0,ASTEP)
            rst['theta'].append([a,b,p,rst_line])
            #if a==0 and b==9:
            #    with open(name,'r') as f:
            #        print(f.read())
    print("theta restraints: %d"%(len(rst['theta'])))


    ########################################################
    # phi: 0..pi
    ########################################################
    nbins = phi.shape[2]-1+4
    bins = np.linspace(-1.5*ASTEP, np.pi+1.5*ASTEP, nbins)
    prob = np.sum(phi[:,:,1:], axis=-1)
    i,j = np.where(prob>PCUT)
    prob = prob[i,j]
    phi = -np.log((phi+MEFF)/(phi[:,:,-1]+MEFF)[:,:,None])
    phi = np.concatenate([np.flip(phi[:,:,1:3],axis=-1),phi[:,:,1:],np.flip(phi[:,:,-2:],axis=-1)], axis=-1)
    for a,b,p in zip(i,j,prob):
        if b!=a:
            name=tmpdir.name+"/%d.%d_phi.txt"%(a+1,b+1)
            with open(name, "w") as f:
                f.write('x_axis'+'\t%.3f'*nbins%tuple(bins)+'\n')
                f.write('y_axis'+'\t%.3f'*nbins%tuple(phi[a,b])+'\n')
                f.close()
            rst_line = 'Angle CA %d CB %d CB %d SPLINE TAG %s 1.0 %.3f %.5f'%(a+1,a+1,b+1,name,1.0,ASTEP)
            rst['phi'].append([a,b,p,rst_line])
            #if a==0 and b==9:
            #    with open(name,'r') as f:
            #        print(f.read())

    print("phi restraints:   %d"%(len(rst['phi'])))

    return rst

def set_random_dihedral(pose):
    nres = pose.total_residue()
    for i in range(1, nres):
        phi,psi=random_dihedral()
        pose.set_phi(i,phi)
        pose.set_psi(i,psi)
        pose.set_omega(i,180)

    return(pose)


#pick phi/psi randomly from:
#-140  153 180 0.135 B
# -72  145 180 0.155 B
#-122  117 180 0.073 B
# -82  -14 180 0.122 A
# -61  -41 180 0.497 A
#  57   39 180 0.018 L
def random_dihedral():
    phi=0
    psi=0
    r=random.random()
    if(r<=0.135):
        phi=-140
        psi=153
    elif(r>0.135 and r<=0.29):
        phi=-72
        psi=145
    elif(r>0.29 and r<=0.363):
        phi=-122
        psi=117
    elif(r>0.363 and r<=0.485):
        phi=-82
        psi=-14
    elif(r>0.485 and r<=0.982):
        phi=-61
        psi=-41
    else:
        phi=57
        psi=39
    return(phi, psi)


def read_fasta(file):
    fasta=""
    with open(file, "r") as f:
        for line in f:
            if(line[0] == ">"):
                continue
            else:
                line=line.rstrip()
                fasta = fasta + line;
    return fasta


def remove_clash(scorefxn, mover, pose):
    for _ in range(0, 5):
        if float(scorefxn(pose)) < 10:
            break
        mover.apply(pose)


def add_rst(pose, rst, sep1, sep2, params, nogly=False):

    pcut=params['PCUT']
    seq = params['seq']

    array=[]

    if nogly==True:
        array += [line for a,b,p,line in rst['dist'] if abs(a-b)>=sep1 and abs(a-b)<sep2 and seq[a]!='G' and seq[b]!='G' and p>=pcut]
        if params['USE_ORIENT'] == True:
            array += [line for a,b,p,line in rst['omega'] if abs(a-b)>=sep1 and abs(a-b)<sep2 and seq[a]!='G' and seq[b]!='G' and p>=pcut+0.5] #0.5
            array += [line for a,b,p,line in rst['theta'] if abs(a-b)>=sep1 and abs(a-b)<sep2 and seq[a]!='G' and seq[b]!='G' and p>=pcut+0.5] #0.5
            array += [line for a,b,p,line in rst['phi'] if abs(a-b)>=sep1 and abs(a-b)<sep2 and seq[a]!='G' and seq[b]!='G' and p>=pcut+0.6] #0.6
    else:
        array += [line for a,b,p,line in rst['dist'] if abs(a-b)>=sep1 and abs(a-b)<sep2 and p>=pcut]
        if params['USE_ORIENT'] == True:
            array += [line for a,b,p,line in rst['omega'] if abs(a-b)>=sep1 and abs(a-b)<sep2 and p>=pcut+0.5]
            array += [line for a,b,p,line in rst['theta'] if abs(a-b)>=sep1 and abs(a-b)<sep2 and p>=pcut+0.5]
            array += [line for a,b,p,line in rst['phi'] if abs(a-b)>=sep1 and abs(a-b)<sep2 and p>=pcut+0.6] #0.6


    if len(array) < 1:
        return

    random.shuffle(array)

    # save to file
    tmpname = params['TDIR']+'/minimize.cst'
    with open(tmpname,'w') as f:
        for line in array:
            f.write(line+'\n')
        f.close()

    # add to pose
    constraints = rosetta.protocols.constraint_movers.ConstraintSetMover()
    constraints.constraint_file(tmpname)
    constraints.add_constraints(True)
    constraints.apply(pose)

    os.remove(tmpname)

def generatePDB(fasta, npz, out):

    ########################################################
    # process inputs
    ########################################################

    # read params
    #scriptdir = os.path.dirname(os.path.realpath(__file__))
    with open('./trRosettaPyRosetta/data/params.json') as jsonfile:
        params = json.load(jsonfile)

    args = {'pcut':params['PCUT'], 'mode':2,  'wdir':params['WDIR'], 'steps':1000, 'use_orient':True, 'fast_relax':True}
    #set FASTA, NPZ, OUT
    args['FASTA'] = fasta
    args['NPZ'] = npz
    args['OUT'] = out
    params['PCUT'] = args['pcut']
    params['USE_ORIENT'] = args['use_orient']
    # init PyRosetta
    init('-hb_cen_soft -relax:default_repeats 5 -default_max_cycles 200 -out:level 100')

    # Create temp folder to store all the restraints
    tmpdir = tempfile.TemporaryDirectory(prefix=args['wdir']+'/')
    params['TDIR'] = tmpdir.name
    print('temp folder:     ', tmpdir.name)

    # read and process restraints & sequence
    npz = np.load(args['NPZ'])
    print (type(npz))
    seq = read_fasta(args['FASTA'])
    L = len(seq)
    print (L)
    params['seq'] = seq
    rst = gen_rst(npz,tmpdir,params)
    seq_polyala = 'A'*len(seq)


    ########################################################
    # Scoring functions and movers
    ########################################################
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

    ########################################################
    # minimization
    ########################################################

    if args['mode'] == 0:

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

    elif args['mode'] == 1:

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

    elif args['mode'] == 2: #default behavior

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


    ########################################################
    # full-atom refinement
    ########################################################

    if args['fast_relax'] == True:

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

    ########################################################
    # save final model
    ########################################################
    pose.dump_pdb(args['OUT'])


def runOnDirectories():
    #generate test PDBs for trRosetta predictions currently done

    #for each MSA, attempt trRosetta structural prediction
    #for cutoffLevel in cutoffs:
    predDir = "./SKEMPI_CHAIN_FILES/trRosettaOutputs/"
    fastaDir = "./SKEMPI_CHAIN_FILES/fastaFiles/"
    pdbDir = "./SKEMPI_CHAIN_FILES/pdbStructures/"

    allFiles = os.listdir(predDir)
    print (len(os.listdir(predDir)))
    failedFiles = []
    for i in range(0,len(allFiles)):
        print ("ON ", i , " OUT OF ", len(allFiles))
        fileName = allFiles[i]
        try:
            nameComps = fileName.replace("_keras_xaa.npz","").split("_") #need recombine all but last list entry
            actualName = "_".join(nameComps[0:len(nameComps)-1])
            pdbName = "_".join(nameComps[0:len(nameComps)])
            print (fileName, actualName)
            fastaFile = fastaDir + actualName + ".fasta"
            npzFile = predDir +fileName
            pdbName = pdbDir + pdbName + ".pdb"
            print ("PDB: ", pdbName)
            print ("NPZ: ", npzFile)
            print ("FASTA: ", fastaFile)
            if not os.path.isfile(pdbName) and os.path.isfile(npzFile):
                #if str(cutoffLevel) in fileName:
                generatePDB(fastaFile, npzFile, pdbName)
            else:
                print ("pdb exists")
        except:
            print ("attempt failed, problem with trRosetta output?")
            failedFiles.append(fileName)

    allPDBs = os.listdir(pdbDir)
    print (len(allPDBs))
    print ("Failed files: ")
    print (failedFiles)

if __name__ == '__main__':
    fasta = "./TR015032_results/negative_279.fasta"
    predNPZ = "./TR015032_results/seq.npz"
    startTime = time.time()
    generatePDB(fasta, predNPZ, "trRosettaPredV1.pdb")
    print ("OVERALL TIME: ")
    print (time.time() - startTime)

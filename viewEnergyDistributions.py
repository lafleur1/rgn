#view the distriution of fa scores, # h-bonds, and percent SS for true and false HLH structures

import pandas as pd
import sys,os,json
import tempfile
import numpy as np
import argparse
from pyrosetta import *
from pyrosetta.rosetta.protocols.minimization_packing import MinMover
import time
from pyrosetta.toolbox import cleanATOM
import seaborn as sns
import matplotlib.pyplot as plt

pyrosetta.init()

def getSSString(pose):
    #returns string of SS positions
    DSSP = pyrosetta.rosetta.protocols.moves.DsspMover()
    DSSP.apply(pose)  # populates the pose's Pose.secstruct, get sec struct for sequence
    return pose.secstruct()

def getHBondInfo(pose):
    #parses get_hbonds output,]
    seqLen = len(pose.sequence())
    hbondset = pose.get_hbonds()
    allhbonds = hbondset.hbonds()
    numberBonds =  len(allhbonds)
    donors = []
    acceptors = []
    donorPairs = [] #(resi,atm)
    acceptorPairs=[] #resi, atm
    for hbond in allhbonds:
        donors.append(hbond.don_res())
        #print ("DONOR")
        #print (hbond.don_res())
        #print (hbond.don_hatm())
        #print (pose.residue(hbond.don_res()))
        #print(pose.residue(hbond.don_res()).atom(hbond.don_hatm()))
        #print(pose.residue(hbond.don_res()).atom_name(hbond.don_hatm()))
        donorPairs.append(hbond.don_res(), pose.residue(hbond.don_res()).atom_name(hbond.don_hatm()))
        print ("ACCEPTOR")
        print (hbond.acc_res())
        print (hbond.acc_atm())
        print(pose.residue(hbond.acc_res()))
        print(pose.residue(hbond.acc_res()).atom(hbond.acc_atm()))
        print(pose.residue(hbond.acc_res()).atom_name(hbond.acc_atm()))
        acceptorPairs.append(hbond.acc_res(), pose.residue(hbond.acc_res()).atom_name(hbond.acc_atm()))
        #print (hbond.don_res(), hbond.acc_res())
    #create donor and acceptor string reps
    donorSet = list(set(donors)).sort()
    acceptorSet = list(set(acceptors)).sort()
    donorString = "" #D if donor, A if acceptor, X if neither
    acceptorStrign = ""
    for i in range(1, seqLen + 1): #residues are 1 indexed in pyrosetta
        if i in donors:
            donorString = donorString + "D"
        else:
            donorString = donorString + "X"
        if i in acceptors:
            acceptorStrign = acceptorStrign + "A"
        else:
            acceptorStrign = acceptorStrign + "X"
    return donorString, acceptorStrign, numberBonds

def getHBondResiduesAndAtomsV1(pose):
    #parses get_hbonds material  into a usage form
    #does locations as (residue, atom)
    seqLen = len(pose.sequence())
    hbondset = pose.get_hbonds()
    allhbonds = hbondset.hbonds()
    donorPairs = [] #(resi,atm)
    acceptorPairs=[] #resi, atm
    for hbond in allhbonds:
        cleanedUpDonorAtom = pose.residue(hbond.don_res()).atom_name(hbond.don_hatm()).replace(" ","")
        donorstr = str(hbond.don_res())+":"+cleanedUpDonorAtom
        donorPairs.append(donorstr)
        cleanedupAccAtom = pose.residue(hbond.acc_res()).atom_name(hbond.acc_atm()).replace(" ","")
        accStr = str(hbond.acc_res()) + ":" + cleanedupAccAtom
        acceptorPairs.append(accStr)
    return donorPairs, acceptorPairs

def getHBondResiduesAndAtomsV2(pose):
    #parses get_hbonds material  into a usage form
    #does locations as (chain, residue, atom)
    seqLen = len(pose.sequence())
    hbondset = pose.get_hbonds()
    allhbonds = hbondset.hbonds()
    donorPairs = [] #(resi,atm)
    acceptorPairs=[] #resi, atm
    for hbond in allhbonds:
        donorlocation = pose.pdb_info().pose2pdb(hbond.don_res()).split(" ") #1st is residue, second is chain ID in pdb
        cleanedUpDonorAtom = pose.residue(hbond.don_res()).atom_name(hbond.don_hatm()).replace(" ","")
        donorstr = donorlocation[1] + ":" + donorlocation[0] +":"+cleanedUpDonorAtom
        donorPairs.append(donorstr)
        acclocation = pose.pdb_info().pose2pdb(hbond.acc_res()).split(" ")
        cleanedupAccAtom = pose.residue(hbond.acc_res()).atom_name(hbond.acc_atm()).replace(" ","")
        accStr = acclocation[1] + ":" + acclocation[0] + ":" + cleanedupAccAtom
        acceptorPairs.append(accStr)
    return donorPairs, acceptorPairs


def processOne(name, scorefxn):
    cleanATOM(name)
    cleanPDBName = name[:-4] +  ".clean.pdb"
    start = pose_from_pdb(cleanPDBName)
    #ss
    ss = getSSString(start)
    fractionH = ss.count("H") / len(ss)
    fractionL = ss.count("L") / len(ss)
    #hbonds
    donors, acceptors, numberBonds = getHBondInfo(start)
    # fa score
    score = scorefxn(start)
    print ("ON ", name)
    print ("Secondary structure: ")
    print (ss)
    print ("fraction H: ", fractionH, " fraction L ", fractionL, " number H bonds: ", numberBonds, " score: ", score)
    return ss, fractionH, fractionL, numberBonds, score


def processAllFolders():
    truePosDir = "./testSetChains/"
    myPosDir = "./positiveHLHTestSet/"
    myNegDir = "./negativeHLHTestSet/"
    scorefxn = get_fa_scorefxn()
    asDict = {'structure': [], 'dataset': [], 'score': [], 'fractionH': [], 'fractionL': [], 'numberHBonds': []}
    for file in os.listdir(myNegDir):
        ss, fractionH, fractionL, numberBonds, score = processOne(myNegDir + file, scorefxn)
        asDict['score'].append(score)
        asDict['fractionH'].append(fractionH)
        asDict['fractionL'].append(fractionL)
        asDict['numberHBonds'].append(numberBonds)
        # asDict['s'].append(ss)
        asDict['dataset'].append("negativeTrRosettaPredicted")
        asDict['structure'].append(file)
    for file in os.listdir(myPosDir):
        ss, fractionH, fractionL, numberBonds, score = processOne(myPosDir + file, scorefxn)
        asDict['score'].append(score)
        asDict['fractionH'].append(fractionH)
        asDict['fractionL'].append(fractionL)
        asDict['numberHBonds'].append(numberBonds)
        # asDict['s'].append(ss)
        asDict['dataset'].append("positiveTrRosettaPredicted")
        asDict['structure'].append(file)
    for file in os.listdir(truePosDir):
        ss, fractionH, fractionL, numberBonds, score = processOne(truePosDir + file, scorefxn)
        asDict['score'].append(score)
        asDict['fractionH'].append(fractionH)
        asDict['fractionL'].append(fractionL)
        asDict['numberHBonds'].append(numberBonds)
        # asDict['s'].append(ss)
        asDict['dataset'].append("positiveOriginalSet")
        asDict['structure'].append(file)
    # graph interesting stuff
    asDF = pd.DataFrame(asDict)  # change to dataframe for sns plotting
    # save csv
    asDF.to_csv("relevantData.csv")
    sns.violinplot(x="dataset", y="score", data=asDF)
    plt.show()
    sns.violinplot(x="dataset", y="numberHBonds", data=asDF)
    plt.show()
    ax = sns.scatterplot(x="fractionH", y="fractionL", hue="dataset", data=asDF)
    plt.show()

if __name__ == '__main__':
    #true positive structures:
    #look at multichain hbonds
    name = "./180312_massive_set_train_pdbs/redesigned_closed_7_7_7_9middlesbobby_1_9_S_237903.pdb_middle1.pdb-7_7_7_9middlesbobby_1_9_S_254850_0001.pdb_middle1.pdb-bobby_1_9_S_246831_padded_0001.pdb"
    cleanATOM(name)
    cleanPDBName = name[:-4] + ".clean.pdb"
    start = pose_from_pdb(cleanPDBName)
    x,y = getHBondResiduesAndAtomsV2(start)
    print (x)
    print (y)
    print (start.get_hbonds())

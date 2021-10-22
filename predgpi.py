#!/usr/bin/env python
import argparse
import os
import json
PREDGPI_HOME=os.environ.get('PREDGPI_HOME')
# moduli miei
import sys,string,numpy
dirbin=os.path.join(PREDGPI_HOME, 'GPIDAT')
from predgpilib.hmm import HMM_IO, algo_HMM
from predgpilib.svm import SVMLike
from predgpilib import utils
try:
    import psyco
    psyco.full()
except ImportError:
    pass


def fitFPR(x):
    ''' fitFPR(x)
          provide an extimate of the False Positive Rate using
          ROC_ok Pierleoni data
          x is the svm prediction
           => it has been chosen so that for a svm value > -0.61 the TPR < 0.005 on the dataset
    '''
    return 1/(1.0+numpy.exp(-6.8*(-x-1.39)))

def runHMM(seq,hmm):
    ''' '''
    (bestpath,bestval)=algo_HMM.viterbi(hmm,seq)
    spos=1
    i=0
    while i < len(bestpath) and hmm.states[bestpath[i]].name != 'Mt' :
        if bestpath[i] in hmm.emits:
            spos+=1
        i+=1
    return len(seq)-spos,algo_HMM.seq_log_Prob(hmm, seq, Scale='Yes')/len(seq)


def mksvmInput(seq,lphmm):
    ''' mksvmInput(seq,lphmm)
        frequenze intera sequenza   (0,20)
        frequenze degli ultimi 40   (20,40)
        frequenze degli ultimi 20   (40,60)
        frequenze dei primi 20      (60,80)
        kd somma/4.5*len last   20  (80,81)
        kd somma/4.5*len primi  21  (81,82)
        logprob hmm degli ultimi 40 (82,83)
    '''
    # DATA that has to be use

    allSeqN=1.0/len(seq)
    n20=1.0/20.0
    n40=1.0/40.0
    nComp=80 # 80 components but only 75 will be used


    kd_scale={'A':1.800, 'R':-4.500, 'N': -3.500, 'D': -3.500,  'C':2.500,
'Q': -3.500,  'E': -3.500, 'G': -0.400,  'H': -3.200,  'I': 4.500,
'L': 3.800,  'K': -3.900,  'M': 1.900,  'F': 2.800,  'P': -1.600,
'S': -0.800,  'T': -0.700,  'W': -0.900,  'Y': -1.300,  'V': 4.200}

    resPos={'A': 10, 'C': 19, 'E': 4, 'D': 3, 'G': 18, 'F': 15, 'I': 13, 'H': 2, 'K': 0, 'M': 16, 'L': 12, 'N': 5, 'Q': 6, 'P': 14, 'S': 7, 'R': 1, 'T': 8, 'W': 17, 'V': 11, 'Y': 9}

    componentToUse=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 34, 35, 38, 39, 40, 41, 42, 43, 44, 45, 46, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 77, 78, 79, 80, 81, 82]
    # end DATA
    tmpVec=[0.0]*nComp
    # whole sequence
    for aa in seq:
        if aa in resPos.keys():
            tmpVec[resPos[aa]]+=allSeqN # whole sequence
    # last 40
    for aa in seq[-40:]:
        if aa in resPos.keys():
            tmpVec[20+resPos[aa]]+=n40
    # last 20
    kdl20=0.0
    for aa in seq[-20:]:
        if aa in resPos.keys():
            tmpVec[40+resPos[aa]]+=n20
            kdl20+=kd_scale[aa]
    kdl20/=20.0*4.5
    kdf20=0.0
    for aa in seq[0:20]:
        if aa in resPos.keys():
            tmpVec[60+resPos[aa]]+=n20
    for aa in seq[0:21]:
        if aa in resPos.keys():
            kdf20+=kd_scale[aa]
    kdf20/=21.0*4.5
    tmpVec+=[kdl20,kdf20,lphmm]
    nv=numpy.zeros(len(componentToUse),float)
#    print 1,
    for i in range(len(componentToUse)):
        nv[i]=tmpVec[componentToUse[i]]
#        print str(i+1)+':'+str(nv[i]),
#    print
    return nv

def predGpipe(sequence,svm,hmmmod):
    ''' predGpipe(sequence,svm,hmmmod) '''
    #print("##",sequence[-40:])
    cut,lprob=runHMM(sequence[-40:],hmmmod)
    svminput=mksvmInput(sequence,-lprob)
    svmout=svm.predict(svminput)
    return lprob,cut,svmout,fitFPR(svmout)


def readFasta(fname):
    ''' readFasta(fname) return a SEQS object
    '''
    import sys
    try:
        lines=open(fname,'r').readlines()
    except:
        sys.stderr.write('cannot open file '+fname+'\n')
        sys.exit()
    seqs={}
    name=''
    seq=''
    for line in lines:
        if line[0]=='>': # assuming fasta file
           if name != '':
               seqs[name]=seq
               seq=''
           name=''.join(line[1:].split())
        else:
           seq+=''.join(line.split())
    if name != '':
        seqs[name]=seq
    return seqs

def printVal(sequence,svm,hmmmod):
    lprob,cut,svmout,fitFPR=predGpipe(sequence,svm,hmmmod)
    print("Extimated False Positive Rate =",fitFPR)
    print("HMM log(probability) =",lprob)
    print("GPI-Anchor length=",cut)
    if fitFPR <= 0.0015: # threshold to accept the sequence as GPI
        print("Highly Probable GPI-Anchored protein")
    elif fitFPR <= 0.005:
        print("Probable GPI-Anchored protein")
    elif fitFPR <= 0.01:
        print("Weakly GPI-Anchored protein")
    else:
        print("Improbable GPI-Anchored protein")



###########

def main():
    DESC = "PredGPI: Prediction of GPI-anchors in proteins"
    parser = argparse.ArgumentParser(description = DESC, prog = "predgpi.py")
    parser.add_argument("-f", "--fasta", help = "The input sequences in FASTA format", dest = "fasta", required = True)
    parser.add_argument("-o", "--output", help = "The output file name", dest="outf", required = "True")
    parser.add_argument("-c", "--conservative", help = "Conservative omega (opt)", dest="conservative", action = "store_true")
    parser.add_argument("-m", "--outfmt", help = "The output format: json or gff3 (default)", choices=['json', 'gff3'], required = False, default = "gff3")
    ns = parser.parse_args()

    sequences = readFasta(ns.fasta)
    if ns.conservative:
        hmmmod=HMM_IO.get_hmm(os.path.join(dirbin, 'PHMM.TOT.ss.mod_CSDGN'))
    else:
        hmmmod=HMM_IO.get_hmm(os.path.join(dirbin, 'PHMM.TOT.ss.mod'))
    svm=SVMLike.getSVMLight(os.path.join(dirbin, 'MOD'))
    ofs = open(ns.outf, 'w')
    protein_jsons = []
    for name in sequences:
        seq=sequences[name]
        seq_t=seq.replace("U","C")
        seq_t=seq_t.replace("Z","A")
        seq_t=seq_t.replace("B","A")
        seq_t=seq_t.replace("X","A")
        if len(seq) > 40:
            lprob,cut,svmout,fitFPR=predGpipe(seq_t,svm,hmmmod)
        else:
            fitFPR = 1
        if fitFPR <= 0.01:
            gpi = True
            cleavage = len(seq) - cut
            if fitFPR <= 0.0015:
                prob = 1.0
            elif fitFPR <= 0.005:
                prob = 0.70
            else:
                prob = 0.55
        else:
            gpi = False
            cleavage = "-"
            prob = 1.0
        if ns.outfmt == "gff3":
            utils.write_gff_output(name, seq, ofs, gpi, cleavage, prob)
        else:
            acc_json = utils.get_json_output(name, seq, gpi, cleavage, prob)
            protein_jsons.append(acc_json)
    if ns.outfmt == "json":
        json.dump(protein_jsons, ofs, indent=5)
    ofs.close()


if __name__=='__main__':
    main()

'''
    Copyright (C) 2007 Piero Fariselli

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Copyright (C) 2007  Piero Fariselli
    contacts: Piero Fariselli
              e-mail: piero@biocomp.unibo.it, piero.fariselli@unibo.it
              Dept. of Computer Science and Engineering
              University of Bologna
              Italy
    This program comes with ABSOLUTELY NO WARRANTY; for details type psCoils.py.
    This is free software, and you are welcome to redistribute it
    under the condition of preserving this text
'''

import sys,string,numpy
import _pickle as cPickle
try:
	import psyco
	psyco.full()
except ImportError:
	pass


class SVMLike:
    ''' This class implement the prediction phase of a SVM
        based on the model of svm light.
        Presently, there are several limitation including:
           - fixed dimension of the input vectors ()
           - !!! only rbf kernel implemente up to now !!!!
    '''

    def __init__(self,sv,ai,b,kernel,params):
        '''__init__(self,sv,ai,b,kerneltype,params)
               sv     -> support vectors
               ai     -> lagrange multipliers
               b      -> shifting threshold
               kernel -> 0: linear (default)
                         1: polynomial (s a*b+c)^d
                         2: radial basis function exp(-gamma ||a-b||^2)
                         3: sigmoid tanh(s a*b + c)
               params -> params['dim']: highest feature index
                         params['-d']:  degree polynomial
                         params['-g']:  gamma in rbf
                         params['-s']:  parameter s in sigmoid/poly kernel
                         params['-r']:  parameter c in sigmoid/poly kernel
        '''
        self.sv=sv
        self.numsv=len(sv)
        self.ai=ai
        self.b=b
        if kernel == 0:
           self.kernel=self._klin
        elif kernel == 1:
           self.kernel=self._kpoly
        elif kernel == 2:
           self.kernel=self._krbf
        elif kernel == 3:
           self.kernel=self._ksig
        else:
           self.kernel=_klin
        self.dim=params['dim']
        self.d=params['-d']
        self.g=params['-g']
        self.s=params['-s']
        self.t=params['-r']

    def _klin(self,v1,v2):
        ''' fake kernel '''
        return self._krbf(v1,v2)

    def _kpoly(self,v1,v2):
        ''' fake kernel '''
        return self._krbf(v1,v2)

    def _krbf(self,v1,v2):
        ''' _krbf(self,v1,v2)'''
        x=v1-v2
        return numpy.exp(-self.g*(numpy.dot(x,x)))

    def predict(self,x):
        kval=0.0
        for i in range(self.numsv):
            kval+=self.ai[i]*self.kernel(self.sv[i],x)
        return kval-self.b

    def svmSave(self,fname):
        ''' svmSave(self,fname) '''
        f=open(fname,'w')
        cPickle.dump(self,f)
        f.close()

#---------------------------------#

def getsvmPickle(fname):
    ''' getsvmPickle(fname) '''
    f=open(fname)
    svm=cPickle.load(f)
    return svm

def unpacksvmVec(x,vecDim):
    ''' unpacksvmVec(x,vecDim) '''
    nv=numpy.zeros(vecDim,float)
    v=x.split('#')[0] # if contains comment
    v=v.split()
    first=float(v[0])
    for e in v[1:]:
        (idx,val)=e.split(':')
        idx,val=int(idx),float(val)
        nv[idx-1]=val
    return first,nv

def getSVMLight(fname):
    ''' getSVMMod(fname) '''
    lines=open(fname).readlines()
    # set kernel type
    i=0
    while lines[i].find('kernel type')<0:
        i+=1
    #kernel=string.atoi(lines[i].split()[0])
    kernel=int(lines[i].split()[0])
    # parameters
    params={}
    while lines[i].find('threshold b')<0: # loop until threshold
        if lines[i].find('kernel parameter') >=0:
            v=lines[i].split()
            if v[-1] != '-u':
               params[v[-1]]=float(v[0])
        elif lines[i].find('highest feature index')>=0:
            params['dim']=int(lines[i].split()[0])
        i+=1

    b=float(lines[i].split()[0])
    # print "threshold found ",i
    i+=1
    ai=[]
    sv=[]
    while i < len(lines):
        a,v=unpacksvmVec(lines[i],params['dim'])
        ai.append(a)
        sv.append(v)
        i+=1
    return SVMLike(sv,ai,b,kernel,params)


if __name__=='__main__':
    pass

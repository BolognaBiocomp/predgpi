'''
This file contains the class trainable object (TR_OBJ).
It consists of the type of the object we want to train with 
our HMM
'''

import math

##########################################
class TR_OBJ:
    ''' define a trainable object
        a TR_OBJ object cositst of:
                                      
	self.len 	# length of the sequence
        self.seq	# the list of symbols || vectors
        self.labels	# the list of the labels
        self.f		# the forward matrix
        self.b		# the backward matrix
        self.eMat       # the precalculate emission matrix
        self.scale 	# the scale factors
        self.prob	# the probability of the sequence
        self.lprob	# the log(probability) of the sequence
        self.name	# the object name
                                                                         
    '''
    def __init__(self,seq,labels=[],f=[],b=[],eMat=[],scale=[],p=None,name=None):
        ''' __init__(self,seq,labels=[],f=[],b=[],scale=[],p=None) 
        self.seq	# the list of symbols || vectors
        self.labels	# the list of the labels
        self.f		# the forward matrix
        self.b		# the backward matrix
        self.scale 	# the scale factors
        self.prob	# the probability of the sequence
        self.lprob	# the log(probability) of the sequence
        self.name	# the object name
        '''
        self.seq=seq		# the list of symbols || vectors
	self.len=len(seq) 	# length of the sequence
        self.labels=labels	# the list of the labels
        self.f=f		# the forward matrix
        self.b=b		# the backward matrix
        self.eMat=eMat          # the precalculate emission matrix
        self.scale=scale 	# the scale factors
        self.prob=p		# the probability of the sequence
        self.lprob=None
        if p:
            self.lprob=math.log(p)	# the log(probability) of the sequence
        self.name=name		# the object name
    


'''
This file contains the HMM definition
'''
try:
   import psyco
   psyco.full()
except:
   pass

import numpy as NUM

ARRAYFLOAT=NUM.float64
ARRAYINT=NUM.int


def _getRand(randint,values,names,norm=10000):
    ''' _getRand(randint,values,names,norm=10000)
    '''
    v=[]
    s=0
    for e in values:
        s=int(norm*e+1)+s
        v.append(s)
    irand=randint(0,v[-1]-1)
    i=0
    while v[i] < irand:
        i=i+1
    return i,names[i]

##########################################
class HMM:
    ''' define a HMM
        a HMM object cositst of:

	self.num_states 		# number of states
        self.states	 		# the states
        self.emission_alphabet
        self.state_names		# the state names
        self.fix_tr 	# the transitions which are not be trained
        self.fix_em	# the emissions which are not trained
        self.topo_order # the topologicrom al ordered state list (in numbers from 0)
        self.emits	# the list of emitting states
        self.nulls	# the list of the null states
        # list of the inlinks and outlinks for each states
        # also restricted to nulls and emitting
        self.out_s	# the outlinks of each state
        self.out_s_e	# the outlinks restricted to emitting states for each state
        self.out_s_n	# the outlinks restricted to null states for each state
        self.in_s	# the inlinks of each state
        self.in_s_e	# the inlinks restricted to emitting states for each state
        self.in_s_n	# the inlinks restricted to null states for each state
        # end states, all, nulls and emitting
        self.end_s	# the list of the states that can end
        self.end_s_e	# the above list restricted to the emittings
        self.end_s_n	# the above list restricted to nulls
        # now we eliminate 0 from null states Begin must be a null state!!
        self.nulls	# the list of the null states without Begin
        self.label_list # the list of the labels of all states
        self.effective_tr # the list of the effective transitions (No node_tr)
        self.effective_em # the list of the effective emissions (No node_em)
        self.effective_tr_pos # dictionary [statename] number of effective transition
        self.effective_em_pos # dictionary [statename] number of effective emissions (No node_em)
        self.hmm_tr # precalculated transitions
        self.em_updatable = intersection({not tied}, {not fixed})
        self.tr_updatable = intersection({not tied}, {not fixed})
        #
        self.labelMusk # dictionary of labels and musk values
        self.mA # fully connected matrix version of transitions
        self.mE # fully connected matrix version of emissions

    '''
    #__metaclass__ = psyco.compacttype
    def __init__(self,states,emission_alphabet,fix_tr,fix_em):
        ''' __init__(self,nstates,states,emission_alphabet) '''
        self.num_states=len(states) # number of states
        self.states=states # the states
        self.emission_alphabet=emission_alphabet
        # create list of states
        self.state_names=[] # the state names
        for i in range(self.num_states):
           self.state_names.append(states[i].name)
        # create the list of the state which have learnable
        # transitions and emissions
        self.fix_tr=[None]*self.num_states
        self.fix_em=[None]*self.num_states
        for i in range(self.num_states):
            self.fix_tr[i]=fix_tr[self.state_names[i]]
            self.fix_em[i]=fix_em[self.state_names[i]]
        from . import Sort_HMM
        # topological order of the states with the corresponding null and emitting states
        (all_links,self.topo_order,self.emits,self.nulls)=Sort_HMM.tolpological_sort(self.state_names,self.states)
        # list of the inlinks and outlinks for each states
        # also restricted to nulls and emitting
        (self.out_s,self.in_s,self.out_s_e,self.out_s_n,self.in_s_e,self.in_s_n)=Sort_HMM.make_links(self.emits,self.nulls,all_links)
        # end states, all, nulls and emitting
        (self.end_s,self.end_s_e,self.end_s_n)=Sort_HMM.make_ends(self.states,self.emits,self.nulls)
        # now we eliminate 0 from null states Begin must be a null state!!
        i=self.nulls.index(0)
        self.nulls=self.nulls[0:i]+self.nulls[i+1:]
        self.label_list=[]
        for s in self.states:
            if (s.label != 'None') and (s.label not in self.label_list):
                self.label_list.append(s.label)
        self.effective_tr=[] # the list of the effective transitions (No node_tr)
        self.effective_em=[] # the list of the effective emissions (No node_em)
        self.effective_tr_pos={} # dictionary [statename] number of effective transition
        self.effective_em_pos={} # dictionary [statename] number of effective emissions (No node_em)
        self.em_updatable = [] # these are used when we update the parameters set_param
        self.tr_updatable = [] # these are used when we update the parameters set_param
        for i in range(self.num_states):
            name=self.states[i].get_tr_name()
            if name not in self.effective_tr:
                self.effective_tr.append(name)
                if not self.fix_tr[i]:
                    self.tr_updatable.append(i)
            name=self.states[i].get_em_name()
            if i in self.emits and name not in self.effective_em:
                self.effective_em.append(name)
                if not self.fix_em[i]:
                    self.em_updatable.append(i)
        for s in self.topo_order:
            self.effective_tr_pos[s]=self.effective_tr.index(self.states[s].get_tr_name())
            if s in self.emits:
                self.effective_em_pos[s]=self.effective_em.index(self.states[s].get_em_name())
        # precalculate transitions
        self.hmm_tr={} # precalculated transitions
        self.hmm_ln_tr={} # precalculated transitions
        for s in self.topo_order:
            for t in self.topo_order:
                self.hmm_tr[(s,t)]=self.states[s].a(self.states[t])
                self.hmm_ln_tr[(s,t)]= self.states[s].ln_a(self.states[t])
        self.mE=NUM.array([[0.0]*len(self.emission_alphabet)]*self.num_states,ARRAYFLOAT)
        self.set_mE()
        self.mA=NUM.array([[0.0]*self.num_states]*self.num_states,ARRAYFLOAT)
        self.set_mA()
        self.set_labelMusk()

    def a(self,i,j):
        ''' a(self,i,j) transition probability between i an j '''
        return self.hmm_tr[(i,j)]
#        return self.states[i].a(self.states[j])

    def ln_a(self,i,j):
        ''' ln_a(self,i,j) log of the transition probability between i an j '''
        return self.hmm_ln_tr[(i,j)]
#        return self.states[i].ln_a(self.states[j])

    def set_a(self,i,j,val):
        ''' set_a(self,i,j,val) sets the transition probability between i an j '''
        self.states[i].set_a(self.states[j],val)
        self.hmm_tr[(i,j)]=val
        self.hmm_ln_tr[(i,j)]=self.states[i].ln_a(self.states[j])
        return

    def e(self,i,x):
        ''' e(self,i,x) emission probability of x in the state i '''
        return self.states[i].e(x)

    def ln_e(self,i,x):
        ''' ln_e(self,i,x) log of the emission probability of x in the state i '''
        return self.states[i].ln_e(x)

    def eV(self,i,x):
        ''' eV(self,i,x) emission probability of vector x in the state i '''
        return self.states[i].eV(x)

    def ln_eV(self,i,x):
        ''' ln_eV(self,i,x) log of the emission probability of vector x in the state i '''
        return self.states[i].ln_eV(x)

    def set_e(self,i,x,val):
        ''' set_e(self,i,x,val) sets with val the emission probability of x in the state i '''
        return self.states[i].set_e(x,val)

    def save(self,file):
        ''' save the hmm into cPickle file '''
        import HMM_IO
        HMM_IO.save_cPickle(self,file)

    def set_mA(self):
        ''' set_mA(self) set sel.mA the fully connected matrix version of the transitions'''
        for s in self.topo_order:
            for t in self.topo_order:
                self.mA[s][t]=self.a(s,t)

    def set_mE(self):
        ''' set_mA(self) set sel.mA the fully connected matrix version of the emissions'''
        for s in self.topo_order:
            i=0
            for c in self.emission_alphabet:
                self.mE[s][i]=self.e(s,c)
                i+=1
    def set_labelMusk(self):
        '''  set_labelMusk(self) set a dictionary of label and allow state values'''
        self.labelMusk={}
        for lab in self.label_list:
            self.labelMusk[lab]=NUM.array([0.0]*self.num_states,ARRAYFLOAT)
            for s in self.topo_order:
                if self.states[s].label == lab :
                    self.labelMusk[lab][s]=1.0

    def load_file(self,file):
        ''' load the hmm from file '''
        import HMM_IO
        obj=HMM_IO.load_cPickle(file)
        self.__dict__=obj.__dict__

    def write_for_humans(self,filename):
        ''' write_for_humans(self,filename)
           write on filename the self in the same format
           it is read by parse_text
        '''
        import HMM_IO
        HMM_IO.write_for_humans(self,filename)


    def randomWalk(self,maxLen):
        '''
           randomWalk(self,maxLen)
           rerutn a list containing tuples of label emission
                  of maximal length = maxLen
                  In reality if not all the states can be ending states
                  the list length can be longer or shorter!
        '''
        from random import randint
        stateNum=0
        state=self.states[0]
        name=state.name
        retlist=[]
        i=0
        while i<maxLen and  self.out_s[stateNum] :
            nextPos,name=_getRand(randint,state.get_transitions(),state.out_links)
            stateNum=self.state_names.index(name)
            state=self.states[stateNum]
            em_letters=state.em_letters
            if em_letters != [] : # not null
                nextPos,emiss=_getRand(randint,state.get_emissions(),em_letters)
                label=state.label
                retlist.append((name,label,emiss))
                i=i+1
        while not state.end_state :
            nextPos,name=_getRand(randint,state.get_transitions(),state.out_links)
            stateNum=self.state_names.index(name)
            state=self.states[stateNum]
            em_letters=state.em_letters
            if em_letters != [] : # not null
                nextPos,emiss=_getRand(randint,state.get_emissions(),em_letters)
                label=state.label
                retlist.append((name,label,emiss))
        return retlist

###########

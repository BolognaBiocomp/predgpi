'''
This file contains the definitions of:

node_tr => transitions class
node_em => emission class
State => state class which contains
            one pointer to a node_tr object
            one pointer to a node_em object
'''

from .Def import DEF
import numpy as NUM
try:
   import psyco
   psyco.full()
except:
   pass

######################################
######## class node_trans ############
######################################
class node_tr:
    ''' this class implement the node transitions '''
    #__metaclass__ = psyco.compacttype
    def __init__(self,name,tr):
        ''' __init__(self,name,tr)
            name = identifier, tr = transition probabilities
        '''
        self.name=name  # "esempio '_tr1','_tr2' ...
        self.len=len(tr)
        self._tr=tr # Vettore di transizioni
        # computes the log of the transitions
        self._ln_tr=[DEF.big_negative]*self.len
        for i in range(self.len):
            if(self._tr[i] > DEF.tolerance):
                self._ln_tr[i]=NUM.log(self._tr[i])

    def tr(self,i):
        '''
           tr(self,i) -> transition probability between self->i-th
        '''
        #assert(i<self.len)
        return self._tr[i]

    def ln_tr(self,i):
        '''
           ln(tr(self,i)) -> transition probability between self->i-th
        '''
        #assert(i<self.len)
        return self._ln_tr[i]

    def set_tr(self,i,value):
        '''
          set_tr(self,i,value) sets the i-th transition of self to value
        '''
        #assert(i<self.len)
        self._tr[i]=value
        if(self._tr[i] > DEF.tolerance):
            self._ln_tr[i]=NUM.log(self._tr[i])
        else:
            self._ln_tr[i]=DEF.big_negative

######################################

######################################
######## class node_em    ############
######################################
class node_em:
    ''' this class implement the node emissions '''
    #__metaclass__ = psyco.compacttype
    def __init__(self,name,em):
        ''' __init__(self,name,em=None)
            name = identifier, em = emission probabilities
        '''
        self.name=name  # "esempio '_tr1','_tr2' ...
        self.len=len(em) # dimensione del vettore
        self._em=em # Vettore emissioni
        # computes the log of the emissions
        self._ln_em=[DEF.big_negative]*self.len
        for i in range(self.len):
            if(self._em[i] > DEF.tolerance):
                self._ln_em[i]=NUM.log(self._em[i])

    def em(self,i):
        '''
           em(self,i) -> transition probability between self->i-th
        '''
        #assert(i<self.len)
        return self._em[i]

    def ln_em(self,i):
        '''
           em(self,i) -> transition probability between self->i-th
        '''
        #assert(i<self.len)
        return self._ln_em[i]

    def set_em(self,i,value):
        '''
           set_em(self,i,value) set the i-th emission of self to value
        '''
        #assert(i<self.len)
        self._em[i]=value
        if(self._em[i] > DEF.tolerance):
            self._ln_em[i]=NUM.log(self._em[i])
        else:
            self._ln_em[i]=DEF.big_negative

######################################

######################################
######## class State      ############
######################################
class State:
    ''' This class implement the state of a HMM '''
    #__metaclass__ = psyco.compacttype
    def __init__(self,name,n_tr,n_em,out_s,in_s,em_let,tied_t,tied_e,end_s,label=None):
        '''
            __init__(self,name,n_tr,n_em,out_s,in_s,em_let,tied_t,tied_e,end_s,label=None)
            name = state name
            n_tr = a node_tr object
            n_em = a node_em object
            out_s = the state outlinks [list of the state names]
            in_s = the state inlinks [list of the state names]
            em_let = emission letter [list in the order given by n_em]
            tied_t = if is tied to a given transition state name of None)
            tied_e = if is tied to a given emission (state name of None)
            end_s = end state flag
            label = classification attribute (None default)
            _idxem={} dictionary name:index
            _idxtr={} dictionary name:index
        '''
        self.name=name
        self._node_tr=n_tr
        self._node_em=n_em
        self.out_links=out_s
        self.in_links=in_s
        self.em_letters=em_let
        self.tied_t=tied_t
        self.tied_e=tied_e
        self.end_state=end_s
        self.label=label
        self._idxem={}
        self._idxtr={}
        for name in self.out_links:
            self._idxtr[name]=self.out_links.index(name)
        for symbol in self.em_letters:
            self._idxem[symbol]=self.em_letters.index(symbol)
        # some tests
        #assert(self._node_tr.len == len(self.out_links))
        #assert(self._node_em.len == len(self.em_letters))

    def get_tr_name(self):
        ''' get_tr_name() -> returns the name of the transitions '''
        return self._node_tr.name

    def get_em_name(self):
        ''' get_em_name() -> returns the name of the emissions '''
        return self._node_em.name

    def get_transitions(self):
        ''' get_transitions() -> returns the value of the transitions '''
        return self._node_tr._tr

    def get_emissions(self):
        ''' get_emissions() -> returns the value of the emissions '''
        return self._node_em._em

    def a(self,state):
        ''' a_{i,j} in durbin et al., 1998
           self.a(state)  -> transition probability between self->state
        '''
        if state.name in self._idxtr:
            return self._node_tr._tr[self._idxtr[state.name]]
        else:
            return(0.0)

    def set_a(self,state,value):
        ''' set the value of a_{i,j} in durbin et al., 1998
          self.a(state,value)  -> self->state = value
        '''
        self._node_tr.set_tr(self._idxtr[state.name],value)

    def e(self,symbol):
        ''' e_{k}(x) in durbin et al., 1998
            self.e(symbol)  -> emission probability in state self of  'symbol'
        '''
        if symbol in self._idxem:
            return self._node_em.em(self._idxem[symbol])
        else:
            return(0.0)

    def eV(self,vec):
        ''' e_{k}(V) in Gigi 2002
          self.e(V)  -> emission probability in state self of Vector
        '''
        #assert len(vec) == self._node_em.len
        return NUM.dot(self._node_em._em,vec)

    def set_e(self,symbol,value):
        ''' set the value of e_{k}(x) in durbin et al., 1998
          self.e(symbol,value)  -> set self.e(symbol)=value
        '''
        self._node_em.set_em(self._idxem[symbol],value)

    def ln_a(self,state):
        ''' ln(a_{i,j}) in durbin et al., 1998
            self.ln_a(state)  -> log(transition probability between self->state)
        '''
        if state.name in self._idxtr:
            return(self._node_tr.ln_tr(self._idxtr[state.name]))
        else:
            return(DEF.big_negative)

    def ln_e(self,symbol):
        ''' ln(e_{k}(x)) in durbin et al., 1998
            self.ln_e(symbol)  -> log(emission probability in state self of  'symbol')
        '''
        ce=self.e(symbol)
        if ce > 0:
            return(NUM.log(ce))
        else:
            return(DEF.big_negative)

    def ln_eV(self,vec,order=None):
        ''' ln(e_{k}(V)) gigi 2002
            self.ln_e(V)  -> log(emission probability in state self of  'vector')
        '''
        ce=self.eV(vec)
        if ce > 0:
           return(NUM.log(ce))
        else:
           return(DEF.big_negative)

########################

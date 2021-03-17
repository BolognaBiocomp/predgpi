#!/usr/bin/env python
'''
This file contains the input/output for the HMM.
The functions that read the human made HMM definitions and
their parsing.
In order to bould an HMM from scratch we need this package
'''


def get_hmm(filename):
    '''
    get_hmm(filename)
    => return the HMM read from filename
    '''
    return Build_HMM(filename).get_HMM()


def save_cPickle(hmm,file):
    ''' save the hmm into cPickle file '''
    import _pickle as cPickle
    fs=open(file,'w')
    cPickle.dump(hmm,fs)
    fs.close()

def load_cPickle(file):
    ''' load the object from file '''
    import _pickle as cPickle
    import copy
    fs=open(file,'r')
    hmm=cPickle.load(fs)
    fs.close()
    return hmm


###########################################
def parse_text(text):
   ''' parse_text(text) legge
       il contenuto dell'hmm passato attraverso text
       che a sua volta e` derivato da un file di tipo mod
   '''
   ret={} # return dictionary
   import sys
   import re
   import string
   curr_name=None # nome corrente
   for line in text:
      if (re.search('^#',line) == None):
         list=line.split()
         if (len(list)>0):
            if(list[0] == 'TRANSITION_ALPHABET' or  \
                          list[0] == 'EMISSION_ALPHABET'):
               ret[list[0]]=list[1:]
            elif(list[0] == 'NAME'):
               ret[list[1]]={}
               curr_name=list[1]
               ret[curr_name].update({'FIX_TR':None})
               ret[curr_name].update({'FIX_EM':None})
               if(curr_name not in ret['TRANSITION_ALPHABET']):
                   sys.exit(curr_name + " not in TRANSITION_ALPHABET ="+str(ret['TRANSITION_ALPHABET']))
            elif(list[0] == 'FIX_TR'):
               ret[curr_name].update({'FIX_TR':'YES'})
            elif(list[0] == 'FIX_EM'):
               ret[curr_name].update({'FIX_EM':'YES'})
            elif(list[0] == 'ENDSTATE'):
               #ret[curr_name].update({list[0]:string.atoi(list[1])})
               ret[curr_name].update({list[0]:int(list[1])})
            elif(list[0] == 'LABEL'):
               ret[curr_name].update({list[0]:list[1]})
            elif(list[0] == 'EMISSION' or list[0] == 'TRANS'):
               if(list[1] == 'None'):
                  tmplist=[]
               else:
                  tmplist=list[1:]
                  for i in range(len(tmplist)):
                     try: tmplist[i]=float(tmplist[i])
                     except: pass
               ret[curr_name].update({list[0]:tmplist})
            elif(list[0] == 'LINK' or list[0]=='EM_LIST'):
               if(list[1] == 'None'):
                  tmplist=[]
               else:
                  tmplist=list[1:]
               ret[curr_name].update({list[0]:tmplist})
   return(ret)
###########################################


def write_for_humans(hmm,filename):
    ''' write_for_humans(hmm,filename)
        write on filename the hmm in the same format
        it is read by parse_text
    '''
    import sys
    try:
       f=open(filename,'w')
    except:
       print("Can't open ", filename)
       sys.exit()
    strPrint="";
    separator="#############################\n"
    strPrint+="# alphabets\n"
    strPrint+="TRANSITION_ALPHABET "
    for i in hmm.topo_order:
       strPrint+=str(hmm.state_names[i])+" "
    strPrint+="\n"
    strPrint+="EMISSION_ALPHABET "+" ".join(hmm.emission_alphabet)+'\n'
# for each state in topological order
    strPrint+=separator
    for i in hmm.topo_order:
       strPrint+="NAME "+hmm.state_names[i]+'\n'
#trransitions
       strPrint+="LINK "
       if(hmm.out_s[i]):
           for j in hmm.out_s[i]:
               strPrint+=hmm.state_names[j]+' '
           strPrint+='\n'
       else:
           strPrint+=" None\n"
       strPrint+="TRANS "
       if(hmm.states[i].tied_t):
           strPrint+="tied "+hmm.states[i].tied_t+'\n'
       elif(hmm.out_s[i]):
           for j in hmm.out_s[i]:
              strPrint+=str(hmm.states[i].a(hmm.states[j]))+' '
           strPrint+='\n'
       else:
           strPrint+=" None\n"
       if(hmm.fix_tr[i]):
           strPrint+="FIX_TR\n"
# end state
       strPrint+="ENDSTATE "+str(hmm.states[i].end_state)+'\n'
#emissions
       if (hmm.states[i].em_letters):
           strPrint+="EM_LIST "+" ".join(hmm.states[i].em_letters)+'\n'
       else:
           strPrint+="EM_LIST None\n"
       strPrint+="EMISSION "
       if(hmm.states[i].tied_e):
           strPrint+="tied "+hmm.states[i].tied_e+'\n'
       elif(hmm.states[i].em_letters):
           v=hmm.states[i].get_emissions()
           for j in v:
               strPrint+=str(j)+' '
           strPrint+='\n'
       else:
           strPrint+=" None\n"
       if(hmm.fix_em[i]):
           strPrint+="FIX_EM\n"
#labels
       if(hmm.states[i].label):
           strPrint+="LABEL "+hmm.states[i].label+'\n'
       else:
           strPrint+="LABEL None\n"
       strPrint+=separator
    f.write(strPrint)
    f.close()

###########################################

class Build_HMM:
   ''' This class build an hmm starting from a file '''
   def __init__(self,file):
      ''' __init__(self,file)
          read from filename '''
      try:
         lines = open(file).readlines()
      except:
         print("Can't open file ",file)
      import string
      from . import State
      from . import HMM
      info=parse_text(lines)
#
      tr_al=info['TRANSITION_ALPHABET']
      if(info['EMISSION_ALPHABET'][0] == 'range'):
         em_al=[]
         for i in range(string.atoi(info['EMISSION_ALPHABET'][1])):
            em_al.append(str(i))
      else:
         em_al=info['EMISSION_ALPHABET']
      tied_t={} # tied transitions None if not tied
      tied_e={} # tied emissions None if not tied
      links={} # temporary list for each state
      states=[]
      label={}
      endstate={}
      in_links={}
      fix_tr={}
      fix_em={}
      #empty transition and emissions
      empty_tr=State.node_tr("_EMPTY_TR_",[])
      empty_em=State.node_em("_EMPTY_EM_",[])
      for name in tr_al: # initialize disctionaries
          in_links[name]=[]
          links[name]=[None,None]
      # create in_link information
      for name in tr_al:
          for in_name in info[name]['LINK']:
             if(name not in in_links[in_name]):
                in_links[in_name].append(name)
      serial=0 # used as incremental internal number for transitions
               # and emissions. It will be used toi set node_tr node_em
      for name in tr_al: # set node_tr
          if(info[name]['TRANS']!=[] and info[name]['TRANS'][0] != 'tied'):
             if(info[name]['TRANS'][0] == 'uniform'): # set uniform
                d=1.0/len(info[name]['LINK'])
                info[name]['TRANS']=[d]*len(info[name]['LINK'])
             obj=State.node_tr("_TR_"+str(serial),info[name]['TRANS'])
             serial=serial + 1
             links[name][0]=obj
             tied_t[name]=None
          if(info[name]['TRANS']==[]):
             links[name][0]=empty_tr
             tied_t[name]=None
      for name in tr_al:
          if(info[name]['TRANS']!=[] and info[name]['TRANS'][0] == 'tied'):
             tmpname=info[name]['TRANS'][1]
             links[name][0]=links[tmpname][0]
             tied_t[name]=tmpname
      serial=0
      for name in tr_al: # set node_em
          if(info[name]['EMISSION']!=[] and info[name]['EM_LIST'][0] == 'all'):
             info[name]['EM_LIST']=em_al
          if(info[name]['EMISSION']!=[] and info[name]['EMISSION'][0] != 'tied'):
             if(info[name]['EMISSION'][0] == 'uniform'): # set uniform
                d=1.0/len(info[name]['EM_LIST'])
                info[name]['EMISSION']=[d]*len(info[name]['EM_LIST'])
             obj=State.node_em("_EM_"+str(serial),info[name]['EMISSION'])
             serial=serial + 1
             links[name][1]=obj
             tied_e[name]=None
          if(info[name]['EMISSION']==[]):
             links[name][1]=empty_em
             tied_e[name]=None
      for name in tr_al:
          if(info[name]['EMISSION']!=[] and info[name]['EMISSION'][0] == 'tied'):
             tmpname=info[name]['EMISSION'][1]
             links[name][1]=links[tmpname][1]
             tied_e[name]=tmpname
      for name in tr_al: # set labels
         # fixed transitions
         if(info[name]['FIX_TR']):
             fix_tr[name]='YES'
         else:
             fix_tr[name]=None
         # fixed emissions
         if(info[name]['FIX_EM']):
             fix_em[name]='YES'
         else:
             fix_em[name]=None
         # LABELS
         if(info[name]['LABEL'] == ['None']) :
            label[name]=None
         else:
            label[name]=info[name]['LABEL']
         # end States
         endstate[name]=info[name]['ENDSTATE'] # set endstates
         # set State[i] and appnd it to the state list
         states.append(State.State(name,links[name][0],links[name][1],info[name]['LINK'],in_links[name],
                          info[name]['EM_LIST'],tied_t[name],tied_e[name],endstate[name],label[name]))
      # set self.hmm => the HMM
      self.hmm=HMM.HMM(states,em_al,fix_tr,fix_em)

   def get_HMM(self):
      ''' get_HMM() returns the HMM '''
      return self.hmm


####################################

if __name__=='__main__':
   print("Test with file mod.mod")
   print("read file and create an hmm")
   hmm=Build_HMM('mod.mod').get_HMM()
   print("inlinks (<-)  - outlinks (->) ")
   for i in range(hmm.num_states):
      print(hmm.states[i].name," -> ",hmm.states[i].out_links, " <- ",hmm.states[i].in_links)
   print("print the tied states")
   for i in range(hmm.num_states):
      n=hmm.states[i]
      if(n.tied_e):
         print("State ",n.name," has tied emission to ",n.tied_e)
      if(n.tied_t):
         print("State ",n.name," has tied transitions to ",n.tied_t)
   print("###############################################################")
   print("# state    tr_name  em_name  # state     tr_name  em_name   tying test       tying test")
   for i in range(hmm.num_states):
      ni=hmm.states[i]
      for j in range(i,hmm.num_states):
         nj=hmm.states[j]
         print("# <",ni.name,"> ",ni._node_tr.name,ni._node_em.name,)
         print("# <",nj.name,"> ",nj._node_tr.name,nj._node_em.name,)
         print(" (TR1==TR2)? ",(hmm.states[i]._node_tr == hmm.states[j]._node_tr),)
         print(" (EM1==EM2)? ",(hmm.states[i]._node_em == hmm.states[j]._node_em))
   print("===================================\nTopo_order=============================")
   for i in hmm.topo_order:
       print(i,hmm.state_names[i])

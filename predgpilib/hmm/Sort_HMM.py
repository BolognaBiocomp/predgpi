'''
This file contains some utilities used
during the construction of an HMM.
A HMM must have sorted its states according to
the fact that the emitting states came first, then
the null (or silent) states should be topologically ordered
after the emitting states.
The States in the HMM, must have the outlinks and inlinks
ordered in the same way, and this package contains the 
routines that take care of that. 
'''

##########################################
def make_ends(states,emits,nulls):
    '''
        o make_ends(states,emits,nulls)
          returns a tuple of lists o indices to states that can be end_states
          return (end_s,end_s_n,end_s_e)

        >>> states = []
    '''
    end_s=[]
    end_s_n=[]
    end_s_e=[]
    for i in emits:
       if(states[i].end_state):
           end_s_e.append(i)
           end_s.append(i)
    for i in nulls:    
       if(states[i].end_state):
           end_s_n.append(i)
           end_s.append(i)
    return(end_s,end_s_e,end_s_n)
##########################################

##########################################
def make_links_nobegin(emits,nulls,to_sort):
    '''
        o _make_links(nulls,emits,to_sort)
          computes the inlinks, outlinks for each states
          according to the order given by nulls (null states) and 
          emits (emit states)
        o returns a tuple of lists o indices to states 
          return (out_s, in_s, out_s_e, out_s_n, in_s_e, in_s_n)
          out_links, in_links for each states and restricted to
          null (out_s_n, in_s_n) and emitting (out_s_e, in_s_e)
          states 
        o the state '0' is the begin and it is not counted in the inlinks and outlinks
          since it is special and it is explicitly used in the algorithm
    '''
    num_states=len(to_sort)
    out_s=[]
    in_s=[]
    out_s_n=[]
    out_s_e=[]
    in_s_n=[]
    in_s_e=[]
    inlinks=[]
    # init links
    for i in range(num_states):
       out_s.append([])
       in_s.append([])
       out_s_n.append([])
       out_s_e.append([])
       in_s_n.append([])
       in_s_e.append([])
       inlinks.append([])
    # outlinks
    for i in range(num_states):
       if(to_sort[i]):
          for k in emits:
             if k in to_sort[i]:
                out_s_e[i].append(k)
                # create inlinks diffenet from begin
                if i not in inlinks[k] and i !=0:
                    inlinks[k].append(i)
          for k in nulls:
             if k in to_sort[i] and k != 0:
                out_s_n[i].append(k)
                # create inlinks diffenet from begin
                if i not in inlinks[k]  and i !=0:
                    inlinks[k].append(i)
       out_s[i]=out_s_e[i]+out_s_n[i] 
    # inlinks
    for i in range(num_states):
       if(inlinks[i]):
          for k in emits:
             if k in inlinks[i]:
                in_s_e[i].append(k)
          for k in nulls:
             if k in inlinks[i]:
                in_s_n[i].append(k)
       in_s[i]=in_s_e[i]+in_s_n[i] 
    return (out_s, in_s, out_s_e, out_s_n, in_s_e, in_s_n)
##########################################

##########################################
def make_sorted_links(emits,nulls,to_sort):
    '''
        o _make_sorted_links(nulls,emits,to_sort)
          computes the inlinks, outlinks for each states
          according to the order given by nulls (null states) and 
          emits (emit states) in sorted way
        o returns a tuple of lists o indices to states 
          return (out_s, in_s, out_s_e, out_s_n, in_s_e, in_s_n)
          out_links, in_links for each states and restricted to
          null (out_s_n, in_s_n) and emitting (out_s_e, in_s_e)
          states 
    '''
    num_states=len(to_sort)
    out_s=[]
    in_s=[]
    out_s_n=[]
    out_s_e=[]
    in_s_n=[]
    in_s_e=[]
    inlinks=[]
    # init links
    for i in range(num_states):
       out_s.append([])
       in_s.append([])
       out_s_n.append([])
       out_s_e.append([])
       in_s_n.append([])
       in_s_e.append([])
       inlinks.append([])
    # outlinks
    for i in range(num_states):
       if(to_sort[i]):
          for k in emits:
             if k in to_sort[i]:
                out_s_e[i].append(k)
                # create inlinks 
                if i not in inlinks[k] :
                    inlinks[k].append(i)
          for k in nulls:
             if k in to_sort[i] :
                out_s_n[i].append(k)
                # create inlinks 
                if i not in inlinks[k]  :
                    inlinks[k].append(i)
       out_s[i]=out_s_e[i]+out_s_n[i] 
    # inlinks
    for i in range(num_states):
       if(inlinks[i]):
          for k in emits:
             if k in inlinks[i]:
                in_s_e[i].append(k)
          for k in nulls:
             if k in inlinks[i]:
                in_s_n[i].append(k)
       in_s[i]=in_s_e[i]+in_s_n[i] 
    return (out_s, in_s, out_s_e, out_s_n, in_s_e, in_s_n)
##########################################
##########################################
def make_links(emits,nulls,to_sort):
    '''
        o _make_links(nulls,emits,to_sort)
          computes the inlinks, outlinks for each states
          according to the order given by nulls (null states) and 
          emits (emit states)
        o returns a tuple of lists o indices to states 
          return (out_s, in_s, out_s_e, out_s_n, in_s_e, in_s_n)
          out_links, in_links for each states and restricted to
          null (out_s_n, in_s_n) and emitting (out_s_e, in_s_e)
          states 
    '''
    num_states=len(to_sort)
    out_s=[]
    in_s=[]
    out_s_n=[]
    out_s_e=[]
    in_s_n=[]
    in_s_e=[]
    inlinks=[]
    # init links
    for i in range(num_states):
       out_s.append([])
       in_s.append([])
       out_s_n.append([])
       out_s_e.append([])
       in_s_n.append([])
       in_s_e.append([])
       inlinks.append([])
    # outlinks
    for i in range(num_states):
       if(to_sort[i]):
          for k in emits:
             if k in to_sort[i]:
                out_s_e[i].append(k)
                # create inlinks 
                if i not in inlinks[k] :
                    inlinks[k].append(i)
          for k in nulls:
             if k in to_sort[i] :
                out_s_n[i].append(k)
                # create inlinks 
                if i not in inlinks[k]  :
                    inlinks[k].append(i)
       out_s[i]=out_s_e[i]+out_s_n[i] 
    # inlinks
    for i in range(num_states):
       if(inlinks[i]):
          for k in emits:
             if k in inlinks[i]:
                in_s_e[i].append(k)
          for k in nulls:
             if k in inlinks[i]:
                in_s_n[i].append(k)
       in_s[i]=in_s_e[i]+in_s_n[i] 
    return (out_s, in_s, out_s_e, out_s_n, in_s_e, in_s_n)
##########################################
def tolpological_sort(names,states):
    ''' 
        tolpological_sort(names,states)
        topological sort using Depth First Search 
        names list of the state names
        states list of the state objects
        We are assuming that the first state is begin
        -> 
          return (all_links,sorted,sorted_emits,sorted_nulls)
           all_links is a list of list
           sorted is the sorted list of states
           sorted_emits and sorted_nulls are the list restricted
              to emitting and silent states 
    '''
    # create the list of list to sort 
    num_states=len(names)
    sorted=[0] # index(begin state) == 0
    all_links=[]
    to_sort={}
    emits=[]
    nulls=[]
    all_links=[None]*num_states  
    all_links[0]=[]  # add links to begin =0
    for s in states[0].out_links:
        all_links[0].append(names.index(s))
    for i in range(1,num_states):
       if not states[i].em_letters: # is null
          nulls.append(i)
          to_sort[i]=[]
          for s in states[i].out_links: # all null outlinks 
              if not states[names.index(s)].em_letters:
                  to_sort[i].append(names.index(s))
       else: # emitting state
          emits.append(i)
       all_links[i]=[]  # add all links
       for s in states[i].out_links:
          all_links[i].append(names.index(s))
    visited={} # flag array
    for k in to_sort.keys():
        visited[k]=None
    sorted_emits=emits # we don't care about the emitting states
    sorted_nulls=[]
    for v in visited.keys():
        if( not visited[v] ):
            __topSort(v,sorted_nulls,visited,to_sort)
    sorted.extend(sorted_emits+sorted_nulls)
    sorted_nulls.insert(0,0)
    return (all_links,sorted,sorted_emits,sorted_nulls)

def __topSort(v, sorted,visited,to_sort):
    ''' iternal topological sort '''
    visited[v]='OK'
    for w in to_sort[v]:
        if(not visited[w]):
            __topSort(w, sorted,visited,to_sort)
    sorted.insert(0,v)

##########################################


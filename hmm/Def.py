'''
Some definitions used all through the HMM
'''
class DEF:
   ''' contains definitions used
       everywhere in the hmm codes
   '''
   big_negative=-10000.0 # maximum negative number instead of log(0) 
   small_positive=1e-30
   tolerance=small_positive # if > tolerance we can take the log

"""
knoxParams.py 

netParams is a dict containing a set of network parameters using a standardized structure

simConfig is a dict containing a set of simulation configurations using a standardized structure

refs:
Destexhe, A., Contreras, D., & Steriade, M. (1998). Mechanisms underlying the 
synchronizing action of corticothalamic feedback through inhibition of thalamic 
relay cells. Journal of neurophysiology, 79(2), 999-1016.

Destexhe, A. (1998). Spike-and-wave oscillations based on the properties of 
GABAB receptors. Journal of Neuroscience, 18(21), 9099-9111.


Contributors: xxx@xxxx.com
"""

from netpyne import specs
from netpyne import sim

netParams = specs.NetParams()   # object of class NetParams to store the network parameters
simConfig = specs.SimConfig()   # object of class SimConfig to store the simulation configuration

import random as rnd
import numpy as np

def smallWorldConn(NPre, NPost, p, K):
    ''' k is smallwordness parameters
    K is ratio of connections from each pre cell to post cells
    if p=0 regular network
    if p between 0 and 1 small-world network and
    if p=1 random network 
    '''
    connMat=[]
    for i in range(NPre):
        #for j in np.arange(-1*int(np.ceil(NPost*K/2)),int(np.ceil(NPost*K/2))+1): # ring-like neighborhood
            #connMat.append([i,(NPost + i + j) % NPost]) # ring like
        for j in np.arange(i-0*1*int(np.ceil(NPost*K/2)),i+0*int(np.ceil(NPost*K/2))+1): # line-like neighborhood
            jbound = j
            #if jbound < 0: jbound = abs(j) - 1
            #if jbound > NPost-1: jbound = 2 * NPost - jbound - 1 
            if (jbound >= 0 and jbound <= NPost-1):
                connMat.append([i,jbound])
            
    if p:
        connects = [x for x in range(len(connMat))]
        rnd_ind = rnd.sample(connects, int(len(connMat)*p))
        for i in rnd_ind:
            connMat[i][1]=rnd.randint(0,NPost-1)
    return connMat

###############################################################################
#
# MPI HH TUTORIAL PARAMS
#
###############################################################################

###############################################################################
# NETWORK PARAMETERS
###############################################################################
N=100; N_PY=N; N_IN=N; N_TC=N; N_RE=N;
#N=100; N_PY=N; N_IN=N/4; N_TC=N/2; N_RE=N/2;

netParams.narrowdiam = 5
netParams.widediam = 10

netParams.xspacing = 20 # um
netParams.yspacing = 100 # um

netParams.axondelay = 0

netParams.defaultThreshold = 0

###############################################################################
# Population parameters
###############################################################################
### Cortical Cells
netParams.popParams['PY'] = {'cellType': 'PY', 'numCells': N_PY, 'cellModel': 'HH_PY', 'ynormRange': [0.1, 0.3]} #, 'yRange': [1*netParams.yspacing,1*netParams.yspacing], 'gridSpacing': netParams.xspacing}
netParams.popParams['IN'] = {'cellType': 'IN', 'numCells': N_IN, 'cellModel': 'HH_IN', 'ynormRange': [0.35, 0.45]} #, 'yRange': [2*netParams.yspacing,2*netParams.yspacing], 'gridSpacing': netParams.xspacing} 

### Thalamic cells    
netParams.popParams['TC'] = {'cellType': 'TC', 'numCells': N_TC, 'cellModel': 'HH_TC', 'ynormRange': [0.65, 0.75]} #, 'yRange': [2+3*netParams.yspacing,2+3*netParams.yspacing], 'gridSpacing': netParams.xspacing}
netParams.popParams['RE'] = {'cellType': 'RE', 'numCells': N_RE, 'cellModel': 'HH_RE', 'ynormRange': [0.8, 0.9]} #, 'yRange': [2+4*netParams.yspacing,2+4*netParams.yspacing], 'gridSpacing': netParams.xspacing}


###############################################################################
# Cell parameters list
###############################################################################

### PY (single compartment)
cellRule = netParams.importCellParams(label='PYrule', conds={'cellType': 'PY', 'cellModel': 'HH_PY'},	fileName='sPY.tem', cellName='sPY')
"""
cellRule['secs']['soma']['mechs']['hh2']={'gnabar': 0.05, 'gkbar': 0.005, 'vtraub': -55.0}
cellRule['secs']['soma']['mechs']['pas']={'g': 1.0e-4, 'e': -70}
cellRule['secs']['soma']['mechs']['im']={'gkbar': 7e-5}
"""
cellRule['secs']['soma']['vinit']=-70
netParams.cellParams['PYrule'] = cellRule

### IN (single compartment)
cellRule = netParams.importCellParams(label='INrule', conds={'cellType': 'IN', 'cellModel': 'HH_IN'},	fileName='sIN.tem', cellName='sIN')
"""
cellRule['secs']['soma']['mechs']['hh2']={'gnabar': 0.05, 'gkbar': 0.01, 'vtraub': -55.0}
cellRule['secs']['soma']['mechs']['pas']={'g': 1.5e-4, 'e': -70}
"""
cellRule['secs']['soma']['vinit']=-70
netParams.cellParams['INrule'] = cellRule

### TC (Destexhe et al., 1996; Bazhenov et al.,2002)
cellRule = netParams.importCellParams(label='TCrule', conds={'cellType': 'TC', 'cellModel': 'HH_TC'}, fileName='TC.tem', cellName='sTC')
"""
cellRule['secs']['soma']['mechs']['hh2']={'gnabar': 0.09, 'gkbar': 0.01, 'vtraub': -25.0}
cellRule['secs']['soma']['mechs']['pas']={'g': 1e-5, 'e': -70-20}
cellRule['secs']['soma']['mechs']['it']={'shift': 2.0, 'gcabar': 2.3e-3, 'q10': 3.0}
cellRule['secs']['soma']['mechs']['iar']={'shift': 0.0, 'ghbar': 1.5e-5}
cellRule['secs']['soma']['mechs']['cad']={'taur': 5.0, 'depth': 1.0, 'kt': 0.0, 'cainf': 2.4e-4, 'kd': 0.0}
cellRule['secs']['soma']['ions']['ca']={'e': 120}
cellRule['secs']['soma']['pointps']['kleak_0']['gmax']= 3e-5 # 0-0.03 mS/cm^2 for TC
"""
cellRule['secs']['soma']['vinit']=-70
netParams.cellParams['TCrule'] = cellRule

### RE (Destexhe et al., 1996; Bazhenov et al.,2002)
cellRule = netParams.importCellParams(label='RErule', conds={'cellType': 'RE', 'cellModel': 'HH_RE'}, fileName='RE.tem', cellName='sRE')
"""
cellRule['secs']['soma']['mechs']['hh2']={'gnabar': 0.1, 'gkbar': 0.01, 'vtraub': -55.0}
cellRule['secs']['soma']['mechs']['pas']={'g': 5e-5, 'e': -77}
cellRule['secs']['soma']['mechs']['it2']={'shift': 2.0, 'gcabar': 2.3e-3, 'qm': 5.0, 'qh': 3.0}
cellRule['secs']['soma']['mechs']['cad']={'taur': 5.0, 'depth': 1.0, 'kt': 0.0, 'cainf': 2.4e-4, 'kd': 0.0}
cellRule['secs']['soma']['ions']['ca']={'e': 120}
cellRule['secs']['soma']['pointps']['kleak_0']['gmax']= 5e-6  # 0.005 mS/cm^2 for RE
cellRule['secs']['soma']['vinit']=-80

cellRule['secs']['soma']['mechs']['it2']={'shift': 2, 'gcabar': 0.003, 'qm': 2.5, 'qh': 2.5}

cellRule['secs']['soma']['mechs']['itrecustom']['shift']=2
cellRule['secs']['soma']['mechs']['itrecustom']['taubase']=85
cellRule['secs']['soma']['mechs']['itrecustom']['gcabar']=0.003*0
"""
cellRule['secs']['soma']['vinit']=-70
netParams.cellParams['RErule'] = cellRule


###############################################################################
# Synaptic mechanism parameters
###############################################################################
# AMPA
netParams.synMechParams['AMPA'] = {'mod': 'ExpSyn', 'tau': 0.1, 'e': 0}

# AMPA_S
#netParams.synMechParams['AMPA_S'] = {'mod': 'Exp2Syn', 'tau1': 0.05, 'tau2': 5.3, 'e': 0}  # AMPA
netParams.synMechParams['AMPA_S'] = {'mod': 'AMPA_S', 'Cmax': 0.5, 'Cdur': 0.3, 'Alpha': 0.94, 'Beta': 0.18, 'Erev': 0} #}

# NMDA
#netParams.synMechParams['NMDA'] = {'mod': 'Exp2Syn', 'tau1': 0.15, 'tau2': 15, 'e': 0}  # NMDA
#netParams.synMechParams['NMDA'] = {'mod': 'NMDA_S', 'Cdur': 1.0, 'Alpha': 0.11, 'Beta': 0.0066, 'Erev': 0, 'mg': 2} #} # Destexhe, 1998
netParams.synMechParams['NMDA'] = {'mod': 'NMDA_S', 'Cdur': 0.3, 'Alpha': 0.11, 'Beta': 0.0066, 'Erev': 0, 'mg': 2} #} # Destexhe, 1998

# GABAa_S
#netParams.synMechParams['GABAA'] = {'mod': 'Exp2Syn', 'tau1': 0.07, 'tau2': 9.1, 'e': -80}  # GABAA
netParams.synMechParams['GABAA'] = {'mod': 'GABAa_S', 'Cmax': 0.5, 'Cdur': 0.3, 'Alpha': 20, 'Beta': 0.162, 'Erev': -85} # }  # GABAA

# GABAb_S
#netParams.synMechParams['GABAB'] = {'mod': 'Exp2Syn', 'tau1': 0.07, 'tau2': 9.1, 'e': -80}  # GABAB
netParams.synMechParams['GABAB'] = {'mod': 'GABAb_S', 'Cmax': 0.5, 'Cdur': 0.3, 'K1': 0.09, 'K2': 0.0012, 'K3': 0.18, 'K4': 0.034, 'KD': 100, 'Erev': -95} # }  # GABAB
#netParams.synMechParams['GABAB'] = {'mod': 'GABAb_S', 'Cmax': 0.5, 'Cdur': 0.3, 'K1': 0.52, 'K2': 0.0045, 'K3': 0.18, 'K4': 0.034, 'KD': 100, 'Erev': -95} # }  # GABAB

# gap
netParams.synMechParams['GAP'] = {'mod': 'GAP_S', 'r': 1.25e6}

###############################################################################
# Stimulation parameters
###############################################################################
netParams.stimSourceParams['bkg'] = {'type': 'NetStim', 'rate': 1, 'noise': 0.5}
#netParams.stimSourceParams['bkg'] = {'type': 'NetStim', 'number': 10, 'noise': 0.1}

netParams.stimTargetParams['bgCrx'] = {'source': 'bkg', 'conds': {'cellType': ['PY', 'IN']}, 
                                            'weight': 0*2*0.5, 'delay': 'uniform(1,5)', 'synMech': 'AMPA_S'}  
netParams.stimTargetParams['bgThl'] = {'source': 'bkg', 'conds': {'cellType': ['TC']}, 
                                            'weight': 0*2*0.5, 'delay': 'uniform(1,5)', 'synMech': 'AMPA_S'}
#netParams.stimTargetParams['bg->sIN'] = {'source': 'bkg', 'conds': {'cellType': 'IN', 'cellModel': 'HH'}, 
#                                            'weight': 1, 'delay': 'uniform(1,5)', 'synMech': 'AMPA_S'}  

#netParams.stimTargetParams['bg->PYR_HH'] = {'source': 'bkg', 'conds': {'cellType': 'PYR', 'cellModel': 'HH'}, 
#                                            'weight': 1, 'synMech': 'AMPA', 'sec': 'dend', 'loc': 1.0, 'delay': 'uniform(1,5)'}

# IClamp PY
netParams.stimSourceParams['Input_1'] = {'type': 'IClamp', 'del': 10050, 'dur': 100, 'amp': 0.7}
# smallPY=1
netParams.stimTargetParams['Input_1->PY'] = {'source': 'Input_1', 'sec':'soma', 'loc': 0.5, 
                                              'conds': {'pop':'PY', 'cellList': [i*(N_PY/5-1)+11 for i in range((N_PY/20))]}}
# largePY=1
#netParams.stimTargetParams['Input_1->PY'] = {'source': 'Input_1', 'sec':'soma', 'loc': 0.5, 
#                                              'conds': {'pop':'PY', 'cellList': [i*(N_PY/20)+2 for i in range((N_PY/5))]}}
# IClamp RE
netParams.stimSourceParams['Input_2'] = {'type': 'IClamp', 'del': 10050, 'dur': 100, 'amp': -0.1*0}
netParams.stimTargetParams['Input_2->RE'] = {'source': 'Input_2', 'sec':'soma', 'loc': 0.5, 
                                              'conds': {'pop':'RE', 'cellList': [i*(N_RE/20-1)+N_RE/10 for i in range((N_RE/5))]}}
# IClamp TC
netParams.stimSourceParams['Input_3'] = {'type': 'IClamp', 'del': 10050, 'dur': 100, 'amp': -0.1*0}
netParams.stimTargetParams['Input_3->TC'] = {'source': 'Input_3', 'sec':'soma', 'loc': 0.5, 
                                              'conds': {'pop':'TC', 'cellList': [i*(N_TC/20-1)+N_TC/10 for i in range((N_TC/5))]}}
###############################################################################
# Connectivity parameters
###############################################################################

####################### intra cortikal projections ############################
p=0*1.0; pCrx=p; pThl=p; pThlCrx=p # small-world-ness param
#K=0.1 # connectivity param
#intraCrxProb=0.1
PY_PY_AMPA_Prob=0.1;PY_IN_AMPA_Prob=0.1;
PY_PY_NMDA_Prob=0.1;PY_IN_NMDA_Prob=0.1;
IN_PY_GABAA_Prob=0.1;IN_PY_GABAB_Prob=0.1;

gabaapercent=1*0.5*2

"""
netParams.connParams['PY->PY_GAP'] = {
    'preConds': {'popLabel': 'PY'}, 
    'postConds': {'popLabel': 'PY'},
    'weight': 1,
    'delay': 0, 
    'loc': 0.5,
    'synMech': 'GAP_S',
    'connList': smallWorldConn(N_PY,N_PY,pCrx,PY_PY_AMPA_Prob)}  
"""


netParams.connParams['PY->PY_AMPA'] = {
    'preConds': {'popLabel': 'PY'}, 
    'postConds': {'popLabel': 'PY'},
    'weight': 0.6/(N_PY*PY_PY_AMPA_Prob+1),            # (Destexhe, 1998)
    #'weight': 0.6,            # (Destexhe, 1998)
    'delay': netParams.axondelay, 
    'loc': 0.5,
    'synMech': 'AMPA_S',
    #'probability': '1.0 if dist_x <= narrowdiam*xspacing else 0.0'}   
    #'probability': PY_PY_AMPA_Prob}
    'connList': smallWorldConn(N_PY,N_PY,pCrx,PY_PY_AMPA_Prob)}   

netParams.connParams['PY->IN_AMPA'] = {
    'preConds': {'popLabel': 'PY'}, 
    'postConds': {'popLabel': 'IN'},
    'weight': 0.2/(N_IN*PY_IN_AMPA_Prob+1),            # (Destexhe, 1998)       
    #'weight': 0.2,            # (Destexhe, 1998)       
    'delay': netParams.axondelay, 
    'loc': 0.5,
    'synMech': 'AMPA_S',
    #'probability': '1.0 if dist_x <= narrowdiam*xspacing else 0.0'}   
    #'probability': PY_IN_AMPA_Prob}
    'connList': smallWorldConn(N_PY,N_IN,pCrx,PY_IN_AMPA_Prob)}   

netParams.connParams['PY->PY_NMDA'] = {
    'preConds': {'popLabel': 'PY'}, 
    'postConds': {'popLabel': 'PY'},
    'weight': 0*0.25*0.6/(N_PY*PY_PY_NMDA_Prob+1),            # (Destexhe, 1998)       
    #'weight': 0*0.25*0.6,            # (Destexhe, 1998)       
    'delay': netParams.axondelay, 
    'loc': 0.5,
    'synMech': 'NMDA',
    #'probability': '1.0 if dist_x <= narrowdiam*xspacing else 0.0'}   
    #'probability': 0*PY_PY_NMDA_Prob}
    'connList': smallWorldConn(N_PY,N_PY,pCrx,PY_PY_NMDA_Prob)}   

netParams.connParams['PY->IN_NMDA'] = {
    'preConds': {'popLabel': 'PY'}, 
    'postConds': {'popLabel': 'IN'},
    'weight': 0*0.25*0.2/(N_IN*PY_IN_NMDA_Prob+1),            # (Destexhe, 1998)        
    #'weight': 0*0.25*0.2,            # (Destexhe, 1998)        
    'delay': netParams.axondelay, 
    'loc': 0.5,
    'synMech': 'NMDA',
    #'probability': '1.0 if dist_x <= narrowdiam*xspacing else 0.0'}   
    #'probability': 0*PY_IN_NMDA_Prob}
    'connList': smallWorldConn(N_PY,N_IN,pCrx,PY_IN_NMDA_Prob)}   

netParams.connParams['IN->PY_GABAA'] = {
    'preConds': {'popLabel': 'IN'}, 
    'postConds': {'popLabel': 'PY'},
    'weight': gabaapercent*0.15/(N_PY*IN_PY_GABAA_Prob+1),         # (Destexhe, 1998)
    #'weight': gabaapercent*0.15,         # (Destexhe, 1998)
    'delay': netParams.axondelay, 
    'loc': 0.5,
    'synMech': 'GABAA',
    #'probability': '1.0 if dist_x <= narrowdiam*xspacing else 0.0'}   
    #'probability': IN_PY_GABAA_Prob}
    'connList': smallWorldConn(N_IN,N_PY,pCrx,IN_PY_GABAA_Prob)}   

netParams.connParams['IN->PY_GABAB'] = {
    'preConds': {'popLabel': 'IN'}, 
    'postConds': {'popLabel': 'PY'},
    'weight': 0.03/(N_PY*IN_PY_GABAB_Prob+1),         # (Destexhe, 1998)
    #'weight': 0.03,         # (Destexhe, 1998)
    'delay': netParams.axondelay, 
    'loc': 0.5,
    'synMech': 'GABAB',
    #'probability': '1.0 if dist_x <= narrowdiam*xspacing else 0.0'}   
    #'probability': IN_PY_GABAB_Prob}
    'connList': smallWorldConn(N_IN,N_PY,pCrx,IN_PY_GABAB_Prob)}   


###################### intra thalamic projections #############################
intraThlProb=0.1
TC_RE_AMPA_Prob=0.1;RE_TC_GABAA_Prob=0.1;
RE_TC_GABAB_Prob=0.1;RE_RE_GABAA_Prob=0.1;

netParams.connParams['TC->RE'] = {
    'preConds': {'popLabel': 'TC'}, 
    'postConds': {'popLabel': 'RE'},
    'weight': 0.2/(N_RE*TC_RE_AMPA_Prob+1),         # (Destexhe, 1998)  
    #'weight': 0.2,         # (Destexhe, 1998)  
    'delay': 0*netParams.axondelay, 
    'sec': 'soma',
    'loc': 0.5,
    'synMech': 'AMPA_S',
    'threshold': 0,
    #'probability': '1.0 if dist_x <= narrowdiam*xspacing else 0.0'}   
    #'probability': TC_RE_AMPA_Prob}
    'connList': smallWorldConn(N_TC,N_RE,pThl,TC_RE_AMPA_Prob)}   

netParams.connParams['RE->TC_GABAA'] = {
    'preConds': {'popLabel': 'RE'}, 
    'postConds': {'popLabel': 'TC'},
    'weight': 0*0.02/(N_TC*RE_TC_GABAA_Prob+1),         # (Destexhe, 1998)
    #'weight': 0.02,         # (Destexhe, 1998)
    'delay': netParams.axondelay, 
    'loc': 0.5,
    'synMech': 'GABAA',
    #'probability': '1.0 if dist_x <= narrowdiam*xspacing else 0.0'}   
    #'probability': RE_TC_GABAA_Prob}
    'connList': smallWorldConn(N_RE,N_TC,pThl,RE_TC_GABAA_Prob)}   

netParams.connParams['RE->TC_GABAB'] = {
    'preConds': {'popLabel': 'RE'}, 
    'postConds': {'popLabel': 'TC'},
    'weight': 0*0.04/(N_TC*RE_TC_GABAB_Prob+1),         # (Destexhe, 1998)
    #'weight': 0.04,         # (Destexhe, 1998)
    'delay': netParams.axondelay, 
    'loc': 0.5,
    'synMech': 'GABAB',
    #'probability': '1.0 if dist_x <= narrowdiam*xspacing else 0.0'}   
    #'probability': RE_TC_GABAB_Prob}
    'connList': smallWorldConn(N_RE,N_TC,pThl,RE_TC_GABAB_Prob)}   

netParams.connParams['RE->RE'] = {
    'preConds': {'popLabel': 'RE'}, 
    'postConds': {'popLabel': 'RE'},
    'weight': 0*0.2/(N_RE*RE_RE_GABAA_Prob+1),            # (Destexhe, 1998)
    #'weight': 0.2,            # (Destexhe, 1998)
    'delay': netParams.axondelay, 
    'loc': 0.5,
    'synMech': 'GABAA',
    #'probability': '1.0 if dist_x <= narrowdiam*xspacing else 0.0'}   
    #'probability': RE_RE_GABAA_Prob}
    'connList': smallWorldConn(N_RE,N_RE,pThl,RE_RE_GABAA_Prob)}   

################# thalamo-cortical projections ################################
#ThlCrxProb=0.2
PY_TC_AMPA_Prob=0.2;PY_RE_AMPA_Prob=0.2;
TC_PY_AMPA_Prob=0.2;TC_IN_AMPA_Prob=0.2;

netParams.connParams['PY->TC'] = {
    'preConds': {'popLabel': 'PY'}, 
    'postConds': {'popLabel': 'TC'},
    'weight': 0*0.01/(N_TC*PY_TC_AMPA_Prob+1),           # (Destexhe, 1998)    
    #'weight': 0.01,           # (Destexhe, 1998)    
    'delay': netParams.axondelay, 
    'loc': 0.5,
    'synMech': 'AMPA_S',
    #'probability': '1.0 if dist_x <= narrowdiam*xspacing else 0.0'}   
    #'probability': PY_TC_AMPA_Prob}
    'connList': smallWorldConn(N_PY,N_TC,pThlCrx,PY_TC_AMPA_Prob)}   

netParams.connParams['PY->RE'] = {
    'preConds': {'popLabel': 'PY'}, 
    'postConds': {'popLabel': 'RE'},
    'weight': 0*1.2/(N_RE*PY_RE_AMPA_Prob+1),           # (Destexhe, 1998)  
    #'weight': 1.2,           # (Destexhe, 1998)  
    'delay': netParams.axondelay, 
    'loc': 0.5,
    'synMech': 'AMPA_S',
    #'probability': '1.0 if dist_x <= narrowdiam*xspacing else 0.0'}   
    #'probability': PY_RE_AMPA_Prob}
    'connList': smallWorldConn(N_PY,N_RE,pThlCrx,PY_RE_AMPA_Prob)}   

netParams.connParams['TC->PY'] = {
    'preConds': {'popLabel': 'TC'}, 
    'postConds': {'popLabel': 'PY'},
    'weight': 0*1.2/(N_PY*TC_PY_AMPA_Prob+1),        # (Destexhe, 1998)   
    #'weight': 1.2,        # (Destexhe, 1998)   
    'delay': netParams.axondelay, 
    'loc': 0.5,
    'synMech': 'AMPA_S',
    #'probability': '1.0 if dist_x <= narrowdiam*xspacing else 0.0'}   
    #'probability': TC_PY_AMPA_Prob}
    'connList': smallWorldConn(N_TC,N_PY,pThlCrx,TC_PY_AMPA_Prob)}   

netParams.connParams['TC->IN'] = {
    'preConds': {'popLabel': 'TC'}, 
    'postConds': {'popLabel': 'IN'},
    'weight': 0*0.4/(N_IN*TC_IN_AMPA_Prob+1),        # (Destexhe, 1998)  
    #'weight': 0.4,        # (Destexhe, 1998)  
    'delay': netParams.axondelay, 
    'loc': 0.5,
    'synMech': 'AMPA_S',
    #'probability': '1.0 if dist_x <= narrowdiam*xspacing else 0.0'}   
    #'probability': TC_IN_AMPA_Prob}
    'connList': smallWorldConn(N_TC,N_IN,pThlCrx,TC_IN_AMPA_Prob)}   

"""
import matplotlib.pyplot as plt 

# intra-cortical
plt.figure()

connMat=netParams.connParams['PY->PY_AMPA']['connList']
xs = [x[0] for x in connMat]
ys = [x[1] for x in connMat]
plt.subplot(221)
plt.scatter(xs,ys)
plt.title('AMPA')
plt.xlabel('PY')
plt.ylabel('PY')

connMat=netParams.connParams['PY->IN_AMPA']['connList']
xs = [x[0] for x in connMat]
ys = [x[1] for x in connMat]
plt.subplot(222)
plt.scatter(xs,ys)
plt.title('AMPA')
plt.xlabel('PY')
plt.ylabel('IN')

connMat=netParams.connParams['IN->PY_GABAA']['connList']
xs = [x[0] for x in connMat]
ys = [x[1] for x in connMat]
plt.subplot(223)
plt.scatter(xs,ys)
plt.title('GABAA')
plt.xlabel('IN')
plt.ylabel('PY')

connMat=netParams.connParams['IN->PY_GABAB']['connList']
xs = [x[0] for x in connMat]
ys = [x[1] for x in connMat]
plt.subplot(224)
plt.scatter(xs,ys)
plt.title('GABAB')
plt.xlabel('IN')
plt.ylabel('PY')

# thalamo-cortical
plt.figure()

connMat=netParams.connParams['PY->TC']['connList']
xs = [x[0] for x in connMat]
ys = [x[1] for x in connMat]
plt.subplot(221)
plt.scatter(xs,ys)
plt.title('AMPA')
plt.xlabel('PY')
plt.ylabel('TC')

connMat=netParams.connParams['PY->RE']['connList']
xs = [x[0] for x in connMat]
ys = [x[1] for x in connMat]
plt.subplot(222)
plt.scatter(xs,ys)
plt.title('AMPA')
plt.xlabel('PY')
plt.ylabel('RE')

connMat=netParams.connParams['TC->PY']['connList']
xs = [x[0] for x in connMat]
ys = [x[1] for x in connMat]
plt.subplot(223)
plt.scatter(xs,ys)
plt.title('AMPA')
plt.xlabel('TC')
plt.ylabel('PY')

connMat=netParams.connParams['TC->IN']['connList']
xs = [x[0] for x in connMat]
ys = [x[1] for x in connMat]
plt.subplot(224)
plt.scatter(xs,ys)
plt.title('AMPA')
plt.xlabel('TC')
plt.ylabel('IN')

"""

###############################################################################
# SIMULATION PARAMETERS
###############################################################################

#------------------------------------------------------------------------------
# SIMULATION CONFIGURATION
#------------------------------------------------------------------------------

# Simulation parameters
simConfig.checkErrors=True
simConfig.trans = 10000
simConfig.Dt = 0.1
simConfig.steps_per_ms = 1/simConfig.Dt
simConfig.npoints = 30000

simConfig.duration = 13*1000 # simConfig.trans + simConfig.npoints * simConfig.Dt # Duration of the simulation, in ms
simConfig.dt = 0.1 # Internal integration timestep to use
simConfig.hParams['celsius'] = 36
simConfig.hParams['v_init'] = -70
simConfig.seeds = {'conn': 1, 'stim': 1, 'loc': 1} # Seeds for randomizers (connectivity, input stimulation and cell locations)
simConfig.verbose = False  # show detailed messages 

# Recording 
simConfig.recordCells = []  # which cells to record from
simConfig.recordTraces = {'V_soma':{'sec':'soma','loc':0.5,'var':'v'}}

simConfig.recordStim = True  # record spikes of cell stims
simConfig.recordStep = 0.1 # Step size in ms to save data (eg. V traces, LFP, etc)

# Saving
simConfig.simLabel = "knox"
simConfig.saveFolder = "data_knox_v1"
simConfig.filename = 'knox_v1'  # Set file output name
simConfig.saveFileStep = 1000 # step size in ms to save data to disk
#simConfig.savePickle = True # Whether or not to write spikes etc. to a .mat file
#simConfig.saveJson = True
#simConfig.saveMat = True
#simConfig.saveDpk = False

# Analysis and plotting 
#simConfig.analysis['plotRaster'] = {'include': ['PY', 'IN', 'TC', 'RE'], 'orderInverse': True} #True # Whether or not to plot a raster
simConfig.analysis['plotRaster'] = {'include': ['RE', 'TC', 'IN', 'PY'], 'orderInverse': False} #True # Whether or not to plot a raster

#simConfig.analysis['plotRaster'] = True  # Plot raster
simConfig.analysis['plotTraces'] = {'include': [('PY',50),('IN',50),('TC',[50]),('RE',50)]} # plot recorded traces for this list of cells

simConfig.analysis['plotRatePSD'] = {'include': ['PY', 'IN', 'TC', 'RE'], 'Fs': 50, 'smooth': 10} # plot recorded traces for this list of cells

#simConfig.addAnalysis('plot2Dnet', {'include': ['PY', 'IN', 'TC', 'RE'],  'showConns': True, 'saveFig': './images/plot2Dnet.png', 'showFig': False})
#simConfig.addAnalysis('plotShape', {'showSyns': True})
#simConfig.addAnalysis('plotConn', {'include': ['allCells'], 'feature': 'strength'})
#simConfig.analysis.plotConn(include=['allCells'], feature='strength', groupBy='pop', figSize=(9,9), showFig=True)
simConfig.analysis['plotConn'] = True           # plot connectivity matrix

# netParams
simConfig.stimtime = 10050
simConfig.randomstim = 0

simConfig.field = 0
simConfig.fieldg = 0
simConfig.ampafield = 0
simConfig.gabaafield = 0
simConfig.gababfield = 0
simConfig.gababTCfield = 0

simConfig.runStopAt = simConfig.duration

sim.createSimulateAnalyze(netParams = netParams, simConfig = simConfig)

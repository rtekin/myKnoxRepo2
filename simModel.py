# -*- coding: utf-8 -*-
"""
Created on Sat Jul 29 15:57:43 2017

@author: rt
"""

#from cellModels import sINpy
import numpy as np
from neuron import h, gui
from math import sin, cos, pi
from matplotlib import pyplot

from neuronpy.graphics import spikeplot
from neuronpy.util import spiketrain

def ncon(diam, ncells):       # function to get the number of connections 
    nc = 2 * diam + 1
    if nc>ncells: nc = ncells
    return nc

def assign_synapses(gRERE_GABAA,gRETC_GABAA,gRETC_GABAB,gTCRE_AMPA,gPYPY_AMPA,gPYIN_AMPA,gINPY_GABAA,gINPY_GABAB,gPYRE_AMPA,gPYTC_AMPA,gTCPY_AMPA,gTCIN_AMPA): # procedure to assign syn conductances 
    # params: 1=intraRE (RERE), 2=GABAa in TC (RETC),
    # 3=GABAb in TC, 4=AMPA in RE (TCRE)
    # 5=PYPY 6=PYIN 7=GABAa INPY
    # 8=GABAb INPY 9=PYRE 10=PYTC
    # 11=TCPY 12=TCIN
    """
    nRERE = int(RE[0].REgabaalist.count())
    nRETCa = int(TC[0].REgabaalist.count())
    nRETCb = int(TC[0].gababpost.count())
    nTCRE = int(RE[0].TClist.count())
    nPYRE = int(RE[0].PYlist.count())
    nPYTC = int(TC[0].PYlist.count())
    nPYPY = int(PY[0].PYlist.count())
    nPYIN = int(IN[0].PYlist.count())
    nINPYa = int(PY[0].INgabaalist.count())
    nINPYb = int(PY[0].gababpost.count())
    nTCPY = int(PY[0].TClist.count())
    nTCIN = int(IN[0].TClist.count())
    """
    print "nRERE:", nRERE, "nRETCa:", nRETCa, "nRETCb:", nRETCb, "nTCRE:", nTCRE, "nPYPY:", nPYPY, "nINPYa:", nINPYa, "nINPYb:", nINPYb, "nPYRE:", nPYRE, "nPYTC:", nPYTC, "nTCPY:", nTCPY, "nTCIN:", nTCIN
    for i in range(nthalamiccells):
        for j in range(nRERE):
            RE[i].REgabaalist[j].weight[0] = gRERE_GABAA / nRERE               
    
        for j in range(nRETCa):
            TC[i].REgabaalist[j].weight[0] = gRETC_GABAA / nRETCa
    
        for j in range(nRETCb):
            TC[i].gababpost.object(j).gmax = gRETC_GABAB / nRETCb
            #TC[i].REgabablist[j].weight[0] = gRETC_GABAB / nRETCb
        
        for j in range(nTCRE):
            RE[i].TClist[j].weight[0] = gTCRE_AMPA / nTCRE
        
        for j in range(nPYRE):   		
            RE[i].PYlist[j].weight[0] = gPYRE_AMPA / nPYRE
        
        for j in range(nPYTC):
            TC[i].PYlist[j].weight[0] = gPYTC_AMPA / nPYTC
        
    for i in range(ncorticalcells):
        for j in range(nPYPY):
            PY[i].PYlist[j].weight[0] = gPYPY_AMPA / nPYPY
        
        for j in range(nPYIN):
            IN[i].PYlist[j].weight[0] = gPYIN_AMPA / nPYIN
        
        for j in range(nINPYa):
            PY[i].INgabaalist[j].weight[0] = gINPY_GABAA / nINPYa
        
        for j in range(nINPYb):
            PY[i].gababpost[j].gmax = gINPY_GABAB / nINPYb
            #PY[i].INgabablist[j].weight[0] = gINPY_GABAB / nINPYb
        
        for j in range(nTCPY):
            PY[i].TClist[j].weight[0]  = gTCPY_AMPA / nTCPY
        
        for j in range(nTCIN):
            IN[i].TClist[j].weight[0] = gTCIN_AMPA / nTCIN


def printWeight(indx):
    # intra-cortical
    print "PYPY-AMPA_weight = ", PY[indx].PYlist.object(0).weight[0]
    print "PYIN-AMPA_weight = ", IN[indx].PYlist.object(0).weight[0]
    print "INPY-GABAA_weight = ", PY[indx].INgabaalist.object(0).weight[0]
    print "INPY-GABAB_weight = ", PY[indx].gababpost.object(0).gmax
    
    # intra-thalamic
    print "TCRE-AMPA_weight = ", RE[indx].TClist.object(0).weight[0]
    print "RETC-GABAA_weight = ", TC[indx].REgabaalist.object(0).weight[0]
    print "RETC-GABAB_weight = ", TC[indx].gababpost.object(0).gmax
    print "RERE-GABAA_weight = ", RE[indx].REgabaalist.object(0).weight[0]
    
    # thalamo-cortical 
    print "PYTC-AMPA_weight = ", TC[indx].PYlist.object(0).weight[0]
    print "PYRE-AMPA_weight = ", RE[indx].PYlist.object(0).weight[0]
    print "TCPY-AMPA_weight = ", PY[indx].TClist.object(0).weight[0]
    print "TCIN-AMPA_weight = ", IN[indx].TClist.object(0).weight[0]

def printREinfo(indx):
    index=indx
    print " "
    print "------ RE Parameter values (index=",index,") ---------"
    print " "
    
    print "diam=",RE[index].soma.diam,"\t L=",RE[index].soma.L," \t Cm=",RE[index].soma.cm," \t Ra=",RE[index].soma.Ra
    print "g_pas=",RE[index].soma.g_pas," \t e_pas=",RE[index].soma.e_pas," \t vinit="
    print "gnabar_hh2=",RE[index].soma.gnabar_hh2," \t ena=", RE[index].soma.ena
    print "gkbar_hh2=",RE[index].soma.gkbar_hh2," \t ek=",RE[index].soma.ek," \t vtraub_hh2=", RE[index].soma.vtraub_hh2
    print "gcabar_it2=",RE[index].soma.gcabar_it2," \t eca=",RE[index].soma.eca," \t cai=", RE[index].soma.cai," \t cao=", RE[index].soma.cao
    print "shift_it2=",RE[index].soma.shift_it2," \t taubase_it2=",RE[index].soma.taubase_it2," \t qm_it2=", RE[index].soma.qm_it2," \t qh_it2=", RE[index].soma.qh_it2
    print "depth_cad=",RE[index].soma.depth_cad," \t taur_cad=",RE[index].soma.taur_cad," \t cainf_cad=", RE[index].soma.cainf_cad," \t kt_cad=", RE[index].soma.kt_cad
    
    print " "
    print "-------- RE Parameter values (end) --------"
    print " "


ncorticalcells = 100
nthalamiccells = 100
narrowdiam = 5
widediam = 10

watchneuron = 50 
axondelay = 0

smallPY = 1

gabaapercent = 0
gababpercent = 1

celsius = 36
v_init = -70	

randvolt = h.Random()
randvolt.uniform(-80,-65)

recncs = []
#----------------------------------------------------------------------------
#  Create Cells
#----------------------------------------------------------------------------
print " "
print "<<==================================>>"
print "<<            CREATE CELLS          >>"
print "<<==================================>>"
print " "

h.xopen("sPY.tem")      # read geometry file
PY = []     # create PY cells
PYVtrace = [] 
for i in range(ncorticalcells):
    cell = h.sPY()
    #cell.soma.v = randvolt.repick()
    PY.append(cell)

h.xopen("sIN.tem")      # read geometry file
IN = []     # create IN cells
INVtrace = [] 
for i in range(ncorticalcells):
    cell = h.sIN()
    #cell.soma.v = randvolt.repick()
    IN.append(cell)

h.xopen("TC.tem")
TC = []     # create TC cells
TCVtrace = [] 
for i in range(nthalamiccells):
    cell = h.sTC()
    #cell.soma.v = randvolt.repick()
    TC.append(cell)

h.xopen("RE.tem")
RE = []     # create RE cells
REVtrace = [] 
for i in range(nthalamiccells):
    cell = h.sRE()
    #cell.soma.v = randvolt.repick()
    RE.append(cell)


#----------------------------------------------------------------------------
# set up recording stuff
#----------------------------------------------------------------------------

PYtimevec = h.Vector()
PYidvec = h.Vector()
INtimevec = h.Vector()
INidvec = h.Vector()
REtimevec = h.Vector()
REidvec = h.Vector()
TCtimevec = h.Vector()
TCidvec = h.Vector()


######################## intra Cortical synapses #############################
# -----------------------------------------------------
#   Glutamate AMPA receptors in synapses from PY to PY
# -----------------------------------------------------
diamPYPY = narrowdiam		# diameter of connectivity for PY->PY
nPYPY = ncon(diamPYPY, ncorticalcells)	# nb of PY cells receiving synapses from one PY cell

for i in range(ncorticalcells):
    for j in range(i-diamPYPY,i+diamPYPY+1):
        jbound = j
        if jbound < 0: jbound = abs(j) - 1
        if jbound > ncorticalcells-1: jbound = 2 * ncorticalcells - jbound - 1 
        src = PY[i]
        tgt = PY[jbound]
        syn = tgt.ampapostPY #h.AMPA_S(tgt.soma(0.5))
        syn.Alpha = 0.94
        syn.Beta = 0.18
        syn.Cmax = 0.5
        syn.Cdur = 0.3
        syn.Erev = 0
        nc = h.NetCon(src.soma(0.5)._ref_v, syn, 0, axondelay, 1, sec=src.soma)
        #nc.weight[i] = 1
        #nc.delay = axondelay
        #nc.threshold = 0
        #nc.record(PYtimevec, PYidvec, i)
        tgt.PYlist.append(nc)
        #print "[%d,%d], "%(i,jbound)

print " "
print "<< ",int(PY[0].PYlist.count())," AMPA-MEDIATED SYNAPTIC CONTACTS FROM PY TO IN >>"
print " "


#----------------------------------------------------------------------------
#  Glutamate AMPA receptors in synapses from PY to IN
#----------------------------------------------------------------------------

diamPYIN = narrowdiam		# diameter of connectivity for PY->IN
nPYIN = ncon(diamPYIN, ncorticalcells)	

for i in range(ncorticalcells):
    for j in range(i-diamPYIN,i+diamPYIN+1):
        jbound = j
        if jbound < 0: jbound = abs(j) - 1
        if jbound > ncorticalcells-1: jbound = 2 * ncorticalcells - jbound - 1 
        src = PY[i]
        tgt = IN[jbound]
        syn = tgt.ampapost
        syn.Alpha = 0.94
        syn.Beta = 0.18
        syn.Cmax = 0.5
        syn.Cdur = 0.3
        syn.Erev = 0
        nc = h.NetCon(src.soma(0.5)._ref_v, syn, 0, axondelay, 1, sec=src.soma)
        tgt.PYlist.append(nc)
        #print "[%d,%d], "%(i,jbound)

print " "
print "<< ",int(IN[0].PYlist.count())," AMPA-MEDIATED SYNAPTIC CONTACTS FROM PY TO IN >>"
print " "

#--------------------------------------------------
#  GABAa receptors in synapses from IN to PY cells
#--------------------------------------------------

diamINPYa = narrowdiam		# diameter of connectivity from IN->PY
nINPYa = ncon(diamINPYa,ncorticalcells)	

for i in range(ncorticalcells):
    for j in range(i-diamINPYa,i+diamINPYa+1):
        jbound = j
        if jbound < 0: jbound = abs(j) - 1
        if jbound > ncorticalcells-1: jbound = 2 * ncorticalcells - jbound - 1 
        src = IN[i]
        tgt = PY[jbound]
        syn = tgt.gabaapost
        syn.Alpha = 20		# from diffusion model
        syn.Beta = 0.162
        syn.Cmax = 0.5
        syn.Cdur = 0.3
        syn.Erev = -85		# Rinzel's Erev
        nc = h.NetCon(src.soma(0.5)._ref_v, syn, 0, axondelay, 1, sec=src.soma)
        tgt.INgabaalist.append(nc)
        #print "[%d,%d], "%(i,jbound)
	
print " "
print "<< ",int(PY[0].INgabaalist.count())," GABAa-MEDIATED SYNAPTIC CONTACTS FROM IN TO PY >>"
print " "

#--------------------------------------------------
#  GABAb receptors in synapses from IN to PY cells
#--------------------------------------------------
diamINPYb = narrowdiam		# diameter of connectivity from IN->PY
nINPYb = ncon(diamINPYb,ncorticalcells)	

for i in range(ncorticalcells):
    for j in range(i-diamINPYb,i+diamINPYb+1):
        jbound = j
        if jbound < 0: jbound = abs(j) - 1
        if jbound > ncorticalcells-1: jbound = 2 * ncorticalcells - jbound - 1 
        src = IN[i]
        tgt = PY[jbound]
        syn = h.GABAb_S(tgt.soma(0.5))
        syn.K1	= 0.09    #	(/ms mM) forward binding to receptor
        syn.K2	= 0.0012  #	(/ms)	backward (unbinding)of receptor
        syn.K3	= 0.18   # 0.098 # (/ms)	rate of G-protein production
        syn.K4	= 0.034  #	(/ms)	rate of G-protein decay  -  larger number = slower decay?
        syn.KD	= 100	#	dissociation constant of K+ channel
        syn.n	= 4	      # of binding sites of G-protein on K+
        syn.Cmax = 0.5		# short pulses
        syn.Cdur = 0.3
        syn.Erev = -95	#	(mV)	reversal potential (E_K)
        tgt.gababpost.append(syn) 
        nc = h.NetCon(src.soma(0.5)._ref_v, syn, 0, axondelay, 1, sec=src.soma)
        tgt.INgabablist.append(nc)
        #print "[%d,%d], "%(i,jbound)

print " "
print "<< ",int(PY[0].INgabablist.count())," GABAb-MEDIATED SYNAPTIC CONTACTS FROM IN TO PY >>"
print " "


######################## intra Thalamic synapses ############################
#----------------------------------------------------------------------------
#  Glutamate AMPA receptors in synapses from TC to RE
#----------------------------------------------------------------------------

diamTCRE = narrowdiam		# diameter of connectivity for TC->RE
nTCRE = ncon(diamTCRE,nthalamiccells)	# nb of RE cells recampapost, PYlist, TClisteiving synapses from one TC cell

for i in range(nthalamiccells):
    for j in range(i-diamTCRE,i+diamTCRE+1):
        jbound = j
        if jbound < 0: jbound = abs(j) - 1
        if jbound > ncorticalcells-1: jbound = 2 * ncorticalcells - jbound - 1 
        src = TC[i]
        tgt = RE[jbound]
        syn = tgt.ampapost 
        syn.Alpha = 0.94
        syn.Beta = 0.18
        syn.Cmax = 0.5
        syn.Cdur = 0.3
        syn.Erev = 0
        nc = h.NetCon(src.soma(0.5)._ref_v, syn, 0, axondelay, 1, sec=src.soma)
        tgt.TClist.append(nc)
        #print "[%d,%d], "%(i,jbound)

print " "
print "<< ",int(RE[0].TClist.count())," AMPA-MEDIATED SYNAPTIC CONTACTS FROM TC TO RE >>"
print " "

#---------------------------------------
#  GABAa receptors in intra-RE synapses
#---------------------------------------

diamRERE = narrowdiam		# diameter of connectivity for RE->RE
nRERE = ncon(diamRERE,nthalamiccells)	# nb of RE cells receiving synapses from one RE cell

for i in range(nthalamiccells):
    for j in range(i-diamRERE,i+diamRERE+1):
        jbound = j
        if jbound < 0: jbound = abs(j) - 1
        if jbound > ncorticalcells-1: jbound = 2 * ncorticalcells - jbound - 1 
        src = RE[i]
        tgt = RE[jbound]
        syn = tgt.gabaapost 
        syn.Alpha = 20		# from diffusion model
        syn.Beta = 0.162
        syn.Cmax = 0.5		# short pulses
        syn.Cdur = 0.3
        syn.Erev = -85		# Rinzel's Erev
        nc = h.NetCon(src.soma(0.5)._ref_v, syn, 0, axondelay, 1, sec=src.soma)
        tgt.REgabaalist.append(nc)
        #print "[%d,%d], "%(i,jbound)
	
print " "
print "<< ",int(RE[0].REgabaalist.count())," GABAa-MEDIATED SYNAPTIC CONTACTS FROM RE TO RE >>"
print " "

#--------------------------------------------------
#  GABAa receptors in synapses from RE to TC cells
#--------------------------------------------------

diamRETCa = narrowdiam		# diameter of connectivity from RE->TC
nRETCa = ncon(diamRETCa,nthalamiccells)	# nb of RE cells receiving synapses from one TC cell

for i in range(nthalamiccells):
    for j in range(i-diamRETCa,i+diamRETCa+1):
        jbound = j
        if jbound < 0: jbound = abs(j) - 1
        if jbound > ncorticalcells-1: jbound = 2 * ncorticalcells - jbound - 1 
        src = RE[i]
        tgt = TC[jbound]
        syn = tgt.gabaapost 
        syn.Alpha = 20		# from diffusion model
        syn.Beta = 0.162
        syn.Cmax = 0.5		# short pulses
        syn.Cdur = 0.3
        syn.Erev = -85		# Rinzel's Erev
        nc = h.NetCon(src.soma(0.5)._ref_v, syn, 0, axondelay, 1, sec=src.soma)
        tgt.REgabaalist.append(nc)
        #print "[%d,%d], "%(i,jbound)
	
print " "
print "<< ",int(TC[jbound].REgabaalist.count())," GABAa-MEDIATED SYNAPTIC CONTACTS FROM RE TO TC >>"
print " "

#--------------------------------------------------
#  GABAb receptors in synapses from RE to TC cells
#--------------------------------------------------

diamRETCb = narrowdiam		# diameter of connectivity from RE->TC
nRETCb = ncon(diamRETCb,nthalamiccells)	# nb of RE cells receiving synapses from one TC cell

#objectvar gababsyn

for i in range(nthalamiccells):
    for j in range(i-diamRETCb,i+diamRETCb+1):
        jbound = j
        if jbound < 0: jbound = abs(j) - 1
        if jbound > ncorticalcells-1: jbound = 2 * ncorticalcells - jbound - 1 
        src = RE[i]
        tgt = TC[jbound]
        syn = h.GABAb_S(tgt.soma(0.5))
        syn.K1	= 0.09    #	(/ms mM) forward binding to receptor
        syn.K2	= 0.0012  #	(/ms)	backward (unbinding)of receptor
        syn.K3	= 0.18   # 0.098 # (/ms)	rate of G-protein production
        syn.K4	= 0.034  #	(/ms)	rate of G-protein decay  -  larger number = slower decay?
        syn.KD	= 100	#	dissociation constant of K+ channel
        syn.n	= 4	      # of binding sites of G-protein on K+
        syn.Cmax = 0.5		# short pulses
        syn.Cdur = 0.3
        syn.Erev = -95	#	(mV)	reversal potential (E_K)
        tgt.gababpost.append(syn)
        nc = h.NetCon(src.soma(0.5)._ref_v, syn, 0, axondelay, 1, sec=src.soma)
        tgt.REgabablist.append(nc)
        #print "[%d,%d], "%(i,jbound)
	
print " "
print "<< ",int(TC[0].REgabablist.count())," GABAb-MEDIATED SYNAPTIC CONTACTS FROM RE TO TC >>"
print " "

######################## Thalamo-Cortical synapses ############################
#----------------------------------------------------------------------------
#  Glutamate AMPA receptors in synapses from PY to TC
#----------------------------------------------------------------------------

diamPYTC = widediam	
nPYTC = ncon(diamPYTC,nthalamiccells)	

for i in range(ncorticalcells):
    for j in range(i-diamPYTC,i+diamPYTC+1):
        jbound = j
        if jbound < 0: jbound = abs(j) - 1
        if jbound > ncorticalcells-1: jbound = 2 * ncorticalcells - jbound - 1 
        src = PY[i]
        tgt = TC[jbound]
        syn = tgt.ampapost 
        syn.Alpha = 0.94
        syn.Beta = 0.18
        syn.Cmax = 0.5
        syn.Cdur = 0.3
        syn.Erev = 0
        nc = h.NetCon(src.soma(0.5)._ref_v, syn, 0, axondelay, 1, sec=src.soma)
        tgt.PYlist.append(nc)
        #print "[%d,%d], "%(i,jbound)

print " "
print "<< ",int(TC[0].PYlist.count())," AMPA-MEDIATED SYNAPTIC CONTACTS FROM PY TO TC >>"
print " "


#----------------------------------------------------------------------------
#  Glutamate AMPA receptors in synapses from PY to RE
#----------------------------------------------------------------------------

diamPYRE = widediam		
nPYRE = ncon(diamPYRE,nthalamiccells)	
divergence = ncorticalcells / nthalamiccells

for i in range(ncorticalcells):
    for j in range(i/divergence-diamPYRE,i/divergence+diamPYRE+1):
        jbound = j
        if jbound < 0: jbound = abs(j) - 1
        if jbound > ncorticalcells-1: jbound = 2 * ncorticalcells - jbound - 1 
        src = PY[i]
        tgt = RE[jbound]
        syn = tgt.ampapost 
        syn.Alpha = 0.94
        syn.Beta = 0.18
        syn.Cmax = 0.5
        syn.Cdur = 0.3
        syn.Erev = 0
        nc = h.NetCon(src.soma(0.5)._ref_v, syn, 0, axondelay, 1, sec=src.soma)
        tgt.PYlist.append(nc)
        #print "[%d,%d], "%(i,jbound)

print " "
print "<< ",int(RE[0].PYlist.count())," AMPA-MEDIATED SYNAPTIC CONTACTS FROM PY TO RE >>"
print " "

#----------------------------------------------------------------------------
#  Glutamate AMPA receptors in synapses from TC to PY
#----------------------------------------------------------------------------

diamTCPY = widediam		
nTCPY = ncon(diamTCPY,ncorticalcells)	
divergence = ncorticalcells / nthalamiccells

for i in range(nthalamiccells):
    for j in np.arange(divergence*(i-diamTCPY),divergence*(i+diamTCPY+1),divergence):
        jbound = j
        if jbound < 0: jbound = abs(j) - divergence
        if jbound > ncorticalcells-1: jbound = 2 * ncorticalcells - jbound - divergence 
        # presynaptic is TC[i], postsynaptic is PY[j] 
        src = TC[i]
        tgt = PY[jbound+divergence/2]
        syn = tgt.ampapostTC 
        syn.Alpha = 0.94
        syn.Beta = 0.18
        syn.Cmax = 0.5
        syn.Cdur = 0.3
        syn.Erev = 0
        nc = h.NetCon(src.soma(0.5)._ref_v, syn, 0, axondelay, 1, sec=src.soma)
        tgt.TClist.append(nc)
        #print "[%d,%d], "%(i,jbound)

print " "
print "<< ",int(PY[0].TClist.count())," AMPA-MEDIATED SYNAPTIC CONTACTS FROM TC TO PY >>"
print " "

#----------------------------------------------------------------------------
#  Glutamate AMPA receptors in synapses from TC to IN
#----------------------------------------------------------------------------

diamTCIN = widediam		# diameter of connectivity for TC->IN
nTCIN = ncon(diamTCIN,ncorticalcells)	

divergence = ncorticalcells / nthalamiccells

for i in range(nthalamiccells):
    for j in np.arange(divergence*(i-diamTCIN),divergence*(i+diamTCIN+1),divergence):
        jbound = j
        if jbound < 0: jbound = abs(j) - divergence
        if jbound > ncorticalcells-1: jbound = 2 * ncorticalcells - jbound - divergence 
        # presynaptic is TC[i], postsynaptic is IN[j] 
        src = TC[i]
        tgt = IN[jbound+divergence/2]
        syn = tgt.ampapost 
        syn.Alpha = 0.94
        syn.Beta = 0.18
        syn.Cmax = 0.5
        syn.Cdur = 0.3
        syn.Erev = 0
        nc = h.NetCon(src.soma(0.5)._ref_v, syn, 0, axondelay, 1, sec=src.soma)
        tgt.TClist.append(nc)
        #print "[%d,%d], "%(i,jbound)

print " "
print "<< ",int(IN[0].TClist.count())," AMPA-MEDIATED SYNAPTIC CONTACTS FROM TC TO IN >>"
print " "


#---------------------------------------------------------------------------
#   Assign potassium leak and Ih current strengths in TC cells 
#---------------------------------------------------------------------------

rgh = h.Random()
rk1 = h.Random()
rgh.normal(17.5,0.0008)        #random number generator behaves weirdly for very small numbers, so multiply by 10^-6 below
rk1.normal(40,0.003)

# setup TC cells in resting mode (no spontaneous oscillation)
for i in range(nthalamiccells):
    TC[i].soma.ghbar_iar = rgh.repick() * 1e-6   
    TC[i].kl.gmax = rk1.repick() * 1e-4         
    print "TC(",i,") gh:",TC[i].soma.ghbar_iar," gmax:",TC[i].kl.gmax

#----------------------------------------------------------------------------
#  set up current stimulus to cells
#----------------------------------------------------------------------------
"""
PYstim = []
stimtime = 500

if smallPY==1:
    #for i in range(ncorticalcells/20):
        #stim = h.IClamp(PY[i*(ncorticalcells/5)+9].soma(0.5))
    for i in range(ncorticalcells/20):
        stim = h.IClamp(0.5, sec=PY[(i*(ncorticalcells/5-1)+11)*0+i].soma)
        stim.delay = stimtime
        stim.dur = 100
        stim.amp = 0.7*5
        PYstim.append(stim)

"""    
stim = h.NetStim() # Make a new stimulator
syn_ = h.ExpSyn(PY[0].soma(0.5))
stim.number = 1
stim.start = 100
ncstim = h.NetCon(stim, syn_)
ncstim.delay = 1
ncstim.weight[0] = 1.5*3 # NetCon weight is a vector.
syn_.tau = 0.1


h("objref nil")
for i in range(ncorticalcells):
    VtracePY = h.Vector()    # Membrane potential vector at soma
    VtracePY.record(PY[i].soma(0.5)._ref_v)
    PYVtrace.append(VtracePY)
    nc = h.NetCon(PY[i].soma(0.5)._ref_v, h.nil, sec=PY[i].soma)
    nc.record(PYtimevec, PYidvec, i+1)
    recncs.append(nc)
    #
    VtraceIN = h.Vector()
    VtraceIN.record(IN[i].soma(0.5)._ref_v)
    INVtrace.append(VtraceIN)
    nc = h.NetCon(IN[i].soma(0.5)._ref_v, h.nil, sec=IN[i].soma)
    nc.record(INtimevec, INidvec, i+1)
    recncs.append(nc)

for i in range(nthalamiccells):
    VtraceTC = h.Vector()
    VtraceTC.record(TC[i].soma(0.5)._ref_v)
    TCVtrace.append(VtraceTC)
    nc = h.NetCon(TC[i].soma(0.5)._ref_v, h.nil, sec=TC[i].soma)
    nc.record(TCtimevec, TCidvec, i+1)
    recncs.append(nc)
    #
    VtraceRE = h.Vector()
    VtraceRE.record(RE[i].soma(0.5)._ref_v)
    REVtrace.append(VtraceRE)
    nc = h.NetCon(RE[i].soma(0.5)._ref_v, h.nil, sec=RE[i].soma)
    nc.record(REtimevec, REidvec, i+1)
    recncs.append(nc)


t_vec = h.Vector()        # Time stamp vector
t_vec.record(h._ref_t)
#syn_i_vec = h.Vector()
#syn_i_vec.record(syn_._ref_i)

h.dt = 0.01
#h.tstop = 1000
tstop = 1000


# run the simulation

def go():
    printREinfo(0)
    
    h.finitialize(v_init)
    
    # params:       RERE, RETCa, RETCb, TCRE, PYPY, PYIN, INPYa,             INPYb,              PYRE,   PYTC,   TCPY,   TCIN
    assign_synapses(0.2,  0.02,  0.04,  0.2,  0.6,  0.2,  gabaapercent*0.15, gababpercent*0.03,  0*1.2,  0*0.01, 0*1.2,  0*0.4)
    printWeight(0)
    
    h.fcurrent()
    h.cvode.re_init()
    h.frecord_init()
    
    while h.t<tstop:
        h.fadvance()

go()

#h.run()

#------------------------------------------
# figures
#------------------------------------------
# voltage traces
fig = pyplot.figure(figsize=(8,4))
ax1 = fig.add_subplot(1,1,1)
PY1 = ax1.plot(t_vec, PYVtrace[0], color='black')
PY2 = ax1.plot(t_vec, PYVtrace[11], color='red')
PY3 = ax1.plot(t_vec, PYVtrace[30], color='blue')
PY4 = ax1.plot(t_vec, PYVtrace[49], color='green')
ax1.legend(PY1 + PY2 + PY3 + PY4, ['PY1', 'PY2', 'PY%d'%(ncorticalcells/3), 'PY%d'%(ncorticalcells/2)])
ax1.set_ylabel('mV')
#ax1.set_xticks([]) # Use ax2's tick labels
ax1.set_xlabel('time (ms)')
"""
ax2 = fig.add_subplot(2,1,2)
syn_plot = ax2.plot(t_vec, syn_i_vec, color='blue')
ax2.legend(syn_plot, ['synaptic current'])
ax2.set_ylabel(h.units('ExpSyn.i'))
ax2.set_xlabel('time (ms)')
pyplot.show()
"""

fig = pyplot.figure(figsize=(8,4))
ax1 = fig.add_subplot(1,1,1)
IN1 = ax1.plot(t_vec, INVtrace[0], color='black')
IN2 = ax1.plot(t_vec, INVtrace[1], color='red')
IN3 = ax1.plot(t_vec, INVtrace[ncorticalcells/3-1], color='blue')
IN4 = ax1.plot(t_vec, INVtrace[ncorticalcells/2-1], color='green')
ax1.legend(IN1 + IN2 + IN3 + IN4, ['IN1', 'IN2', 'IN%d'%(ncorticalcells/3), 'IN%d'%(ncorticalcells/2)])
ax1.set_ylabel('mV')
#ax1.set_xticks([]) # Use ax2's tick labels
ax1.set_xlabel('time (ms)')

fig = pyplot.figure(figsize=(8,4))
ax1 = fig.add_subplot(1,1,1)
TC1 = ax1.plot(t_vec, TCVtrace[0], color='black')
TC2 = ax1.plot(t_vec, TCVtrace[1], color='red')
TC3 = ax1.plot(t_vec, TCVtrace[nthalamiccells/3-1], color='blue')
TC4 = ax1.plot(t_vec, TCVtrace[nthalamiccells/2-1], color='green')
ax1.legend(TC1 + TC2 + TC3 + TC4, ['TC1', 'TC2', 'TC%d'%(nthalamiccells/3), 'TC%d'%(nthalamiccells/2)])
ax1.set_ylabel('mV')
#ax1.set_xticks([]) # Use ax2's tick labels
ax1.set_xlabel('time (ms)')

fig = pyplot.figure(figsize=(8,4))
ax1 = fig.add_subplot(1,1,1)
RE1 = ax1.plot(t_vec, REVtrace[0], color='black')
RE2 = ax1.plot(t_vec, REVtrace[1], color='red')
RE3 = ax1.plot(t_vec, REVtrace[nthalamiccells/3-1], color='blue')
RE4 = ax1.plot(t_vec, REVtrace[nthalamiccells/2-1], color='green')
ax1.legend(RE1 + RE2 + RE3 + RE4, ['RE1', 'RE2', 'RE%d'%(nthalamiccells/3), 'RE%d'%(nthalamiccells/2)])
ax1.set_ylabel('mV')
#ax1.set_xticks([]) # Use ax2's tick labels
ax1.set_xlabel('time (ms)')

#----------------------------------------------------------------------------
# Raster plot
#----------------------------------------------------------------------------
#sp = spikeplot.SpikePlot(spikes,label='spikes1',savefig=True,marker='.',markerscale=1)
#sp.plot_spikes(spikes)
sp = spikeplot.SpikePlot()
sp.set_markerscale(1)
sp.set_marker('.')

PYspikes = spiketrain.netconvecs_to_listoflists(PYtimevec, PYidvec)
sp.set_markercolor('black')
sp.plot_spikes(PYspikes, draw=False, label='spikes')

INspikes = spiketrain.netconvecs_to_listoflists(INtimevec, INidvec)
sp.set_markercolor('blue')
sp.plot_spikes(INspikes, cell_offset=len(PYspikes))

TCspikes = spiketrain.netconvecs_to_listoflists(TCtimevec, TCidvec)
sp.set_markercolor('red')
sp.plot_spikes(TCspikes, cell_offset=len(INspikes))

REspikes = spiketrain.netconvecs_to_listoflists(REtimevec, REidvec)
sp.set_markercolor('green')
sp.plot_spikes(REspikes, cell_offset=len(TCspikes))



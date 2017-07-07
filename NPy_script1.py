# -*- coding: utf-8 -*-
"""
Created on Wed Jul 05 17:13:12 2017

@author: rttr
"""

import neuron
#from neuron import h

ncorticalcells = 100
nthalamiccells = 100
narrowdiam = 5
widediam = 10

trans = 10000*0
Dt = 0.1
dt = 0.1			# must be submultiple of Dt
npoints = 30000+10000*10
stimtime = 10050			
randomstim = 0

fieldlower = 30
fieldupper = 70
fielddist = 50

watchneuron = 50 
axondelay = 0

smallPY = 1
mediumPY = 0
largePY = 0
largePYhole = 0
smallPYhole = 0
mediumPYoffset = 7
largePYoffset = 2

gabaapercent = 1*0.5*2
gababpercent = 1

#h.load_file('TC.tem')
#h('objectvar PY[ncorticalcells]')
#d=h.objectvar()
#PY[ncorticalcells] = h.objectvar()


neuron.hoc.execute('load_file("Fspikewave.oc")')
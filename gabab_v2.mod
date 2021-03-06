: $Id: gababS.mod,v 1.4 2004/06/17 18:38:24 billl Exp $
TITLE simple GABAb receptors

COMMENT
-----------------------------------------------------------------------------

	Kinetic model of GABA-B receptors
	=================================

  MODEL OF SECOND-ORDER G-PROTEIN TRANSDUCTION AND FAST K+ OPENING
  WITH COOPERATIVITY OF G-PROTEIN BINDING TO K+ CHANNEL

  PULSE OF TRANSMITTER

  SIMPLE KINETICS WITH NO DESENSITIZATION

	Features:

  	  - peak at 100 ms; time course fit to Tom Otis' PSC
	  - SUMMATION (psc is much stronger with bursts)


	Approximations:

	  - single binding site on receptor	
	  - model of alpha G-protein activation (direct) of K+ channel
	  - G-protein dynamics is second-order; simplified as follows:
		- saturating receptor
		- no desensitization
		- Michaelis-Menten of receptor for G-protein production
		- "resting" G-protein is in excess
		- Quasi-stat of intermediate enzymatic forms
	  - binding on K+ channel is fast


	Kinetic Equations:

	  dR/dt = K1 * T * (1-R-D) - K2 * R

	  dG/dt = K3 * R - K4 * G

	  R : activated receptor
	  T : transmitter
	  G : activated G-protein
	  K1,K2,K3,K4 = kinetic rate cst

  n activated G-protein bind to a K+ channel:

	n G + C <-> O		(Alpha,Beta)

  If the binding is fast, the fraction of open channels is given by:

	O = G^n / ( G^n + KD )

  where KD = Beta / Alpha is the dissociation constant

-----------------------------------------------------------------------------

  Parameters estimated from patch clamp recordings of GABAB PSP's in
  rat hippocampal slices (Otis et al, J. Physiol. 463: 391-407, 1993).

-----------------------------------------------------------------------------

  PULSE MECHANISM

  Kinetic synapse with release mechanism as a pulse.  

  Warning: for this mechanism to be equivalent to the model with diffusion 
  of transmitter, small pulses must be used...

  see details at http://cns.iaf.cnrs-gif.fr

  Written by A. Destexhe, 1995
  27-11-2002: the pulse is implemented using a counter, which is more
	stable numerically (thanks to Yann LeFranc)

-----------------------------------------------------------------------------
ENDCOMMENT



INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	POINT_PROCESS GABAb_S_v2
	POINTER pre
	RANGE C, R, G, g, gmax, lastrelease, TimeCount
	NONSPECIFIC_CURRENT i
	:GLOBAL Cmax, Cdur, Prethresh, Deadtime
	RANGE Cmax, Cdur, Prethresh, Deadtime
	:GLOBAL K1, K2, K3, K4, KD, Erev
	RANGE K1, K2, K3, K4, KD, n, Erev
}
UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(umho) = (micromho)
	(mM) = (milli/liter)
}

PARAMETER {
	dt		(ms)
	Cmax	= 0.5	(mM)		: max transmitter concentration
	Cdur	= 0.3	(ms)		: transmitter duration (rising phase)
	Prethresh = 0 			: voltage level nec for release
	Deadtime = 1	(ms)		: mimimum time between release events
:
:	From Kfit with long pulse (5ms 0.5mM)
:
	K1	= 0.52	(/ms mM)	: forward binding rate to receptor
	K2	= 0.0013 (/ms)		: backward (unbinding) rate of receptor
	K3	= 0.098 (/ms)		: rate of G-protein production
	K4	= 0.033 (/ms)		: rate of G-protein decay
	KD	= 100			: dissociation constant of K+ channel
	n	= 4			: nb of binding sites of G-protein on K+
	Erev	= -95	(mV)		: reversal potential (E_K)
	gmax		(umho)		: maximum conductance
}

ASSIGNED {
	v		(mV)		: postsynaptic voltage
	i 		(nA)		: current = g*(v - Erev)
	g 		(umho)		: conductance
	C		(mM)		: transmitter concentration
	Gn
	pre 				: pointer to presynaptic variable
	lastrelease	(ms)		: time of last spike
	TimeCount	(ms)		: time counter
}

STATE {
	R				: fraction of activated receptor
	G				: fraction of activated G-protein
}


INITIAL {
	C = 0
	lastrelease = -9e9

	R = 0
	G = 0
	TimeCount=-1
}

BREAKPOINT {
	SOLVE bindkin METHOD cnexp
	Gn = G^n
	g = gmax * Gn / (Gn+KD)
	i = g*(v - Erev)
}

DERIVATIVE bindkin {
	R' = K1 * C * (1-R) - K2 * R
	G' = K3 * R - K4 * G
}

NET_RECEIVE(weight ,on) {
  if (on) { : at end of Cdur pulse so turn off
    on = 0
    C = 0
  } else {
    if (weight==0) { : do nothing
      C=0 
    } else {
      C=Cmax
      on=1
      net_send(Cdur, on)
    }
  }
}







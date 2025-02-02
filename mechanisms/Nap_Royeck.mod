TITLE nap

COMMENT
	
	File copied from Royeck et al., 2008, ModelDB #115356
	Converted into pS for easier implementation
	half changed from -52.3 to -62.3

	Made threadsafe (CCohen)

ENDCOMMENT

NEURON {
	
	SUFFIX nap_roy
	USEION na READ ena WRITE ina
	RANGE  gbar, g, timestauh, timestaum, shifttaum, shifttauh
	GLOBAL minf, mtau	
	THREADSAFE
}

PARAMETER {
	
	:gbar = 5   	(pS/um2)
	gbar = 1 (S/cm2)
	
	:q10m=3.1
	:q10h=2.3
	
	timestauh=1
	timestaum=1
	shifttaum=1
	shifttauh=1
	
	
	ena	= 55	(mV)		: Golomb et al. //need to check RMP
	:ena		(mV)		: must be explicitly def. in hoc
	celsius 	(degC)
	v 			(mV)
}


UNITS {
	
	(mA) = (milliamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(um) = (micron)
} 

ASSIGNED {
	g 			(S/cm2)
	ina 		(mA/cm2)
	:inap		(mA/cm2)
	minf 		:hinf 		
	mtau (ms)	:htau (ms) 	
}

STATE { m }

UNITSOFF

BREAKPOINT {
    
    SOLVE states METHOD cnexp
	mtau = 1
	minf = (1/(1+exp(-(v+62.3)/6.8))) :midpoint -56.3, slope 7.4	        		
	g = gbar * m
	ina = g * (v - ena)
} 

INITIAL {
	
	mtau = 1
	minf = (1/(1+exp(-(v+62.3)/6.8))) :midpoint - 52.3 slope 6.8	5.5
	m=minf  
}

DERIVATIVE states {   
  
	mtau = 1
	minf = (1/(1+exp(-(v+62.3)/6.8))) :midpoint - 52.3 (47) slope 6.8	         	
	m' = (minf-m)/mtau
}

UNITSON

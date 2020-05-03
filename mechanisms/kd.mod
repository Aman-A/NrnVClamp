TITLE K-D
: K-D current for prefrontal cortical neuron ------Yuguo Yu  2007

NEURON {
	SUFFIX kd
	USEION k READ ek WRITE ik
	RANGE  gkbar, ik, ek, gk
	GLOBAL minf, mtau, hinf, htau, q10
}

PARAMETER {
	gkbar = 1   	(S/cm2)										
	ek = -100	(mV)            : must be explicitly def. in hoc
	v 		(mV)
	vhalfm=-43  (mV)
	km=8
	vhalfh=-67  (mV) 
      kh=7.3
	q10=2.3
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(um) = (micron)
} 

ASSIGNED {
	ik 		(mA/cm2)
	gk (S/cm2)
	minf 		mtau (ms)	 	
	hinf 		htau (ms)	 	
}
 

STATE { m h}

BREAKPOINT {
	SOLVE states METHOD cnexp
	gk = gkbar * m*h
	ik = gk*(v-ek)
} 

INITIAL {
	trates(v)
	m=minf  
	h=hinf  
}

DERIVATIVE states {   
        trates(v)      
        m' = (minf-m)/mtau
        h' = (hinf-h)/htau
}

PROCEDURE trates(v) {  
	LOCAL qt
        qt=q10^((celsius-22)/10)
        minf=1-1/(1+exp((v-vhalfm)/km))
        hinf=1/(1+exp((v-vhalfh)/kh))

  	 mtau = 0.6
	 htau = 1500
}


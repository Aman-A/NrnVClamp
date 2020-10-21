TITLE Axonal Kv7-current

COMMENT

Model of cortical pyramidal neuron Kv7/M currents. The kinetic parameters, voltage-dependence and reversal potentials are estimated by fitting and analysis of axonal K7/M-currents at 33 degrees celsius. MHP Kole, Canberra, 2008 and Amsterdam, 2011.

Made threadsafe (CCohen)

Modified by Aman Aberra

ENDCOMMENT

UNITS {
	
	(mA) = (milliamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(um) = (micron)
}

:INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

PARAMETER {	
	
	dt	(ms)
	v 	(mV)
	celsius	(degC)
	Ca	=	0.0388
	Cb	= 	0.00168
	za	=	0.90979
	zb	=	1.23645
	:gbar	= 	20 	(pS/um2)	: 0.002 mho/cm2
	gbar	= 	0.000002 (S/cm2)	: 0.002 mho/cm2
	temp	= 	34	(degC)		: original temp
	q10  	= 	3.0				: temperature sensitivity
}

NEURON {
	
	SUFFIX Kv7
	USEION k READ ek WRITE ik
	RANGE gbar, g, ik
	GLOBAL q10
	THREADSAFE
}

STATE { m }

ASSIGNED {
	
	ik (mA/cm2)
	g (S/cm2)
	ek (mV)
	tadj
}

INITIAL {
	
	m=alpha(v)/(beta(v)+alpha(v))
	tadj = q10^((celsius - temp)/10)
}

BREAKPOINT {
	
	SOLVE state METHOD cnexp
	:tadj = q10^((celsius - temp)/10)	:this repeated calculation allows changes in temperature during the simulation
	:ik = (1e-4) * gbar * m * (v-ek)
	g = gbar * m
	ik = g * (v-ek)
}

FUNCTION alpha(v(mV)) {
	
	alpha = tadj*Ca*exp(za*v*0.037788)
}

FUNCTION beta(v(mV)) {
	
	beta = tadj*Cb*exp(-zb*v*0.037788)
}

DERIVATIVE state {
	
	m' = (1-m)*alpha(v) - m*beta(v)
}

:Reference :Colbert and Pan 2002
:comment: took the NaTa and shifted both activation/inactivation by 6 mv

NEURON	{
	SUFFIX NaTs2_t
	USEION na READ ena WRITE ina
	RANGE gNaTs2_tbar, gNaTs2_t, ina
	GLOBAL q10
}

UNITS	{
	(S) = (siemens)
	(mV) = (millivolt)
	(mA) = (milliamp)
}

PARAMETER	{
	gNaTs2_tbar = 1 (S/cm2)
	q10  = 2.3			: temperature sensitivity
	vShift = 0 (mV) : try -8 mV to match Schmidt-Hieber 2010
}

ASSIGNED	{
	v	(mV)
	ena	(mV)
	ina	(mA/cm2)
	gNaTs2_t	(S/cm2)
	mInf
	mTau
	mAlpha
	mBeta
	hInf
	hTau
	hAlpha
	hBeta
}

STATE	{
	m
	h
}

BREAKPOINT	{
	SOLVE states METHOD cnexp
	gNaTs2_t = gNaTs2_tbar*m*m*m*h
	ina = gNaTs2_t*(v-ena)
}

DERIVATIVE states	{
	rates(q10)
	m' = (mInf-m)/mTau
	h' = (hInf-h)/hTau
}

INITIAL{
	rates(q10)
	m = mInf
	h = hInf
}

PROCEDURE rates(q10){
  LOCAL qt
  :qt = 2.3^((34-21)/10)
  qt = q10^((celsius-21)/10)	
	UNITSOFF
    if(v-vShift == -32){
    	v = v+0.0001
    }
		mAlpha = (0.182 * ((v-vShift)- -32))/(1-(exp(-((v-vShift)- -32)/6)))
		mBeta  = (0.124 * (-(v-vShift) -32))/(1-(exp(-(-(v-vShift) -32)/6)))
		mInf = mAlpha/(mAlpha + mBeta)
		mTau = (1/(mAlpha + mBeta))/qt

    if(v-vShift == -60){
      v = v + 0.0001
    }
		hAlpha = (-0.015 * ((v-vShift)- -60))/(1-(exp(((v-vShift)- -60)/6)))
		hBeta  = (-0.015 * (-(v-vShift) -60))/(1-(exp((-(v-vShift) -60)/6)))
		hInf = hAlpha/(hAlpha + hBeta)
		hTau = (1/(hAlpha + hBeta))/qt
	UNITSON
}
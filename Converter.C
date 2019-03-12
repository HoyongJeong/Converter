#include "Converter.h"

#include "TMath.h"


//////////////////////////////////////////////////
//   Default convstructor                       //
//////////////////////////////////////////////////
Converter::Converter()
{
	m_MPri = 0.0;
	m_MTar = 0.0;
}


//////////////////////////////////////////////////
//   Convstructor 2                             //
//////////////////////////////////////////////////
Converter::Converter(double pri, double tar)
{
	m_MPri = pri;
	m_MTar = tar;
}


//////////////////////////////////////////////////
//   Destructor                                 //
//////////////////////////////////////////////////
Converter::~Converter()
{
}


//////////////////////////////////////////////////
//   Set primary mass                           //
//////////////////////////////////////////////////
void Converter::SetMPri(double pri)
{
	m_MPri = pri;
}


//////////////////////////////////////////////////
//   Set primary mass                           //
//////////////////////////////////////////////////
void Converter::SetMTar(double tar)
{
	m_MTar = tar;
}


//////////////////////////////////////////////////
//   Convert CM angle to momentum transfer      //
//////////////////////////////////////////////////
void Converter::ThetaCMToQ(double &q, double &eq, double beamE, double thetaCM, double eBeamE = 0.0, double eThetaCM = 0.0)
{
	// Deg to rad
	thetaCM  *= TMath::DegToRad();
	eThetaCM *= TMath::DegToRad();

	// Momentum of primary in CM
	double p = m_MTar
	         * TMath::Sqrt( (2.0*m_MPri*beamE+beamE*beamE)
			              / ((m_MPri+m_MTar)*(m_MPri+m_MTar)+2.0*m_MTar*beamE) );

	q = 2.0 * p * TMath::Sin(0.5*thetaCM);

	double dqdE = 2.0 * TMath::Sin(0.5*thetaCM)
	            * m_MTar * m_MTar / p
				* ((m_MPri+beamE)*(m_MPri+m_MTar)*(m_MPri+m_MTar)+2*m_MTar*beamE*beamE)
				/ ((m_MPri)*(m_MTar)+2.0*m_MTar*beamE)
				/ ((m_MPri)*(m_MTar)+2.0*m_MTar*beamE);
	double dqdT = p * TMath::Cos(0.5*thetaCM);

	eq = TMath::Sqrt(dqdE*dqdE*eBeamE*eBeamE + dqdT*dqdT*eThetaCM*eThetaCM);
}


//////////////////////////////////////////////////
//   Convert CM angle to momentum transfer      //
// squared                                      //
//////////////////////////////////////////////////
void Converter::ThetaCMToT(double &t, double &et, double beamE, double thetaCM, double eBeamE = 0.0, double eThetaCM = 0.0)
{
	// Deg to rad
	thetaCM  *= TMath::DegToRad();
	eThetaCM *= TMath::DegToRad();

	// Momentum of primary in CM
	double p = m_MTar
	         * TMath::Sqrt( (2.0*m_MPri*beamE+beamE*beamE)
			              / ((m_MPri+m_MTar)*(m_MPri+m_MTar)+2.0*m_MTar*beamE) );
	double q = 2.0 * p * TMath::Sin(0.5*thetaCM);

	t = - q * q;

	double dtdE = 8.0 * p * TMath::Sin(0.5*thetaCM) * TMath::Sin(0.5*thetaCM)
	            * m_MTar * m_MTar / p
				* ((m_MPri+beamE)*(m_MPri+m_MTar)*(m_MPri+m_MTar)+2*m_MTar*beamE*beamE)
				/ ((m_MPri)*(m_MTar)+2.0*m_MTar*beamE)
				/ ((m_MPri)*(m_MTar)+2.0*m_MTar*beamE);
	double dtdT = 2.0 * p * p * TMath::Sin(thetaCM);

	et = TMath::Sqrt(dtdE*dtdE*eBeamE*eBeamE + dtdT*dtdT*eThetaCM*eThetaCM);
}


//////////////////////////////////////////////////
//   Convert CM angle to Lab angle              //
//////////////////////////////////////////////////
void Converter::ThetaCMToThetaLab(double &tL, double &etL, double beamE, double thetaCM, double eBeamE = 0.0, double eThetaCM = 0.0)
{
	// Deg to rad
	thetaCM  *= TMath::DegToRad();
	eThetaCM *= TMath::DegToRad();

	double A = TMath::Sqrt((m_MPri+m_MTar)*(m_MPri+m_MTar) + 2.0*m_MTar*beamE)
	         / (m_MPri + m_MTar + beamE);
	double B = (m_MPri*m_MPri + m_MPri*m_MTar + m_MTar*beamE)
	         / m_MTar / (m_MPri + m_MTar + beamE);

	double sint = TMath::Sin(thetaCM);
	double cost = TMath::Cos(thetaCM);


	double tan = A * sint / (cost + B);
	double aTan = TMath::ATan(tan);
	
	if ( aTan < 0.0 ) aTan += TMath::Pi();

	tL = aTan;


	double dAdE = - (m_MPri*m_MPri + m_MPri*m_MTar + m_MTar*beamE)
	            / (m_MPri + m_MTar + beamE) / (m_MPri + m_MTar + beamE)
		        / TMath::Sqrt((m_MPri+m_MTar)*(m_MPri+m_MTar) + 2.0*m_MTar*beamE);
	double dBdE = (m_MTar*m_MTar - m_MPri*m_MPri) / m_MTar
	            / (m_MPri + m_MTar + beamE) / (m_MPri + m_MTar + beamE);

	double dtdE = (dAdE*sint*(cost+B) - dBdE*A*sint)
	            / TMath::Cos(tan) / TMath::Cos(tan)
				/ (cost + B) / (cost + B);
	double dtdt = A * (cost*(cost+B) - sint*sint)
	            / TMath::Cos(tan) / TMath::Cos(tan)
				/ (cost + B) / (cost + B);

	etL = TMath::Sqrt(dtdE*dtdE*eBeamE*eBeamE + dtdt*dtdt*eThetaCM*eThetaCM);

	// Rad to deg
	tL  *= TMath::RadToDeg();
	etL *= TMath::RadToDeg();
}


//////////////////////////////////////////////////
//   Convert momentum transfer to CM angle      //
//////////////////////////////////////////////////
void Converter::QToThetaCM(double &tC, double &etC, double beamE, double q, double eBeamE = 0.0, double eq = 0.0)
{
	// Momentum of primary in CM
	double p = m_MTar
	         * TMath::Sqrt( (2.0*m_MPri*beamE+beamE*beamE)
			              / ((m_MPri+m_MTar)*(m_MPri+m_MTar)+2.0*m_MTar*beamE) );

	tC = 2.0 * TMath::ASin(0.5 * q / p);


	double dCdq = 1.0 / TMath::Cos(0.5*tC);
	double dCdE = q * dCdq / p;

	etC = TMath::Sqrt(dCdq*dCdq*eq*eq + dCdE*dCdE*eBeamE*eBeamE);


	// Rad to deg
	tC  *= TMath::RadToDeg();
	etC *= TMath::RadToDeg();
}


//////////////////////////////////////////////////
//   Convert momentum transfer to Lab angle      //
//////////////////////////////////////////////////
void Converter::QToThetaLab(double &tL, double &etL, double beamE, double q, double eBeamE = 0.0, double eq = 0.0)
{
	// Momentum of primary in CM
	double p = m_MTar
	         * TMath::Sqrt( (2.0*m_MPri*beamE+beamE*beamE)
			              / ((m_MPri+m_MTar)*(m_MPri+m_MTar)+2.0*m_MTar*beamE) );
	double tC = 2.0 * TMath::ASin(0.5 * q / p);

	double A = TMath::Sqrt((m_MPri+m_MTar)*(m_MPri+m_MTar) + 2.0*m_MTar*beamE)
	         / (m_MPri + m_MTar + beamE);
	double B = (m_MPri*m_MPri + m_MPri*m_MTar + m_MTar*beamE)
	         / m_MTar / (m_MPri + m_MTar + beamE);

	double sint = TMath::Sin(tC);
	double cost = TMath::Cos(tC);

	double tan = A * sint / (cost + B);
	double aTan = TMath::ATan(tan);
	
	if ( aTan < 0.0 ) aTan += TMath::Pi();

	tL = aTan;


	double dCdq = 1.0 / TMath::Cos(0.5*tC);
	double dCdE = q * dCdq / p;
	double etC = TMath::Sqrt(dCdq*dCdq*eq*eq + dCdE*dCdE*eBeamE*eBeamE);

	double dAdE = - (m_MPri*m_MPri + m_MPri*m_MTar + m_MTar*beamE)
	            / (m_MPri + m_MTar + beamE) / (m_MPri + m_MTar + beamE)
		        / TMath::Sqrt((m_MPri+m_MTar)*(m_MPri+m_MTar) + 2.0*m_MTar*beamE);
	double dBdE = (m_MTar*m_MTar - m_MPri*m_MPri) / m_MTar
	            / (m_MPri + m_MTar + beamE) / (m_MPri + m_MTar + beamE);

	double dtdE = (dAdE*sint*(cost+B) - dBdE*A*sint)
	            / TMath::Cos(tan) / TMath::Cos(tan)
				/ (cost + B) / (cost + B);
	double dtdt = A * (cost*(cost+B) - sint*sint)
	            / TMath::Cos(tan) / TMath::Cos(tan)
				/ (cost + B) / (cost + B);

	etL = TMath::Sqrt(dtdE*dtdE*eBeamE*eBeamE + dtdt*dtdt*etC*etC);


	// Rad to deg
	tL  *= TMath::RadToDeg();
	etL *= TMath::RadToDeg();
}


//////////////////////////////////////////////////
//   Convert Lab angle to momentum transfer     //
//////////////////////////////////////////////////
void Converter::ThetaLabToThetaCM(double &tC, double &etC, double beamE, double tL, double eBeamE = 0.0, double etL = 0.0)
{
	// Deg to rad
	tL  *= TMath::DegToRad();
	etL *= TMath::DegToRad();

	double sign = 1;
	if ( tL >= 0.5*TMath::Pi() ) sign = -1;

	double A = TMath::Sqrt((m_MPri+m_MTar)*(m_MPri+m_MTar) + 2.0*m_MTar*beamE)
	         / (m_MPri + m_MTar + beamE);
	double B = (m_MPri*m_MPri + m_MPri*m_MTar + m_MTar*beamE)
	         / m_MTar / (m_MPri + m_MTar + beamE);
	double T = TMath::Tan(tL);

	double cost = (sign*A*TMath::Sqrt(T*T+A*A-B*B*T*T) - B*T*T)
	            / (T*T + A*A);
	
	tC = TMath::ACos(cost);
	etC = 0.0; // To be updated

	tC  *= TMath::RadToDeg();
	etC *= TMath::RadToDeg();
}


//////////////////////////////////////////////////
//   Convert Lab angle to momentum transfer     //
//////////////////////////////////////////////////
void Converter::ThetaLabToQ(double &q, double &eq, double beamE, double tL, double eBeamE = 0.0, double etL = 0.0)
{
	// Deg to rad
	tL  *= TMath::DegToRad();
	etL *= TMath::DegToRad();

	double sign = 1;
	if ( tL >= 0.5*TMath::Pi() ) sign = -1;

	double A = TMath::Sqrt((m_MPri+m_MTar)*(m_MPri+m_MTar) + 2.0*m_MTar*beamE)
	         / (m_MPri + m_MTar + beamE);
	double B = (m_MPri*m_MPri + m_MPri*m_MTar + m_MTar*beamE)
	         / m_MTar / (m_MPri + m_MTar + beamE);
	double T = TMath::Tan(tL);

	double cost = (sign*A*TMath::Sqrt(T*T+A*A-B*B*T*T) - B*T*T)
	            / (T*T + A*A);
	
	double tC = TMath::ACos(cost);
	double etC = 0.0; // To be updated

	// Momentum of primary in CM
	double p = m_MTar
	         * TMath::Sqrt( (2.0*m_MPri*beamE+beamE*beamE)
			              / ((m_MPri+m_MTar)*(m_MPri+m_MTar)+2.0*m_MTar*beamE) );

	q = 2.0 * p * TMath::Sin(0.5*tC);

	double dqdE = 2.0 * TMath::Sin(0.5*tC)
	            * m_MTar * m_MTar / p
				* ((m_MPri+beamE)*(m_MPri+m_MTar)*(m_MPri+m_MTar)+2*m_MTar*beamE*beamE)
				/ ((m_MPri)*(m_MTar)+2.0*m_MTar*beamE)
				/ ((m_MPri)*(m_MTar)+2.0*m_MTar*beamE);
	double dqdT = p * TMath::Cos(0.5*tC);

	eq = TMath::Sqrt(dqdE*dqdE*eBeamE*eBeamE + dqdT*dqdT*etC*etC);
}

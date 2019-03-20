#ifndef CONVERTER_H
#define CONVERTER_H

////////////////////////////////////////////////////////////////////////////////
//   Header file of Converter class.                                          //
//                                                                            //
//   This class is for conversion between variables such as lab scattering    //
// angle, CM scattering angle, momentum transfer and so forth. Use only unit  //
// set of following.                                                          //
//                                                                            //
//   {deg, GeV}                                                               //
//                                                                            //
//   Note! This class requires ROOT (https://root.cern.ch/) package.          //
//                                                                            //
//                                   - Hoyong Jeong (hyjeong@hep.korea.ac.kr) //
////////////////////////////////////////////////////////////////////////////////
class Converter
{
  public:
	Converter();
	Converter(double pri, double tar);
	~Converter();

	void SetMPri(double);
	void SetMTar(double);

	void ThetaCMToQ(double &q, double &eq, double beamE, double thetaCM, double eBeamE = 0.0, double eThetaCM = 0.0);
	void ThetaCMToT(double &t, double &et, double beamE, double thetaCM, double eBeamE = 0.0, double eThetaCM = 0.0);
	void ThetaCMToThetaLab(double &tL, double &etL, double beamE, double thetaCM, double eBeam = 0.0, double eThetaCM = 0.0);
	void QToThetaCM(double &thetaCM, double &eThetaCM, double beamE, double q, double eBeamE = 0.0, double eq = 0.0);
	void QToThetaLab(double &thetaLab, double &eThetaLab, double beamE, double q, double eBeamE = 0.0, double eq = 0.0);
	void ThetaLabToThetaCM(double &thetaCM, double &eThetaCM, double beamE, double thetaLab, double eBeamE = 0.0, double eThetaLab = 0.0);
	void ThetaLabToQ(double &q, double &eq, double beamE, double thetaLab, double eBeamE = 0.0, double eThetaLab = 0.0);

  private:
	double m_MPri;
	double m_MTar;
};

#endif

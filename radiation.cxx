
#include <globals.hxx> // for mesh object
#include <output.hxx>
#include <boutexception.hxx>
#include <utils.hxx>

#include "radiation.hxx"

#include <fstream>
#include <string>
#include <sstream>

const Field3D RadiatedPower::power(const Field3D &Te, const Field3D &Ne, const Field3D &Ni) {
  Field3D result;
  result.allocate();
  
  for(const auto& i : result) {
    result[i] = power(Te[i], Ne[i], Ni[i]);
  }
  
  return result;
}

InterpRadiatedPower::InterpRadiatedPower(const std::string &filename) {
  std::ifstream file(filename.c_str());
  
  output.write("Loading data from file: %s\n", filename.c_str());

  if(!file.is_open())
    throw BoutException("InterpRadiatedPower: Couldn't open file %s\n", filename.c_str());

  std::string line;
  int linenr = 1;
  while( std::getline(file, line) ) {
    // Expecting either a comment, blank line, or two numbers
    // Remove comments, then whitespace from left and right
    std::string strippedline = trim( trimComments( line ) );
    
    if (strippedline.length() == 0)
      continue;
    
    std::stringstream ss(strippedline);
    BoutReal t, p;
    if( !(ss >> t) )
      throw BoutException("InterpRadiatedPower: file '%s' line %d: %s\n",
                          filename.c_str(), linenr, line.c_str());
    
    if( !(ss >> p) )
      throw BoutException("InterpRadiatedPower: file '%s' line %d: %s\n",
                          filename.c_str(), linenr, line.c_str());
    
    te_array.push_back(t);
    p_array.push_back(p);
  
    linenr++;
  }
  

  file.close();
}

BoutReal InterpRadiatedPower::power(BoutReal Te, BoutReal ne, BoutReal ni) {
  return 0.0;
}

////////////////////////////////////////////////////////////////
// 

BoutReal HydrogenRadiatedPower::power(BoutReal Te, BoutReal ne, BoutReal ni) {
  throw BoutException("HydrogenRadiatedPower::power not implemented");
}

// Collision rate coefficient <sigma*v> [m3/s]
BoutReal HydrogenRadiatedPower::ionisation(BoutReal T) {
  double fION;	//collision rate coefficient <sigma*v> [m3/s]
  double TT,X,S;

  TT=T;
  
  if (TT<1.0) TT=1.0;
  X=log10(TT);
  
  if (TT>=20.0) 
    S=-0.5151*X-2.563/X-5.231;
  else
    S=-3.054*X-15.72*exp(-X)+1.603*exp(-X*X);
  
  fION=pow(10.0,S-6.0);           
  
  return fION;
}

//<sigma*v> [m3/s]
BoutReal HydrogenRadiatedPower::recombination(BoutReal n, BoutReal T) {
  double fREC;	//<sigma*v> [m3/s]
  double TT,RDNE,RTE,DNE,E,RN,RT,RNJ,RTI,suma;
  int i,j,i1,j1;
  
  if(n < 1e3) // Log(n) used, so prevent NaNs
    return 0.0;
  
  double MATA[9][9]= {
    {-2.855728479302E+01, -7.664042607917E-01, -4.930424003280E-03, \
     -5.386830982777E-03, -1.626039237665E-04,  6.080907650243E-06, \
     2.101102051942E-05, -2.770717597683E-06,  1.038235939800E-07,}, 
    { 3.488563234375E-02, -3.583233366133E-03, -3.620245352252E-03, \
      -9.532840484460E-04,  1.888048628708E-04, -1.014890683861E-05, \
      2.245676563601E-05, -4.695982369246E-06,  2.523166611507E-07,}, 
    {-2.799644392058E-02, -7.452514292790E-03,  6.958711963182E-03, \
     4.631753807534E-04,  1.288577690147E-04, -1.145028889459E-04, \
     -2.245624273814E-06,  3.250878872873E-06, -2.145390398476E-07,},  
    { 1.209545317879E-02,  2.709299760454E-03, -2.139257298118E-03, \
      -5.371179699661E-04, -1.634580516353E-05,  5.942193980802E-05, \
      -2.944873763540E-06, -9.387290785993E-07,  7.381435237585E-08,},
    {-2.436630799820E-03, -7.745129766167E-04,  4.603883706734E-04, \
     1.543350502150E-04, -9.601036952725E-06, -1.211851723717E-05, \
     1.002105099354E-06,  1.392391630459E-07, -1.299713684966E-08,},
    { 2.837893719800E-04,  1.142444698207E-04, -5.991636837395E-05, \
      -2.257565836876E-05,  3.425262385387E-06,  1.118965496365E-06, \
      -1.291320799814E-07, -1.139093288575E-08,  1.265189576423E-09,},
    {-1.886511169084E-05, -9.382783518064E-06,  4.729262545726E-06, \
     1.730782954588E-06, -4.077019941998E-07, -4.275321573501E-08, \
     7.786155463269E-09,  5.178505597480E-10, -6.854203970018E-11,},
    { 6.752155602894E-07,  3.902800099653E-07, -1.993485395689E-07, \
      -6.618240780594E-08,  2.042041097083E-08,  3.708616111085E-10, \
      -2.441127783437E-10, -9.452402157390E-12,  1.836615031798E-12,},
    {-1.005893858779E-08, -6.387411585521E-09,  3.352589865190E-09, \
     1.013364275013E-09, -3.707977721109E-10,  7.068450112690E-12, \
     3.773208484020E-12, -4.672724022059E-14, -1.640492364811E-14,},
  };
	
  RDNE=n;
  RTE=T;
	
  DNE=RDNE;
  TT=RTE;
  E=DNE*1.0E-14;
	
  if (TT<1.0) TT=1.0;
	
  RN=log(E);
  RT=log(TT);
	
  suma=0.0;
  for (i=1;i<=9;i++)
    {
      i1=i-1;
      for (j=1;j<=9;j++)
        {
          j1=j-1;
          RNJ=pow(RN,j1);
          if ((RN==0.0) && (j1==0)) RNJ=1.0;
          RTI=pow(RT,i1);
          if ((RT==0.0) && (i1==0)) RTI=1.0;		
          suma=suma+MATA[j-1][i-1]*RNJ*RTI;
        }
    }

  fREC=exp(suma)*1.0E-6;

  return fREC;
}

// <sigma*v> [m3/s]
BoutReal HydrogenRadiatedPower::chargeExchange(BoutReal Te) {
  double fCX;		//<sigma*v> [m3/s]
  double TT,S;
  
  TT=Te;
  
  if (TT<1.0) TT=1.0;    
  S=-14.0+log10(TT)/3.0;
  
  fCX=pow(10.0,S); 
  
  return fCX;
}

// <sigma*v> [m3/s]
BoutReal HydrogenRadiatedPower::excitation(BoutReal Te) {
  double fEXC;	//<sigma*v> [m3/s]
  double TT,Y;
  
  TT=Te;
  
  if (TT<1.0) TT=1.0;
  Y=10.2/TT;
  
  fEXC=49.0E-14/(0.28+Y)*exp(-Y)*sqrt(Y*(1.0+Y));
  
  return fEXC;
}


/////////////////////////////////////////////////////////////////////////////


BoutReal UpdatedRadiatedPower::power(BoutReal Te, BoutReal ne, BoutReal ni) {
  throw BoutException("UpdatedRadiatedPower::power not implemented");
}

// Collision rate coefficient <sigma*v> [m3/s]
BoutReal UpdatedRadiatedPower::ionisation(BoutReal T) {
    double fION; // Rate coefficient
    double TT;

    if (T < 0.025) {
      T = 0.025; // 300K
    }
    
    TT = T;

    double ioncoeffs[9] = {-3.271397E1, 1.353656E1, -5.739329, 1.563155, \
			   -2.877056E-1, 3.482560e-2, -2.631976E-3, \
			   1.119544E-4, -2.039150E-6};
    
    double lograte = 0.0;
    for (int i=0;i<=8;i++)
      {
	lograte = lograte + ioncoeffs[i]*pow(log(TT),i);
      }

    fION = exp(lograte)*1.0E-6;

    return fION;
}

//<sigma*v> [m3/s]
BoutReal UpdatedRadiatedPower::recombination(BoutReal n, BoutReal T) {
  double TT, RDNE, RTE, DNE, E, RN, RT, RNJ, RTI, suma, fHAV, fRAD;
  int i, j, i1, j1;

  if (n < 1e3) // Log(n) used, so prevent NaNs
    return 0.0;

  if (T < 0.025) {
    T = 0.025; // 300K
  }

  double MATA[9][9] = {
      {
          -2.855728479302E+01, -7.664042607917E-01, -4.930424003280E-03,
          -5.386830982777E-03, -1.626039237665E-04, 6.080907650243E-06,
          2.101102051942E-05, -2.770717597683E-06, 1.038235939800E-07,
      },
      {
          3.488563234375E-02, -3.583233366133E-03, -3.620245352252E-03,
          -9.532840484460E-04, 1.888048628708E-04, -1.014890683861E-05,
          2.245676563601E-05, -4.695982369246E-06, 2.523166611507E-07,
      },
      {
          -2.799644392058E-02, -7.452514292790E-03, 6.958711963182E-03,
          4.631753807534E-04, 1.288577690147E-04, -1.145028889459E-04,
          -2.245624273814E-06, 3.250878872873E-06, -2.145390398476E-07,
      },
      {
          1.209545317879E-02, 2.709299760454E-03, -2.139257298118E-03,
          -5.371179699661E-04, -1.634580516353E-05, 5.942193980802E-05,
          -2.944873763540E-06, -9.387290785993E-07, 7.381435237585E-08,
      },
      {
          -2.436630799820E-03, -7.745129766167E-04, 4.603883706734E-04,
          1.543350502150E-04, -9.601036952725E-06, -1.211851723717E-05,
          1.002105099354E-06, 1.392391630459E-07, -1.299713684966E-08,
      },
      {
          2.837893719800E-04, 1.142444698207E-04, -5.991636837395E-05,
          -2.257565836876E-05, 3.425262385387E-06, 1.118965496365E-06,
          -1.291320799814E-07, -1.139093288575E-08, 1.265189576423E-09,
      },
      {
          -1.886511169084E-05, -9.382783518064E-06, 4.729262545726E-06,
          1.730782954588E-06, -4.077019941998E-07, -4.275321573501E-08,
          7.786155463269E-09, 5.178505597480E-10, -6.854203970018E-11,
      },
      {
          6.752155602894E-07, 3.902800099653E-07, -1.993485395689E-07,
          -6.618240780594E-08, 2.042041097083E-08, 3.708616111085E-10,
          -2.441127783437E-10, -9.452402157390E-12, 1.836615031798E-12,
      },
      {
          -1.005893858779E-08, -6.387411585521E-09, 3.352589865190E-09,
          1.013364275013E-09, -3.707977721109E-10, 7.068450112690E-12,
          3.773208484020E-12, -4.672724022059E-14, -1.640492364811E-14,
      },
  };

  RDNE = n;
  RTE = T;

  DNE = RDNE;
  TT = RTE;
  E = DNE * 1.0E-14;

  RN = log(E);
  RT = log(TT);

  suma = 0.0;
  for (i = 1; i <= 9; i++) {
    i1 = i - 1;
    for (j = 1; j <= 9; j++) {
      j1 = j - 1;
      RNJ = pow(RN, j1);
      if ((RN == 0.0) && (j1 == 0))
        RNJ = 1.0;
      RTI = pow(RT, i1);
      if ((RT == 0.0) && (i1 == 0))
        RTI = 1.0;
      suma = suma + MATA[j - 1][i - 1] * RNJ * RTI;
    }
  }

  fHAV = exp(suma) * 1.0E-6 / (1.0 + 0.125 * TT);

  double A = 3.92E-20;
  double B = 3.0E-124 * pow(1.6E-19, -4.5);
  double Ry = 13.60569;
  double chi = 0.35;

  fRAD = A * pow(Ry, 1.5) * 1.0 / (sqrt(TT) * (Ry + chi * TT));

  return fHAV + fRAD + B * DNE * pow(TT, -5.0);
}

// <sigma*v> [m3/s]
BoutReal UpdatedRadiatedPower::chargeExchange(BoutReal T) {
  if (T < 0.025) {
    T = 0.025; // 300K
  }

  double E = 10.0;

  double cxcoeffs[9][9] = {
      {
          -1.829079582E1, 1.640252721E-1, 3.364564509E-2, 9.530225559E-3,
          -8.519413900E-4, -1.247583861E-3, 3.014307546E-4, -2.499323170E-5,
          6.932627238E-7,
      },
      {
          2.169137616E-1, -1.106722014E-1, -1.382158680E-3, 7.348786287E-3,
          -6.343059502E-4, -1.919569450E-4, 4.075019352E-5, -2.850044983E-6,
          6.966822400E-8,
      },
      {
          4.307131244E-2, 8.948693625E-3, -1.209480567E-2, -3.675019470E-4,
          1.039643391E-3, -1.553840718E-4, 2.670827249E-6, 7.695300598E-7,
          -3.783302282E-8,
      },
      {
          -5.754895093E-4, 6.062141761E-3, 1.075907882E-3, -8.119301728E-4,
          8.911036876E-6, 3.175388950E-5, -4.515123642E-6, 2.187439284E-7,
          -2.911233952E-9,
      },
      {
          -1.552077120E-3, -1.210431588E-3, 8.297212634E-4, 1.361661817E-4,
          -1.008928628E-4, 1.080693990E-5, 5.106059414E-7, -1.299275586E-7,
          5.117133050E-9,
      },
      {
          -1.876800283E-4, -4.052878752E-5, -1.907025663E-4, 1.141663042E-5,
          1.775681984E-5, -3.149286924E-6, 3.105491555E-8, 2.274394089E-8,
          -1.130988251E-9,
      },
      {
          1.125490271E-4, 2.875900436E-5, 1.338839629E-5, -4.340802793E-6,
          -7.003521917E-7, 2.318308730E-7, -6.030983538E-9, -1.755944926E-9,
          1.005189187E-10,
      },
      {
          -1.238982763E-5, -2.616998140E-6, -1.171762874E-7, 3.517971869E-7,
          -4.928692833E-8, 1.756388999E-10, -1.446756796E-10, 7.143183138E-11,
          -3.989884106E-12,
      },
      {
          4.163596197E-7, 7.558092849E-8, -1.328404104E-8, -9.170850254E-9,
          3.208853884E-9, -3.952740759E-10, 2.739558476E-11, -1.693040209E-12,
          6.388219930E-14,
      },
  };

  double lograte = 0.0;
  for (int i = 0; i <= 8; i++) {
    for (int j = 0; j <= 8; j++) {
      lograte = lograte + cxcoeffs[i][j] * pow(log(T), i) * pow(log(E), j);
    }
  }

  return 1.0E-6 * exp(lograte);
}

// <sigma*v> [m3/s]
BoutReal UpdatedRadiatedPower::excitation(BoutReal Te) {
  double fEXC;	//<sigma*v> [m3/s]
  double TT,Y;
  
  TT=Te;
  
  if (TT<1.0) TT=1.0;
  Y=10.2/TT;
  
  fEXC=49.0E-14/(0.28+Y)*exp(-Y)*sqrt(Y*(1.0+Y));
  
  return fEXC;
}

///////// Updated by Yulin//////////

                      ///////// New reactions   /////////
//  1.Dissicaiation ionization rate coefficient (DI); T_range = [3eV~30eV].
//    e- + H2 -> e- + e- + H2+ -> H+ + H(1s) + 2e-  //

//  2.Moelecules dissociated by collision with atoms (Diss_H_H2)
//    H(1s) + H2(v) -> H(1s) + 2H(1s)

//  3.Molecures dissciated by collision with molecures(Diss_H2_H2), T_range=[2eV~30eV].
//    H2(v)+H2(v) -> H2(v)+2H(1s)

//  4.Dissociation Attachment(DA_1 and DA_2),Channel 1 and 2
//    channel 1: e- + H2(v) -> H2(-1) - > e- + H(1s) + H(1s)
//    channel 2: e- + H2(v) -> H2(-1) - > e- + H(1s) + H(1s)

//  5.Dissociation
//    e- + H2(v) -> H2(-1) - > H- + H(1s)
//  6.recombination bt collision between H- and H
//    H- + H -> H2 + e

//  7.Molecule ionized by collision with electrons
//    e + H2 -> 2e + H2+
//  8.H2+ ionized by collision with electrons
//    e + H2+ -> 2e + 2H+
  
// <sigma*v> [m3/s] for Dissicaiation ionization (DI)
BoutReal UpdatedRadiatedPower::dissociation_ion(BoutReal T) {
  if (T < 3.0) {
    T = 3.0; // 300K
  }
  double c1,c2,c3,c4,c5,DI_ln_rate,DI_RateCoefficient,TT;
  TT=T;
  double DIcoeffs[14][5] = {
{-2.1196e2,1.0022,-20.350,-4.5201,1.0773e-2,},
{-2.0518e2,0.99226,-19.905,-3.3364,1.1725e-2,},
{-1.9936e2,0.98837,-19.600,-3.0891,1.2838e-2,},
{-1.9398e2,0.98421,-19.457,-3.1386,1.3756e-2,},
{-1.8893e2,0.97647,-19.397,-3.2807,1.4833e-2,},
{-1.8422e2,0.96189,-19.310,-3.2609,1.6030e-2,},
{-1.7903e2,0.94593,-19.170,-3.0592,1.7254e-2,},
{-1.7364e2,0.93986,-19.052,-2.9880,1.8505e-2,},
{-1.6960e2,0.93507,-18.908,-2.7334,1.8810e-2,},
{-1.6664e2,0.92602,-18.723,-2.2024,1.8055e-2,},
{-1.6521e2,0.92124,-18.549,-1.6895,1.6245e-2,},
{-1.6569e2,0.93366,-18.479,-1.6311,1.5194e-2,},
{-1.6464e2,0.94682,-18.440,-1.7259,1.5304e-2,},
{-1.6071e2,0.95533,-18.405,-1.8938,1.6254e-2,}};

double delt_Ev[14]={30.6,27.4,25.2,23.4,21.9,20.6,19.6,18.5,17.6,16.7,16.0,15.4,14.8,14.2,};

double T_Kelvin=TT*0.001*11600.0; //in Kelvin * 0.001

DI_RateCoefficient = 0.0;

for (int i = 0;i <= 13;i++) {
     c1=DIcoeffs[i][0];
     c2=DIcoeffs[i][1];
     c3=DIcoeffs[i][2];
     c4=DIcoeffs[i][3];
     c5=DIcoeffs[i][4];
     DI_ln_rate= c1/pow(T_Kelvin,c2)+c3+c4*exp(-c5*T_Kelvin);
     DI_RateCoefficient += exp(DI_ln_rate-delt_Ev[i]/TT)*1.0e-6;
}
return DI_RateCoefficient;
}

////  2.Moelecules dissociated by collision with atoms (Diss_H_H2)
//    H(1s) + H2(v) -> H(1s) + 2H(1s)
BoutReal UpdatedRadiatedPower::Diss_H_H2(BoutReal T) {
  if (T < 3.0) {
    T = 3.0; // 300K
  }
  double a1,a2,a3,a4,a5,TT,pow_a3,pow_a5,Diss_H_H2_ln_rate,Diss_H_H2_RateCoefficient;
  TT = T;
  double Diss_H_H2_coeffs[15][5] = {
{2.06964e1, 7.32149e7, 1.74660, 4.75874e3, -9.42775e-1,},
{2.05788e1, 4.32679e7, 1.68520, 1.91812e3, -8.16838e-1,},
{2.05183e1, 5.15169e7, 1.73345, 3.09006e3, -8.88414e-1,},
{2.04460e1, 1.87116e8, 1.87951, 9.04442e3, -9.78327e-1,},
{2.03608e1, 4.93688e8, 1.99499, 2.32656e4, -1.06294,},
{2.02426e1, 1.80194e8, 1.92249, 1.28777e4, -1.02713,},
{2.00161e1, 2.96945e5, 1.31044, 9.55214e2, -1.07546,},
{1.98954e1, 4.53104e5, 1.37055, 3.88065e2, -8.71521e-1,},
{1.97543e1, 5.13174e5, 1.39819, 3.54272e2, -8.07563e-1,},
{1.97464e1, 9.47230e4, 1.24048, 2.28283e2, -8.51591e-1,},
{1.95900e1, 6.43990e4, 1.22211, 1.16196e2, -7.35645e-1,},
{1.94937e1, 3.49017e4, 1.20883, 1.26329e2, -8.15130e-1,},
{1.90708e1, 1.05971e5, 9.91646e-1, 1.05518e2, -1.93837e-1,},
{1.89718e1, 7.76046e5, 7.84577e-1, 1.31409e3, -1.00479e-2,},
{1.87530e1, 5.81508e5, 7.35904e-1, 1.69328e3, 4.47757e-3,}};

double Diss_H_H2_delt_Ev[15] = {30.6,27.4,25.2,23.4,21.9,20.6,19.6,18.5,17.6,16.7,16.0,15.4,14.8,14.2,13.7,}; // in eV

double T1=TT*11600; //in Kelvin

Diss_H_H2_RateCoefficient = 0.0;
for (int i = 0;i <= 14;i++) {
     a1=Diss_H_H2_coeffs[i][0];
     a2=Diss_H_H2_coeffs[i][1];
     a3=Diss_H_H2_coeffs[i][2];
     a4=Diss_H_H2_coeffs[i][3];
     a5=Diss_H_H2_coeffs[i][4];
     pow_a3=pow(T1,a3);
     pow_a5=a4*pow(T1,a5)+1;
     Diss_H_H2_ln_rate = -a1-a2/(pow_a3*pow_a5);
     Diss_H_H2_RateCoefficient += exp(Diss_H_H2_ln_rate-Diss_H_H2_delt_Ev[i]/TT)*1.0e-6;
}

return Diss_H_H2_RateCoefficient;
}


//  7.Molecule ionized by collision with electrons
//    e + H2 -> 2e + H2+
BoutReal UpdatedRadiatedPower::Ion_H2(BoutReal T) {
  if (T < 3.2) {
    T = 3.2; 
  }
  double Ion_H2_ln_rate,Ion_H2_RateCoefficient,rate_c,TT;
  TT = T;
  double Ion_H2_coeffs[9]= {-3.56864029E1, 1.73346898E1, -7.76746936, 2.211579405,-4.169840174E-1, 5.088289820e-2, -4.832737518E-3,1.612863120E-4, -2.8933919044E-6};
  Ion_H2_RateCoefficient=0.0;
  for (int i=0;i<=8;i++){
   Ion_H2_RateCoefficient += Ion_H2_coeffs[i]*pow(log(TT),i);
}
rate_c = exp(Ion_H2_RateCoefficient)*1.0e-6;
return rate_c;
}


//  8.H2+ ionized by collision with electrons
//    e + H2+ -> 2e + 2H+

BoutReal UpdatedRadiatedPower::Ion_charged_H2(BoutReal T) {
  if (T < 3.2) {
    T = 3.2; 
  }
  double Ion_charged_H2_ln_rate,Ion_charged_H2_RateCoefficient,rate_i,TT;
  TT = T;
  double Ion_charged_H2_coeffs[9]= {-3.7461923E1, 1.55935503E1, -6.6932383, 1.98170029,-4.04482088E-1, 5.35239162e-2, -4.31745184E-3,1.91849987E-4, -3.5917797E-6};
  Ion_charged_H2_RateCoefficient=0.0;
  for (int i=0;i<=8;i++){
   Ion_charged_H2_RateCoefficient += Ion_charged_H2_coeffs[i]*pow(log(TT),i);
}
rate_i = exp(Ion_charged_H2_RateCoefficient)*1.0e-6; 
return rate_i;
}







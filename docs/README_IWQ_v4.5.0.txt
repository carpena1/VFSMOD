3 -11.5142 0.5949 0.4892 -0.3753 0.2039 ;IWQPRO CSAB(I)
1 1000 1.                               ; IKD (Kd or Koc) (%OC)
40.                                     ; %Clay content in field soil
1                                       ; IDG
0 30 0.28 100. 2. 0.05 0.582368E+00     ; ndgday dgHalf FC dgPin dgML dgLD dgmres0
11.1 13.4 14.9 14.1 16.9 18.9 20.4 17.4 16.1 14.8 14.1           ; dgT(i)
0.274 0.272 0.272 0.276 0.267 0.257 0.247 0.238 0.23 0.226 0.223 ; dgTheta(i)
1                                       ; IMOB

------------------------------------------------------------------
IWQPRO  : Pesticide trapping: 1=Sabbagh;2= Sabbagh(refit);3=mech.mass bal.;4=Chen
CSAB(I) : Coefficients for refitted Sabbagh equation (used when IWQPRO=2)
IKD     : Sorption type: 0, Kd(L/Kg); 1, Koc (L/Kg)
Kd Koc  : Sorption coefficient, distribution Kd (IKD= 0) or Koc (IKD=1) (L/Kg)
%OC     : Source soil organic carbon, only read when IKD=1 (Koc) (%)
IDG     : Degradation type: 1: EU-FOCUS k=Kref.k(T).k(theta); 2: US-EPA k=Kref;
         3: k=Kref.k(T);4: k=Kref.k(1theta; 0: No degradation
ndgday  : no. of days between events (d)
dgHalf  : t0.5, pesticide half-life (d), kref=Ln2/t0.5
FC      : top soil field capacity (m3/m3). This can be taken from FOCUS R1-R4 scenario
         parameters used by the PRZM model
dgPin   : Pin, pesticide mass entering filter for event over source area (mg/m2)
dgML    : Surface mixing layer thickness (cm, standard= 2cm PRZM)
dgLD    : lambda, dispersion length of chemical (m). This can be taken as  
          0.05m from FOCUS-Pearl (Default)
dgmres0 : Pesticide residue on VFS surface (mixing layer) when event starts (mg/m2)
dgT(i)  : T, daily air temperatures (C) for period between events, PRZM weather
dgTheta(i) : theta, Topsoil daily volumetric moisture (-) for period between events
IMOB    : Residues remobilization: 1(or none): partial (recomm); 2:full; 3:no remob.

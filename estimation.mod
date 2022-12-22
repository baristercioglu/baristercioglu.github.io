var
    TO TS Y GDP PiY pE pER pEF pZ pW LH LL pLH pLL pL pO_RoW
    W Zd Ze Ze_RoW CH CL CG IPZ IPR IPF IGZ IGR IGF Y_RoW %9
    A BP BG BW KZ KR KF KPZ KPR KPF KGZ KGR KGF %23
    PiH PiL lamH lamL Lam Lam_RoW QZ QR QF pSZ pSR pSF vph pZ_RoW RER dNER PiY_RoW %40
    R RB RS RW %44
    T TCH TCL TA TR TH TL TSZ TSR TSF TPZ TPR TPF TEF %52
    NL LZ LR LF ELe ELi %58
    g DZ BZ DR BR DF BF uZ EZ EY ER EF SZ SF O GamW 
    Vg VWL VWH VR VR_RoW VpO_RoW %77
    VC VN VGamB VGamS VPZ VL VZ VK VCG VIGZ VIGR VIGF VGamW VY_RoW VpZ_RoW VER VERT VEF VY VW VWI VW_RoW VEFT VPiY_RoW%88
%    psiC psiNH psiNL TLL lam tauZ1 gamW alpW_RoW alpL alpSZ alpW alpZ alpER alpEF alpE alpSF alpY alpKZ alpKR alpKF
    taxbase_TC taxbase_TA  taxbase_TR  taxbase_TH  taxbase_TL taxbase_TS taxbase_TP taxbase_TEF taxbase_TO 
    CAPEXshare_Z OPEXshare_Z DIVshare_Z 
    CAPEXshare_R OPEXshare_R DIVshare_R
    CAPEXshare_F OPEXshare_F DIVshare_F 
    Con IG IP Inv Ex Im GDP_m E TB EF_ann ER_ann EY_ann EZ_ann O_ann
    tOd tEFd tSd tLd tHd tCd tPd tRd tAd tOd_add tEFd_add tSd_add tLd_add tHd_add tCd_add tPd_add tRd_add tAd_add
    dlnGDP dlnC dlnIG dlnIP dlnCG dlnEx dlnIm dlnEL dlnpL dlnE dlnpE dlnpZ dlnpEF dlnY_RoW dlnpZ_RoW dlnpO_RoW dlnEF dlnER dlnEY dlnEZ dlnO 
    dlnEF_ann dlnER_ann dlnEY_ann dlnEZ_ann dlnO_ann
    devGDP devC devIG devIP devCG devEx devIm devEL devpL devR devRB devRW devRER devPiY devPiY_RoW devE devu devpE devpZ devpEF devY_RoW devpZ_RoW devpO_RoW
    devEF_ann devER_ann devEY_ann devEZ_ann devO_ann
 % GDP and components (percent change from pre-policy level)
 lGDP lConP lConH lConL lIVP lIVG lExp lImp lImp_WO ppTB    

% GDP and components (normalized by GDP under pre-policy prices)
 ppU ConP_GDP ConH_GDP ConL_GDP IVP_GDP IVG_GDP IVG_Z_GDP IVG_R_GDP IVG_F_GDP Exp_GDP Imp_GDP  

% Sectoral output (percent change from pre-policy level)
 lZ lER lEF lEY lEZ lE lW 

 % Sectoral output (normalized by GDP under pre-policy prices)
 Z_GDP ER_GDP EF_GDP EY_GDP EZ_GDP E_GDP W_GDP lLZ lLR lLF    

% prices
 lpZ lpER lpEF lpE lpW lpZ_RoW lpLH lpLL lpL     

% private Investment (percent change from pre-policy level)
 lIPZ lIPR lIPF lKPZ lKPR lKPF lpSZ lpSR lpSF ppuZ 

% Tax revenues
ppT ppTO ppTEF ppTS ppTL ppTH ppTC ppTP ppTR ppTA 

% Prices, inflation, emissions, wealth
ppPiY ppR ppRP ppRER lO lptO ppAY ppBGY ppBPY ppBWY ppGamW  
;



varexo epsg epsWL epsWH epsR epsRW epsC epsN epsRB epsRS epsPZ epsZ epsK epsCG epsIG epsIGR epsY_RoW epspZ_RoW epsGamW epsER epsERT epsEF epsY epsW epsWI epsW_RoW epsEFT epsL epspO_RoW epsPiY_RoW 
       epsGDP_m
       tOd_shock tEFd_shock tSd_shock tLd_shock tHd_shock tCd_shock tPd_shock tRd_shock tAd_shock
;

parameters
    bet psiA GamB GamS phiW  
    kap vrh
    etaH etaL
    tauW tauP tauI tauZ2 tauO tauZe alpO psiu
    muP muW delZ delR delF sigSZ sigSF sigL sigW sigKZ sigKR sigKF sigZ sigE sigY sigW_RoW sigER sigEF
    tC tA tR tH tL tS tP tEF tO
    phiWE phiWP phiCY phiIY phiRP phiRY phidRP phidRY xiELe
    rhog rhoWL rhoWH rhoR rhoRW rhoW rhoWI rhoW_RoW rhopO_RoW
    rhoC rhoN rhoRB rhoRS rhoPZ rhoL rhoZ rhoK rhoCG rhoIG rhoIGR rhoY_RoW rhopZ_RoW rhoER rhoERT rhoEF rhoY rhoEFT rhoPiY_RoW 
    psiC psiNH psiNL TLL lam tauZ1 gamW alpW_RoW alpL alpSZ alpW alpZ alpER alpEF alpE alpSF alpY alpKZ alpKR alpKF
    g_bar PiY_bar PiH_bar PiL_bar R_bar CG_bar IGZ_bar IGR_bar IGF_bar Y_RoW_bar pZ_RoW_bar RERBWY_bar RW_bar pO_RoW_bar PiY_RoW_bar
    GDP_bar ELe_bar BGY_bar BPY_bar uZ_bar LHLL_bar RER_bar CHCL_bar EYEZ_bar ZeRoWY_bar EREF_bar OY_bar EY_bar WYZ_bar WYE_bar LRLF_bar IPZY_bar IPRY_bar IPFY_bar
    RBbar Cbar IGbar IPbar Exbar Imbar RERbar pLbar pEbar pZbar Ebar EFbar ERbar EYbar EZbar Obar pEFbar
    ;



% Household sector:
bet      = 0.99;
psiA     = 1;	  %given endogenous psiC and psiNH, this parameter has not effect. 
GamB     = 0.001; %10 year bond rate over 2010-2016 is 0.38% above 3-month interest rate
GamS     = 0.013; %Return on shares over 2010-2016 is 5.25% above 10 year bond rate
phiW     = 1;	  %should be estimated
kap      = 0.6517;
vrh      = 0.6883;
etaH     = 1.5140;

% Firm sector:
muP      = 0.2;	  %needs review
muW      = 0.2;	  %needs review
delZ     = 0.025; %Argentiero (2017)
delR     = 0.03;  %Argentiero (2017)
delF     = 0.02;  %Argentiero (2017)

%Input shares
alpO     = 0.87;  %87% of carbon mining is imported

%Elasticities
sigY     = 0.3;  %elasticity of substitution between core consumption and energy; Medina Soto 2005 estimate 0.65 as elast betw oil and core cons. An and Kang 2011 estimate 0.3 for south korea btw oil and core cons
sigSZ    = 0.5;  %van der Werf (2008) finds mean elasticity 0.4 and reports values of 0.4-0.5 in the literature
sigL     = 2;  %david autor 2018 lectures notes and Angrist’s (1995) state sigma=2;
sigW     = 2;    %Bajzik et al. 2020 estimate trade elasticity between 2.5 and 5.1 with a median of 3.8 % LARGE ELASTICITY DRIVES NEGATIVE CARBON TAX MULTIPLIER
sigKZ    = 3;    %should be around 3 for EME. https://www.imf.org/-/media/Files/Publications/WP/2019/wpiea2019232-print-pdf.ashx estimated for Adv countries by Zidong An, Alvar Kangur, and Chris Papageorgiou IMF 2019; they also confirm Cobb Douglas for capital-labor as parameter is insignificant! Zidong An, Alvar Kangur, and Chris Papageorgiou IMF 2019 find parameters that are all over the place and insignificant. With the substitution parameter (not elasticity!) insignificant, elasticity = 1 and we are Cobb Douglas!
sigKR    = 3;
sigKF    = 3;
sigZ     = 0.9;	 %NAWM, GM, E-QUEST and SW assume Cobb-Douglas but Gechert et al. (2022) report lower values. 
sigSF    = 0.4;  %0.5393;
sigE     = 6;    %E-QUEST
sigW_RoW = 2;
sigER    = 0.3;  %should be estimated
sigEF    = 0.3;  %should be estimated

% Adjustment cost scaling paramters:
tauW     = 55.2007;
tauP     = 61.7619;
tauI     = 3.3819;
tauZ2    = 0.9665;
tauO     = 1.7848;
tauZe    = 4.1079;

%Scaling parameters
psiu     = 1;

% Adjustment elasticities:
phiWE    = 0.4645;
phiWP    = 0.6809;
xiELe    = 0.8659;
phiCY    = -0.1656;
phiIY    = -0.0316;
phiRP    = 1.7332; %Try with lower phiRP and not lagged!
phiRY    = 0.1720;
phidRP   = 0.3120;
phidRY   = 0.1375;

% Taxes, transfers, and tax rates: OECD Revenue Statistics 2017: Tax/GDP is 25% https://data.oecd.org/tax/tax-revenue.htm#indicator-chart
tC       = 0.16;    %Taxes on goods and services are 11% of GDP
tA       = 0.0003;  %Property taxes are 1% of GDP
tR       = 0.008;   %PIT is 4% of GDP
tH       = 0.08;    %PIT is 4% of GDP
tL       = 0.06;    %PIT is 4% of GDP
tS       = 0.17;    %Social Security Payroll taxes are 7% of GDP
tP       = 0.045;   %Corporate profit taxes are 2% of GDP
tEF      = 0.025;   %Energy taxes are 1.412% of GDP https://data.oecd.org/envpolicy/environmental-tax.htm
tO       = 0.00001; 

% Persistence parameters of endogenous and exogenous shock processes:
rhog      = 0.7254;
rhoWL     = 0.3906;
rhoWH     = 0.5;	%why not estimated? Cannot find mode!
rhoR      = 0.6259;	%should be higher for simulation
rhoRW     = 0.5;	%why not estimated? Cannot find mode!
rhoC      = 0.4464;
rhoN      = 0.5;	%why not estimated? Cannot find mode!
rhoRB     = 0.5;	%why not estimated? Cannot find mode!
rhoRS     = 0.5;	%why not estimated? Cannot find mode!
rhoPZ     = 0;
rhoL      = 0.7825;
rhoZ      = 0.6124;
rhoK      = 0.5303;
rhoCG     = 0.6256;
rhoIG     = 0.6866;
rhoIGR    = 0.6866;	%assumed to be equal to rhoIG
rhoY_RoW  = 0.6862;
rhopZ_RoW = 0.6466;
rhoER     = 0.5745;
rhoERT    = 0.8518;
rhoEF     = 0.5575;
rhoEFT    = 0.7401;
rhoY      = 0.8768;
rhoW      = 0.7216;	
rhoWI     = 0.8241;
rhoW_RoW  = 0.3592;
rhopO_RoW = 0.7262;
rhoPiY_RoW= 0.5;	%why not estimated? Cannot find mode!


% Exogenous steady-state values (no parameter needs to be restricted to realize the steady-state value)
g_bar     = 1.0083;   %3.32% annual GDP/LF growth over 2010-2016 (FRED)
PiY_bar   = 1.0199;   %7.96% mean inflation rate over 2010-2016 (FRED)
PiH_bar   = g_bar*PiY_bar;
PiL_bar   = g_bar*PiY_bar;
R_bar     = 1.0219;   %8.76% short-term interest rate over 2010-2016(FRED)
RW_bar    = 1.0036;   %1.45% mean 10-year government bond rate for Germany over 2010-2016
CG_bar    = 0.14;     %World Bank
IGZ_bar   = 0.036;    %for 2017 government gross investment/GDP is 0.04 (Turkish Statistical Institute + OECD)
IGR_bar   = 0.001;    %public energy investment is 0.004 https://www.sbb.gov.tr/yatirimlar/yatirimlarin-sektorel-dagilimi/
IGF_bar   = 0.003;    %expert opinion is asked for the split
RERBWY_bar= -1.6;     %https://en.wikipedia.org/wiki/Net_international_investment_position
Y_RoW_bar = 1;
pZ_RoW_bar= 1;
PiY_RoW_bar = 1.0031; %1.2% mean inflation rate for the Euro Zone over 2010-2016
pO_RoW_bar= 0.07;    %I chose pO such that total carbon input into energy sector is 3% of GDP

% Endogenous steady-state values (some parameter needs to be restricted to realize the steady-state value)
GDP_bar   = 1;        %psiC
ZeRoWY_bar= 0.26;     %alpW 0.26, averaged over 2000-2019, import share %with Ybar=1 and RER=1, 1-alpY is the steady-state imports/GDP ratio. ZeRoWY is the ratio of the value of imports in domestic currency and deflated by local GDP deflator and domestic real GDP.
LHLL_bar  = 1;        %psiNH US data (No Turkish data available for wealth distribution)
ELe_bar   = 0.9;      %psiNL unemployment rate is 10% (2006-2020) https://fred.stlouisfed.org/series/LRUN64TTTRQ156S
BGY_bar   = 1.32;     %TLL Debt securities: 20% domestic 12% foreign for 2020 https://stats.bis.org/statx/srs/table/c3?c=TR&f=pdf The rest is loans - Page 8 - https://www.tcmb.gov.tr/wps/wcm/connect/960e0b6e-2fd0-43e9-85ef-4fdd363d425e/2017-I.Quarter.pdf?MOD=AJPERES&CACHEID=ROOTWORKSPACE-960e0b6e-2fd0-43e9-85ef-4fdd363d425e-m5lU.P0
BPY_bar   = 2.3;      %lam 
IPZY_bar  = 0.223;    %alpKZ total private investment/GDP averaged over 2002-2019 is 23%
IPRY_bar  = 0.003;   %alpKR private energy investment averaged over 2002-2019 is 0.007 of GDP (currently 0.006 because of alpKER and alpKEF)
IPFY_bar  = 0.004;   %alpKF we need data about the split
WYZ_bar   = 0.39;     %alpZ  Gloria - 2014 - the ratio of labor income in core goods to core goods (value-added+energy bill)
WYE_bar   = 0.08;     %alpEF Gloria - 2014 - the ratio of labor income in energy to energy sectors' (value-added+carbon bill)
LRLF_bar  = 0.15;     %alpER - GDLD and GTAP-Power - the ratio of renewable employment to fossil employment
EY_bar    = 0.067;    %alpY GLORIA - 2014 - the share of (value-added of carbon mining+energy sector+total carbon bill of the economy) in GDP 
EYEZ_bar  = 0.94;     %alpSZ, average over 2000-2019, IEA data, the share of household energy consumption (residential+other final consumption+half of transportation) to firm energy consumption (industry+commercial/public+half of transportation) 
EREF_bar  = 0.14;     %alpE The ratio between renewable energy and fossil energy, average of 2000-2019.Page 18 - (https://iea.blob.core.windows.net/assets/cc499a7b-b72a-466c-88de-d792a9daff44/Turkey_2021_Energy_Policy_Review.pdf) 
CHCL_bar  = 2;        %alpL TUIK data on distribution of consumption by income quintiles
OY_bar    = 0.46;     %alpSF CO2 emission, average of 2000-2019 https://data.worldbank.org/indicator/EN.ATM.CO2E.KD.GD?locations=TR 
uZ_bar    = 1;        %tauZ1
RER_bar   = 1;        %alpW_RoW

% All of these parameters are overwritten in the steadystate file to satisfy the endogenous steady-state calibration targets above

                   
psiC     = 0;
psiNH    = 0;
psiNL    = 0;
TLL      = 0;
lam      = 0;
gamW     = 0;
alpW_RoW = 0;
alpSZ    = 0;
alpZ     = 0;
alpW     = 0;
alpS     = 0;
alpL     = 0;
alpE     = 0;
alpER    = 0;
alpEF    = 0;
alpY     = 0;
alpSF    = 0;
tauZ1    = 0;
alpKZ    = 0;
alpKR    = 0;
alpKF    = 0;

RBbar    = 0;
Cbar     = 0;
IGbar    = 0;
IPbar    = 0;
Exbar    = 0;
Imbar    = 0;
RERbar   = 0;
pLbar    = 0;
pEbar    = 0;
pZbar    = 0;
Ebar     = 0;
EFbar    = 0;
ERbar    = 0;
EYbar    = 0;
EZbar    = 0;
Obar     = 0;
pEFbar   = 0;

etaL     = 0;
sigKR    = 0;
sigKF    = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

model;

% High-skilled labor households (continuum of differentiated labor services plus perf. competitive homogenizer)
(1+tCd)*lamH = VC*psiC*(CH-kap*CH(-1)/g)^(-vrh);
(1+tAd)*lamH - psiA = bet*(RB-tRd*(RB-1))/PiY(+1)*lamH(+1)/g(+1);
(1+tAd)*lamH - psiA = bet*(1-GamW)*(RW-tRd*(RW-1))/PiY(+1)*dNER(+1)*lamH(+1)/g(+1);
dNER=RER/RER(-1)*PiY/PiY_RoW;
GamW = gamW/phiW*(VGamW^(1/gamW)*exp(RER*BW/GDP)^phiW-1);
(1+tAd)*lamH - psiA = bet*(1-GamS*VGamS^(1/GamS))*(RS-tRd*(RS-1))/PiY(+1)*lamH(+1)/g(+1);
(1-tHd)*pLH*LH = VN*psiNH*LH^(1+etaH)/lamH*(1+VWH*muW) - tauW*VWH*muW*((PiH-PiH_bar)*PiH - bet*lamH(+1)/lamH*(PiH(+1)-PiH_bar)*PiH(+1));
CH + A + TCH + TA + TH + TR = pLH*LH + RB(-1)/PiY*(BG(-1)+BP(-1))/g + RW(-1)/PiY_RoW*RER*BW(-1)/g + RS(-1)/PiY*(pSZ(-1)+pSR(-1)+pSF(-1))/g + (1-alpO)*RER*pO_RoW*O;
pSZ = Lam(+1)*g(+1)*(pSZ(+1)+DZ(+1));
pSR = Lam(+1)*g(+1)*(pSR(+1)+DR(+1));
pSF = Lam(+1)*g(+1)*(pSF(+1)+DF(+1));
Lam = (RS/PiY(+1))^(-1);

% Low-skilled labor households
(1+tCd)*lamL = VC*psiC*(CL-kap*CL(-1)/g)^(-vrh);
VN*psiNL*(ELi*NL)^etaL = (1-tLd)*lamL*pLL;
CL + TCL + TL = pLL*LL + TLL;

% Retailers (combining cars with fuel)
W = VY^(sigY-1)*alpY*pW^(-sigY)*Y;
EY = VY^(sigY-1)*(1-alpY)*pE^(-sigY)*Y;
1 = VY^(-1)*(alpY*pW^(1-sigY) + (1-alpY)*pE^(1-sigY))^(1/(1-sigY));

% Wholesalers (combining domestic and foreign unfueled cars)
Zd = VW^(sigW-1)*alpW*(pZ/pW)^(-sigW)*W;
Ze_RoW = VWI^(sigW-1)*(1-alpW)*(RER*pZ_RoW/pW)^(-sigW)*W;
pW = (alpW*VW^(sigW-1)*pZ^(1-sigW) + (1-alpW)*VWI^(sigW-1)*(RER*pZ_RoW)^(1-sigW))^(1/(1-sigW));

% Intermediate goods firms (continuum of differentiated good producers plus perf. competitive homogenizer)
g = g(-1)^rhog*g_bar^(1-rhog)*Vg;
vph*(1+VPZ*muP)*1/pZ - (1-tPd)*tauP*VPZ*muP*((pZ/pZ(-1)*PiY-PiY_bar)*PiY/pZ(-1)-Lam(+1)*(pZ(+1)/pZ*PiY(+1)-PiY_bar)/pZ*pZ(+1)/pZ*PiY(+1)*(Zd(+1)+Ze(+1))*g(+1)/(Zd+Ze)) = 1-tPd;
(1-tPd)*(1+tSd)*pL = vph*alpSZ^(1/sigSZ)*((Zd+Ze)/SZ)^(1/sigSZ)*(1-alpZ)^(1/sigZ)*VZ^((sigZ-1)/sigZ)*(SZ/LZ)^(1/sigZ);
1 = QZ*VK - (1-tPd)*QZ*VK*(tauI/2*(IPZ/IPZ(-1)-1)^2+tauI*(IPZ/IPZ(-1)-1)*IPZ/IPZ(-1)) + Lam(+1)*(1-tPd)*QZ(+1)*VK(+1)*tauI*(IPZ(+1)/IPZ-1)*(IPZ(+1)/IPZ)^2*g(+1);
(1-tPd)*(tauZ1+tauZ2*(uZ-1))*KZ(-1)/g = vph*alpSZ^(1/sigSZ)*((Zd+Ze)/SZ)^(1/sigSZ)*alpZ^(1/sigZ)*VZ^((sigZ-1)/sigZ)*(KZ(-1)/g)^((sigZ-1)/sigZ)*(SZ/uZ)^(1/sigZ);
QZ = Lam(+1)*(vph(+1)*alpSZ^(1/sigSZ)*((Zd(+1)+Ze(+1))/SZ(+1))^(1/sigSZ)*alpZ^(1/sigZ)*VZ(+1)^((sigZ-1)/sigZ)*uZ(+1)^((sigZ-1)/sigZ)*(SZ(+1)/(KZ/g(+1)))^(1/sigZ)*alpKZ^(1/sigKZ)*(KZ/KPZ)^(1/sigKZ) + (1-(1-tPd)*delZ)*QZ(+1)) + lam*(1-delZ)*PiY(+1)*(1/RB-1/RS);
(1-tPd)*pE = vph*(1-alpSZ)^(1/sigSZ)*((Zd+Ze)/EZ)^(1/sigSZ);
KPZ = VK*(1-tauI/2*(IPZ/IPZ(-1)-1)^2)*IPZ + (1-delZ)*KPZ(-1)/g;
Zd+Ze = (alpSZ^(1/sigSZ)*SZ^((sigSZ-1)/sigSZ)+(1-alpSZ)^(1/sigSZ)*EZ^((sigSZ-1)/sigSZ))^(sigSZ/(sigSZ-1));
SZ = VZ*(alpZ^(1/sigZ)*(uZ*KZ(-1)/g)^((sigZ-1)/sigZ)+(1-alpZ)^(1/sigZ)*LZ^((sigZ-1)/sigZ))^(sigZ/(sigZ-1));
KZ = (alpKZ^(1/sigKZ)*KPZ^((sigKZ-1)/sigKZ)+(1-alpKZ)^(1/sigKZ)*KGZ^((sigKZ-1)/sigKZ))^(sigKZ/(sigKZ-1));
DZ = pZ*(Zd+Ze) - pL*LZ - pE*EZ - TPZ - TSZ - IPZ + BZ - RB(-1)/PiY*BZ(-1)/g;
BZ = lam*QZ(+1)*(1-delZ)*KPZ*(RB/PiY(+1))^(-1);

% Labor firm (combining homogenized high and low skilled labor)
LH = VL^(sigL-1)*alpL*(pLH/pL)^(-sigL)*(LZ+LR+LF);
LL = VL^(sigL-1)*(1-alpL)*(pLL/pL)^(-sigL)*(LZ+LR+LF);
pL = VL^(-1)*(alpL*pLH^(1-sigL) + (1-alpL)*pLL^(1-sigL))^(1/(1-sigL));
pLH/pLH(-1)*g = PiH/PiY;
pLL/pLL(-1)*g = PiL/PiY;

% Energy firms (combining renewable energy and fossil energy)
ER = VERT^(sigE-1)*alpE*(pER/pE)^(-sigE)*(EY+EZ);
EF = VEFT^(sigE-1)*(1-alpE)*((pEF+tEFd)/pE)^(-sigE)*(EY+EZ);
pE = (alpE*VERT^(sigE-1)*pER^(1-sigE) + (1-alpE)*VEFT^(sigE-1)*(pEF+tEFd)^(1-sigE))^(1/(1-sigE));

% Renewable energy firms
(1+tSd)*pL = pER*(1-alpER)^(1/sigER)*VER^((sigER-1)/sigER)*(ER/LR)^(1/sigER);
1 = QR*VK - (1-tPd)*QR*VK*(tauI/2*(IPR/IPR(-1)-1)^2+tauI*(IPR/IPR(-1)-1)*IPR/IPR(-1)) + Lam(+1)*(1-tPd)*QR(+1)*VK(+1)*tauI*(IPR(+1)/IPR-1)*(IPR(+1)/IPR)^2*g(+1);
QR = Lam(+1)*((1-tPd)*pER(+1)*alpER^(1/sigER)*VER(+1)^((sigER-1)/sigER)*(ER(+1)/(KR/g(+1)))^(1/sigER)*alpKR^(1/sigKR)*(KR/KPR)^(1/sigKR) + (1-(1-tPd)*delR)*QR(+1)) + lam*(1-delR)*PiY(+1)*(1/RB-1/RS);
KPR = VK*(1-tauI/2*(IPR/IPR(-1)-1)^2)*IPR + (1-delR)*KPR(-1)/g;
ER = VER*(alpER^(1/sigER)*(KR(-1)/g)^((sigER-1)/sigER)+(1-alpER)^(1/sigER)*LR^((sigER-1)/sigER))^(sigER/(sigER-1));
KR = (alpKR^(1/sigKR)*KPR^((sigKR-1)/sigKR)+(1-alpKR)^(1/sigKR)*KGR^((sigKR-1)/sigKR))^(sigKR/(sigKR-1));
DR = pER*ER - pL*LR - TPR - TSR - IPR + BR - RB(-1)/PiY*BR(-1)/g;
BR = lam*QR(+1)*(1-delR)*KPR*(RB/PiY(+1))^(-1);

% Fossil energy firms
(1+tSd)*pL = pEF*alpSF^(1/sigSF)*(EF/SF)^(1/sigSF)*(1-alpEF)^(1/sigEF)*VEF^((sigEF-1)/sigEF)*(SF/LF)^(1/sigEF);
1 = QF*VK - (1-tPd)*QF*VK*(tauI/2*(IPF/IPF(-1)-1)^2+tauI*(IPF/IPF(-1)-1)*IPF/IPF(-1)) + Lam(+1)*(1-tPd)*QF(+1)*VK(+1)*tauI*(IPF(+1)/IPF-1)*(IPF(+1)/IPF)^2*g(+1);
QF = Lam(+1)*((1-tPd)*pEF(+1)*alpSF^(1/sigSF)*(EF(+1)/SF(+1))^(1/sigSF)*alpEF^(1/sigEF)*VEF(+1)^((sigEF-1)/sigEF)*((SF(+1))/(KF/g(+1)))^(1/sigEF)*alpKF^(1/sigKF)*(KF/KPF)^(1/sigKF) + (1-(1-tPd)*delF)*QF(+1)) + lam*(1-delF)*PiY(+1)*(1/RB-1/RS);
tOd + RER*pO_RoW + tauO*(O/O(-1)-1)*1/g - Lam(+1)*tauO*(O(+1)/O-1)*O(+1)/O = pEF*(1-alpSF)^(1/sigSF)*(EF/O)^(1/sigSF);
KPF = VK*(1-tauI/2*(IPF/IPF(-1)-1)^2)*IPF + (1-delF)*KPF(-1)/g;
EF = (alpSF^(1/sigSF)*SF^((sigSF-1)/sigSF)+(1-alpSF)^(1/sigSF)*O^((sigSF-1)/sigSF))^(sigSF/(sigSF-1));
SF = VEF*(alpEF^(1/sigEF)*(KF(-1)/g)^((sigEF-1)/sigEF)+(1-alpEF)^(1/sigEF)*LF^((sigEF-1)/sigEF))^(sigEF/(sigEF-1));
KF = (alpKF^(1/sigKF)*KPF^((sigKF-1)/sigKF)+(1-alpKF)^(1/sigKF)*KGF^((sigKF-1)/sigKF))^(sigKF/(sigKF-1));
DF = pEF*EF - pL*LF - RER*pO_RoW*O - TPF - TSF - TO - IPF + BF - RB(-1)/PiY*BF(-1)/g;
BF = lam*QF(+1)*(1-delF)*KPF*(RB/PiY(+1))^(-1);

% Policy
PiL/PiL_bar = (PiL(-1)/PiL_bar)^rhoWL*(ELe/ELe_bar)^(phiWE*(1-rhoWL))*(PiY(-1)/PiY_bar)^(phiWP*(1-rhoWL))*VWL;
ELe/ELe_bar = (ELe(+1)/ELe_bar)^(bet/(1+bet))*(ELe(-1)/ELe_bar)^(1/(1+bet))*((LL/NL)/(ELe/ELe_bar))^((1-bet*xiELe)*(1-xiELe)/((1+bet)*xiELe));
ELi=1;
T = TCH + TCL + TA + TR + TH + TL + TSZ + TSR + TSF + TPZ + TPR + TPF + TEF + TO;
TCH = tCd*CH;
TCL = tCd*CL;
TA = tAd*A;
TR = tRd*((RB(-1)-1)/PiY*(BG(-1)+BP(-1))/g + (RW(-1)-1)/PiY_RoW*RER*BW(-1)/g + (RS(-1)-1)/PiY*(pSZ(-1)+pSR(-1)+pSF(-1))/g);
TH = tHd*pLH*LH;
TL = tLd*pLL*LL;
TS = TSZ+TSR+TSF;
TSZ = tSd*pL*LZ;
TSR = tSd*pL*LR;
TSF = tSd*pL*LF;
TPZ = tPd*(pZ*(Zd+Ze) - pL*LZ - tSd*pL*LZ - QZ*delZ*KPZ(-1)/g - VK*tauI/2*(IPZ/IPZ(-1)-1)^2*IPZ);
TPR = tPd*(pER*ER - pL*LR - tSd*pL*LR - delR*KPR(-1)/g - VK*tauI/2*(IPR/IPR(-1)-1)^2*IPR);
TPF = tPd*(pEF*EF - pL*LF - tSd*pL*LF - TO - RER*pO_RoW*O - delF*KPF(-1)/g - VK*tauI/2*(IPF/IPF(-1)-1)^2*IPF);
TEF = tEFd*EF;
TO = tOd*O;
CG/CG_bar = (GDP/GDP_bar)^(-phiCY)*VCG;
IGZ/IGZ_bar = (GDP/GDP_bar)^(-phiIY)*VIGZ;
IGR/IGR_bar = (GDP/GDP_bar)^(-phiIY)*VIGR;
IGF/IGF_bar = (GDP/GDP_bar)^(-phiIY)*VIGF;
KGZ = VK*(1-tauI/2*(IGZ/IGZ(-1)-1)^2)*IGZ + (1-delZ)*KGZ(-1)/g;
KGR = VK*(1-tauI/2*(IGR/IGR(-1)-1)^2)*IGR + (1-delR)*KGR(-1)/g;
KGF = VK*(1-tauI/2*(IGF/IGF(-1)-1)^2)*IGF + (1-delF)*KGF(-1)/g;
CG + IGZ + IGR + IGF + TLL + RB(-1)/PiY*BG(-1)/g = T + BG;
R/R_bar = (R(-1)/R_bar)^rhoR*(pW/pW(-1)*PiY/PiY_bar)^(phiRP*(1-rhoR))*(GDP/GDP_bar)^(phiRY*(1-rhoR))*((pW/pW(-1)*PiY)/(pW(-1)/pW(-2)*PiY(-1)))^phidRP*(GDP/GDP(-1))^phidRY*VR;
R = (1-GamB*VGamB^(1/GamB))*RB;

% Rest of the world
pZ/RER*Ze - (pZ_RoW*Ze_RoW + alpO*pO_RoW*O) = BW - RW(-1)/PiY_RoW*BW(-1)/g; 
Y_RoW = Y_RoW(-1)^rhoY_RoW*Y_RoW_bar^(1-rhoY_RoW)*VY_RoW;
Ze = VW_RoW^(sigW_RoW-1)*(1-alpW_RoW)*(((pZ+tauZe*(Ze/Ze(-1)-1)*1/g-Lam_RoW(+1)*tauZe*(Ze(+1)/Ze-1)*Ze(+1)/Ze)/RER)/pZ_RoW)^(-sigW_RoW)*Y_RoW;
Lam_RoW = (RW/PiY_RoW)^(-1);
pZ_RoW = pZ_RoW(-1)^rhopZ_RoW*pZ_RoW_bar^(1-rhopZ_RoW)*VpZ_RoW;
RW/RW_bar = (RW(-1)/RW_bar)^rhoRW*VR_RoW;
PiY_RoW = PiY_RoW(-1)^rhoPiY_RoW*PiY_RoW_bar^(1-rhoPiY_RoW)*VPiY_RoW;
pO_RoW = pO_RoW(-1)^rhopO_RoW*pO_RoW_bar^(1-rhopO_RoW)*VpO_RoW;

% Market clearing
Y = CH + CL + CG + IPZ + IPR + IPF + IGZ + IGR + IGF;
BP = BZ+BR+BF;

% Indicator variables
GDP = Y+pZ*Ze-RER*(pZ_RoW*Ze_RoW + alpO*pO_RoW*O);

% Exogenous processes
VC      = VC(-1)      ^rhoC     *exp(epsC);
VN      = VN(-1)      ^rhoN     *exp(epsN);
VGamB   = VGamB(-1)   ^rhoRB    *exp(epsRB);
VGamW   =                        exp(epsGamW);
VGamS   = VGamS(-1)   ^rhoRS    *exp(epsRS);
VWH     = VWH(-1)     ^rhoWH    *exp(epsWH);
VWL     =                        exp(epsWL);
VL      = VL(-1)      ^rhoL     *exp(epsL)     *1^(1-rhoL);
VY      = VY(-1)      ^rhoY     *exp(epsY)     *1^(1-rhoY);
VW      = VW(-1)      ^rhoW     *exp(epsW)     *1^(1-rhoW);
VWI     = VWI(-1)     ^rhoWI    *exp(epsWI)    *1^(1-rhoWI);
VW_RoW  = VW_RoW(-1)  ^rhoW_RoW *exp(epsW_RoW) *1^(1-rhoW_RoW);
VZ      = VZ(-1)      ^rhoZ     *exp(epsZ)     *1^(1-rhoZ);
VEFT    = VEFT(-1)    ^rhoEFT   *exp(epsEFT)   *1^(1-rhoEFT);
VER     = VER(-1)     ^rhoER    *exp(epsER)    *1^(1-rhoER);
VERT    = VERT(-1)    ^rhoERT   *exp(epsERT)   *1^(1-rhoERT);
VEF     = VEF(-1)     ^rhoEF    *exp(epsEF)    *1^(1-rhoEF);
VPZ     = VPZ(-1)     ^rhoPZ    *exp(epsPZ);
VK      = VK(-1)      ^rhoK     *exp(epsK);
Vg      =                        exp(epsg);
VCG     = VCG(-1)     ^rhoCG    *exp(epsCG);
VIGZ    = VIGZ(-1)    ^rhoIG    *exp(epsIG);
VIGR    = VIGR(-1)    ^rhoIGR   *exp(epsIGR);
VIGF    = VIGF(-1)    ^rhoIG    *exp(epsIG);
VR      =                        exp(epsR);
VR_RoW  =                        exp(epsRW);
VpZ_RoW =                        exp(epspZ_RoW);
VY_RoW  =                        exp(epsY_RoW);
VPiY_RoW=                        exp(epsPiY_RoW);
VpO_RoW =                        exp(epspO_RoW);



tOd_add = 0.99999*tOd_add(-1) + tOd_shock;
tOd = tO + tOd_add;

tEFd_add = 0.99999*tEFd_add(-1) + tEFd_shock;
tEFd = tEF + tEFd_add;

tSd_add = 0.99999*tSd_add(-1) + tSd_shock;
tSd = tS + tSd_add;

tHd_add = 0.99999*tHd_add(-1) + tHd_shock;
tHd = tH + tHd_add;

tLd_add = 0.99999*tLd_add(-1) + tLd_shock;
tLd = tL + tLd_add;

tCd_add = 0.99999*tCd_add(-1) + tCd_shock;
tCd = tC + tCd_add;

tPd_add = 0.99999*tPd_add(-1) + tPd_shock;
tPd = tP + tPd_add;

tRd_add = 0.99999*tRd_add(-1) + tRd_shock;
tRd = tR + tRd_add;

tAd_add = 0.99999*tAd_add(-1) + tAd_shock;
tAd = tA + tAd_add;

taxbase_TC = (TCH(1)+TCL(1))/tC;
taxbase_TA = TA(1)/tA;
taxbase_TR = TR(1)/tR;
taxbase_TH = TH(1)/tH;
taxbase_TL = TL(1)/tL;
taxbase_TS = TS(1)/tS;
taxbase_TP = (TPZ(1)+TPR(1)+TPF(1))/tP;
taxbase_TEF = TEF(1)/tEF;
taxbase_TO = TO(1)/tO;

CAPEXshare_Z = 100*(IPZ+IGZ)/(pZ*(Zd+Ze));
CAPEXshare_R = 100*(IPR+IGR)/(pER*ER);
CAPEXshare_F = 100*(IPF+IGF)/(pEF*EF);
OPEXshare_Z = 100*(pL*LZ + pE*EZ + TPZ + TSZ - BZ + RB/PiY*BZ/g)/(pZ*(Zd+Ze));
OPEXshare_R = 100*(pL*LR + TPR + TSR - BR + RB/PiY*BR/g)/(pER*ER);
OPEXshare_F = 100*(pL*LF + RER*pO_RoW*O + TPF + TSF + TO - BF + RB/PiY*BF/g)/(pEF*EF);
DIVshare_Z = 100*DZ/(pZ*(Zd+Ze));
DIVshare_R = 100*DR/(pER*ER);
DIVshare_F = 100*DF/(pEF*EF);


% GDP and components (percent change from pre-policy level)
 lGDP    = 100*log(GDP);
 lConP   = 100*log(CH+CL);
 lConH   = 100*log(CH);
 lConL   = 100*log(CL);
 lIVP    = 100*log(IPZ+IPR+IPF);
 lIVG    = 100*log(IGZ+IGR+IGF);
 lExp    = 100*log(pZ*Ze);
 lImp    = 100*log(RER*(pZ_RoW*Ze_RoW + alpO*pO_RoW*O));
 lImp_WO = 100*log(RER*pZ_RoW*Ze_RoW);
 ppTB    = 100*(pZ*Ze - RER*(pZ_RoW*Ze_RoW + alpO*pO_RoW*O)); 


% GDP and components (normalized by GDP under pre-policy prices)
 ppU        = 100*(1-ELe);
 ConP_GDP   = 100*(CH+CL);
 ConH_GDP   = 100*(CH);
 ConL_GDP   = 100*(CL);
 IVP_GDP    = 100*(IPZ+IPR+IPF);
 IVG_GDP    = 100*(IGZ+IGR+IGF);
 IVG_Z_GDP  = 100*(IGZ);
 IVG_R_GDP  = 100*(IGR);
 IVG_F_GDP  = 100*(IGF); 
 Exp_GDP    = 100*(pZ*Ze);
 Imp_GDP    = 100*(RER*(pZ_RoW*Ze_RoW + alpO*pO_RoW*O));

% Sectoral output (percent change from pre-policy level)
 lZ  = 100*log(Zd+Ze);
 lER = 100*log(ER);
 lEF = 100*log(EF);
 lEY = 100*log(EY);
 lEZ = 100*log(EZ);
 lE  = 100*log(EY+EZ);
 lW  = 100*log(W);
 lLZ = 100*log(LZ);
 lLR = 100*log(LR);
 lLF = 100*log(LF);

 % Sectoral output (normalized by GDP under pre-policy prices)
 Z_GDP  = 100*(Zd+Ze);
 ER_GDP = 100*(ER);
 EF_GDP = 100*(EF);
 EY_GDP = 100*(EY);
 EZ_GDP = 100*(EZ);
 E_GDP  = 100*(EY+EZ);
 W_GDP  = 100*(W);


% prices
 lpZ     = 100*log(pZ);
 lpER    = 100*log(pER);
 lpEF    = 100*log(pEF);
 lpE     = 100*log(pE);
 lpW     = 100*log(pW);
 lpZ_RoW = 100*log(RER*pZ_RoW);
 lpLH    = 100*log(pLH);
 lpLL    = 100*log(pLL);
 lpL     = 100*log(pL);

% private Investment (percent change from pre-policy level)
 lIPZ = 100*log(IPZ);
 lIPR = 100*log(IPR);
 lIPF = 100*log(IPF);
 lKPZ = 100*log(KPZ(-1));
 lKPR = 100*log(KPR(-1));
 lKPF = 100*log(KPF(-1));
 lpSZ = 100*log(pSZ);
 lpSR = 100*log(pSR);
 lpSF = 100*log(pSF);
 ppuZ = 100*(uZ);

% Tax revenues

ppT  = 100*T;
ppTO = 100*TO;
ppTEF= 100*TEF;
ppTS = 100*TS;
ppTL = 100*TL;
ppTH = 100*TH;
ppTC = 100*(TCH+TCL);
ppTP = 100*(TPZ+TPR+TPF);
ppTR = 100*TR;
ppTA = 100*TA;

% Prices, inflation, emissions

ppPiY   = 400*PiY;
ppR     = 400*R;
ppRP    = 400*(R/PiY(+1));
ppRER   = 100*RER;
lO      = 100*log(O);
lptO    = 100*log(tOd+RER*pO_RoW);
ppAY    = 100*A/GDP/4;
ppBGY   = 100*BG/GDP/4;
ppBPY   = 100*BP/GDP/4;
ppBWY   = 100*BW/GDP/4;
ppGamW  = 100*GamW;


Con=CH+CL;
IG=IGR+IGF+IGZ;
IP=IPR+IPF+IPZ;
Inv=IG+IP;
Ex=pZ*Ze;
Im=RER*(pZ_RoW*Ze_RoW + alpO*pO_RoW*O);
TB=(pZ*Ze-RER*(pZ_RoW*Ze_RoW + alpO*pO_RoW*O)); 
E=EY+EZ;
EF_ann = EF+EF(-1)+EF(-2)+EF(-3);
ER_ann = ER+ER(-1)+ER(-2)+ER(-3);
EY_ann = EY+EY(-1)+EY(-2)+EY(-3);
EZ_ann = EZ+EZ(-1)+EZ(-2)+EZ(-3);
O_ann  = O+O(-1)+O(-2)+O(-3);

% Measurement error

GDP_m=GDP*exp(epsGDP_m);

% Observed variables

dlnGDP=(GDP_m/GDP_m(-1)-1)*100;
dlnC=(Con/Con(-1)-1)*100;
dlnIG=(IG/IG(-1)-1)*100;
dlnIP=(IP/IP(-1)-1)*100;
dlnCG=(CG/CG(-1)-1)*100;
dlnEx=(Ex/Ex(-1)-1)*100;
dlnIm=(Im/Im(-1)-1)*100;
dlnEL=(ELe/ELe(-1)-1)*100;
dlnpL=(pL/pL(-1)-1)*100;
dlnE=(E/E(-1)-1)*100;
dlnpE=(pE/pE(-1)-1)*100;
dlnpZ=(pZ/pZ(-1)-1)*100;
dlnpEF=((pEF+tEFd)/(pEF(-1)+tEFd(-1))-1)*100;
dlnY_RoW=(Y_RoW/Y_RoW(-1)-1)*100;
dlnpZ_RoW=(pZ_RoW/pZ_RoW(-1)-1)*100;
dlnpO_RoW=(pO_RoW/pO_RoW(-1)-1)*100;
dlnEF=(EF/EF(-1)-1)*100;
dlnER=(ER/ER(-1)-1)*100;
dlnEY=(EY/EY(-1)-1)*100;
dlnEZ=(EZ/EZ(-1)-1)*100;
dlnO=(O/O(-1)-1)*100;

dlnEF_ann=dlnEF+dlnEF(-1)+dlnEF(-2)+dlnEF(-3);
dlnER_ann=dlnER+dlnER(-1)+dlnER(-2)+dlnER(-3);
dlnEY_ann=dlnEY+dlnEY(-1)+dlnEY(-2)+dlnEY(-3);
dlnEZ_ann=dlnEZ+dlnEZ(-1)+dlnEZ(-2)+dlnEZ(-3);
dlnO_ann=dlnO+dlnO(-1)+dlnO(-2)+dlnO(-3);

devGDP=(log(GDP_m)-log(GDP_bar))*100;
devC=(log(Con)-log(Cbar))*100;
devIG=(log(IG)-log(IGbar))*100;
devIP=(log(IP)-log(IPbar))*100;
devCG=(log(CG)-log(CG_bar))*100;
devEx=(log(Ex)-log(Exbar))*100;
devIm=(log(Im)-log(Imbar))*100;
devEL=(log(ELe)-log(ELe_bar))*100;
devpL=(log(pL)-log(pLbar))*100;
devPiY=(PiY-PiY_bar)*100;
devPiY_RoW=(PiY_RoW-PiY_RoW_bar)*100;
devR=(R-R_bar)*100;
devRB=(RB-RBbar)*100;
devRW=(RW-RW_bar)*100;
devRER=(RER-RERbar)*100;
devpE=(log(pE)-log(pEbar))*100;
devpZ=(log(pZ)-log(pZbar))*100;
devpEF=(log(pEF)-log(pEFbar))*100;
devY_RoW=(log(Y_RoW)-log(Y_RoW_bar))*100;
devpZ_RoW=(log(pZ_RoW)-log(pZ_RoW_bar))*100;
devpO_RoW=(log(pO_RoW)-log(pO_RoW_bar))*100;
devE=(log(E)-log(Ebar))*100;
devu=psiu*(log(uZ)-log(uZ_bar))*100;

devEF_ann=(log(EF_ann)-log(EFbar))*100;
devER_ann=(log(ER_ann)-log(ERbar))*100;
devEY_ann=(log(EY_ann)-log(EYbar))*100;
devEZ_ann=(log(EZ_ann)-log(EZbar))*100;
devO_ann=(log(O_ann)-log(Obar))*100;

end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

steady;
check;

% ESTIMATION
  
varobs dlnGDP dlnC dlnIG dlnIP dlnCG devR devRB devPiY dlnpL dlnEx dlnIm dlnEL dlnpZ dlnpEF dlnpE devRER dlnEF_ann dlnER_ann dlnEY_ann dlnEZ_ann dlnO_ann dlnpZ_RoW dlnpO_RoW dlnY_RoW devu devPiY_RoW;

estimated_params;

kap,0.53,,,               beta_pdf,              0.5,   0.15;
vrh,1,,,                  gamma_pdf,             1.5,   0.5;
etaH,1.9,,,               gamma_pdf,             1.75,  0.5;
sigSZ,1.94,,,             gamma_pdf,             0.9,   0.15;
sigL,1.63,,,              gamma_pdf,             2,     0.5;
sigY,0.03,,,              gamma_pdf,             0.3,   0.1;
sigSF,0.27,,,             gamma_pdf,             0.4,   0.15;
sigKZ,                    gamma_pdf,             3,     0.25;
sigZ,0.9,,,               gamma_pdf,             0.9,   0.05;
sigER,0.28,,,             gamma_pdf,             0.3,   0.1;
sigEF,0.38,,,             gamma_pdf,             0.3,   0.1;
sigW,1.84,,,              gamma_pdf,             2,     0.25;
sigW_RoW,1.95,,,          gamma_pdf,             2,     0.25;
sigE,14,,,                gamma_pdf,             6,     4;
tauW,37,,,                gamma_pdf,             50,    20;
tauP,85,,,                gamma_pdf,             50,    20;
tauI,4.5,,,               gamma_pdf,             5,     1;
tauZ2,0.35,,,             gamma_pdf,             0.5,   0.1;
tauO,2.75,,,              gamma_pdf,             5,     1;
tauZe,4.6,,,              gamma_pdf,             5,     1;
psiu,3.64,,,              normal_pdf,            2,     1;
xiELe,0.95,,,             beta_pdf,              0.5,   0.15;
phiWE,0.20,,,             normal_pdf,            0.4,   0.25;
phiWP,0.20,,,             normal_pdf,            0.4,   0.25;
phiCY,-0.15,,,            normal_pdf,            0,     0.1;
phiIY,-0.03,,,            normal_pdf,            0,     0.1;
phiRP,                    normal_pdf,            1.7,   0.1;
phidRP,-0.03,,,           normal_pdf,            0.3,   0.1;
phiRY,0.2,,,              normal_pdf,            0.125, 0.05;
phidRY,0.03,,,            normal_pdf,            0.0625,0.05;
rhog,0.65,,,              beta_pdf,              0.5,   0.1;
rhoW,0.75,,,              beta_pdf,              0.5,   0.1;
rhoWI,0.74,,,             beta_pdf,              0.5,   0.1;
rhoW_RoW,0.33,,,          beta_pdf,              0.5,   0.1;
rhoR,0.91,,,              beta_pdf,              0.8,   0.1;
%rhoRW,0.55,,,             beta_pdf,              0.5,   0.05;
rhoC,0.43,,,              beta_pdf,              0.5,   0.1;
rhoN,0.43,,,              beta_pdf,              0.5,   0.1;
rhoRB,0.4,,,              beta_pdf,              0.5,   0.1;
%rhoRS,                    beta_pdf,              0.5,   0.05;
rhoZ,0.73,,,              beta_pdf,              0.5,   0.1;
rhoEF,0.75,,,             beta_pdf,              0.5,   0.1;
rhoER,0.78,,,             beta_pdf,              0.5,   0.1;
rhoERT,0.81,,,            beta_pdf,              0.5,   0.1;
rhoK,0.57,,,              beta_pdf,              0.5,   0.1;
rhoCG,0.63,,,             beta_pdf,              0.5,   0.1;
rhoIG,0.70,,,             beta_pdf,              0.5,   0.1;
rhoY_RoW,0.89,,,          beta_pdf,              0.75,   0.05;
rhopZ_RoW,0.86,,,         beta_pdf,              0.75,   0.05;
rhoWL,0.53,,,             beta_pdf,              0.5,   0.1;
%rhoWH,                    beta_pdf,              0.5,   0.1;
rhoEFT,0.75,,,            beta_pdf,              0.5,   0.1;
rhoL,0.63,,,              beta_pdf,              0.5,   0.1;
rhoY,0.89,,,              beta_pdf,              0.5,   0.1;
rhopO_RoW,0.86,,,         beta_pdf,              0.75,   0.05;
rhoPiY_RoW,0.44,,,        beta_pdf,              0.5,   0.05;
stderr epsg,0.03,,,       inv_gamma_pdf,         0.1,   2;
stderr epsWL,0.04,,,      inv_gamma_pdf,         0.1,   2;
%stderr epsWH,             inv_gamma_pdf,         0.1,   2;
stderr epsR,0.005,,,      inv_gamma_pdf,         0.01,  0.1;
%stderr epsRW,0.0005,,,    inv_gamma_pdf,         0.001, 0.01;
stderr epsC,0.08,,,       inv_gamma_pdf,         0.1,   2;
stderr epsN,0.96,,,       inv_gamma_pdf,         0.1,   2;
stderr epsRB,0.02,,,      inv_gamma_pdf,         0.01,  0.1;
%stderr epsRS,             inv_gamma_pdf,         0.1,   2;
stderr epsPZ,1.41,,,      inv_gamma_pdf,         0.1,   2;
stderr epsZ,0.07,,,       inv_gamma_pdf,         0.1,   2;
stderr epsER,0.16,,,      inv_gamma_pdf,         0.1,   2;
stderr epsERT,0.79,,,     inv_gamma_pdf,         0.1,   2;
stderr epsEF,0.25,,,      inv_gamma_pdf,         0.1,   2;
stderr epsEFT,0.07,,,     inv_gamma_pdf,         0.1,   2;
stderr epsL,0.02,,,       inv_gamma_pdf,         0.1,   2;
stderr epsK,0.18,,,       inv_gamma_pdf,         0.1,   2;
stderr epsCG,0.04,,,      inv_gamma_pdf,         0.1,   2;
stderr epsIG,             inv_gamma_pdf,         0.1,   2;
stderr epsY_RoW,0.006,,,  inv_gamma_pdf,         0.01,  0.1;
stderr epspZ_RoW,0.007,,, inv_gamma_pdf,         0.01,  0.1;
stderr epsGamW,0.3,,,     inv_gamma_pdf,         0.1,   2;
stderr epsY,0.03,,,       inv_gamma_pdf,         0.1,   2;
stderr epsW,0.04,,,       inv_gamma_pdf,         0.1,   2;
stderr epsWI,0.07,,,      inv_gamma_pdf,         0.1,   2;
stderr epsW_RoW,0.37,,,   inv_gamma_pdf,         0.1,   2;
stderr epspO_RoW,0.08,,,  inv_gamma_pdf,         0.1,   2;
stderr epsPiY_RoW,0.006,,,inv_gamma_pdf,         0.01,  0.1;
stderr epsGDP_m,0.02,,,   inv_gamma_pdf,         0.01,  0.1;

end;

identification(ar=1);

estimation(datafile=data,mode_compute = 5,mh_replic=50000,mode_check,mh_jscale=0.225,bayesian_irf) devGDP devC devIG devIP devCG devR devRB devRW devPiY devpL devEL devEx devIm devRER devpE devpZ devpEF devY_RoW devpZ_RoW devpO_RoW devEF_ann devER_ann devEY_ann devEZ_ann devO_ann devu devPiY_RoW;

%estimation(datafile=data,mode_compute =0,mode_file='estimation_mode.mat', mh_replic=50000,mh_jscale=0.225,bayesian_irf) devGDP devC devIG devIP devCG devR devRB devRW devPiY devpL devEL devEx devIm devRER devpE devpZ devpEF devE devY_RoW devpZ_RoW devpO_RoW devEF_ann devER_ann devEY_ann devEZ_ann devO_ann devu devPiY_RoW;

%estimation(datafile=data,mode_compute =0,mode_file='estimation_mode.mat',load_mh_file,mh_replic=0,nograph,nodisplay) devGDP devC devIG devIP devCG devR devRB devRW devPiY devpL devEL devEx devIm devRER devpE devpZ devpEF devE devY_RoW devpZ_RoW devpO_RoW devEF_ann devER_ann devEY_ann devEZ_ann devO_ann devu devPiY_RoW;

%smoother2histval(outfile = 'Smoothed_Variables.mat');

%conditional_forecast_paths;
%var O;
%periods 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132;
%values 0.4375,0.43,0.4235,0.4185,0.415,0.4125,0.41,0.4075,0.405,0.4025,0.40,0.3975,0.395,0.3925,0.39,0.3875,0.385,0.3825,0.38,0.3775,0.375,0.3725,0.37,0.3675,0.365,0.3625,0.36,0.3575,0.355,0.3525,0.35,0.3475,0.345,0.3425,0.34,0.3375,0.335,0.3325,0.33,0.3275,0.325,0.3225,0.32,0.3175,0.315,0.3125,0.31,0.3075,0.305,0.3025,0.3,0.2975,0.295,0.2925,0.29,0.2875,0.285,0.2825,0.28,0.2775,0.275,0.2725,0.27,0.2675,0.265,0.2625,0.26,0.2575,0.255,0.2525,0.25,0.2475,0.245,0.2425,0.24,0.2375,0.235,0.2325,0.23,0.2275,0.225,0.2225,0.22,0.2175,0.215,0.2125,0.21,0.2075,0.205,0.2025,0.2,0.1975,0.195,0.1925,0.19,0.1875,0.185,0.1825,0.18,0.1775,0.175,0.1725,0.17,0.1675,0.165,0.1625,0.16,0.1575,0.155,0.1525,0.15,0.1475,0.145,0.1425,0.14,0.1375,0.135,0.1325,0.13,0.1275,0.125,0.1225,0.12,0.1175,0.115,0.1125,0.11,0.1075,0.1050,0.1025,0.10,0.0975;
%periods 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40;
%values 0.4375,0.43,0.4225,0.4175,0.415,0.4125,0.41,0.4075,0.405,0.4025,0.40,0.3975,0.395,0.3925,0.39,0.3875,0.385,0.3825,0.38,0.3775,0.375,0.3725,0.37,0.3675,0.365,0.3625,0.36,0.3575,0.355,0.3525,0.35,0.3475,0.345,0.3425,0.34,0.3375,0.335,0.3325,0.33,0.3275;
%values 0.4325,0.4225,0.415,0.4075,0.40,0.3925,0.385,0.3775,0.37,0.3625,0.355,0.3475,0.34,0.3325,0.325,0.3175,0.31,0.305,0.30,0.295,0.29,0.285,0.28,0.275,0.27,0.265,0.26,0.255,0.25,0.245,0.24,0.235,0.23,0.225,0.22,0.215,0.21,0.205,0.2,0.195; 
%var tOd;
%periods 1:40;
%values 0.1;
%var PiY;
%periods 1:40;
%values 1.0199;
%var BG;
%periods 1:40;
%values 1.3546;
%end;

%conditional_forecast(parameter_set=posterior_mean, controlled_varexo=(tOd_shock), periods=40, conf_sig=0.68);

%plot_conditional_forecast (periods=40) 
   
% GDP and components (normalized by GDP under pre-policy prices)
% O BG VIGR tOd lGDP ppU ConP_GDP ConH_GDP ConL_GDP IVP_GDP IVG_GDP Exp_GDP Imp_GDP IVG_R_GDP PiY PiL ppT

% Sectoral output (normalized by GDP under pre-policy prices)
% Z_GDP ER_GDP EF_GDP EY_GDP EZ_GDP E_GDP W_GDP       

% private Investment (percent change from pre-policy level)
% IPZ IPR IPF KPZ KPR KPF pSZ pSR pSF ppuZ 

% Tax revenues
% ppT ppTO ppTEF ppTS ppTL ppTH ppTC ppTP ppTR ppTA 

% Prices, inflation, emissions, wealth
% ppPiY ppR ppRP ppRER lO lptO ppAY ppBGY ppBPY ppBWY 
%;
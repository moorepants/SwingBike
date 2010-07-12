%---------------------------------------------------------------------%
%    FILE: SW_01.AL
%    DATE: SEPTEMBER 23, 2008
%  AUTHOR: JASON MOORE
% PROBLEM: GENERATES THE LINEARIZED EQUATIONS OF MOTION FOR A RIGID 
%          RIDER NO-HANDS SWING BICYCLE MADE UP OF FIVE RIGID BODIES.
%   NOTES: OUTPUTS CODE FOR NONLINEAR SIMULATION
%---------------------------------------------------------------------%
%         DEFAULT SETTINGS
%---------------------------------------------------------------------%

AUTOZ ON
AUTORHS OFF
OVERWRITE ALL
BEEPSOUND OFF

%---------------------------------------------------------------------%
%         NEWTONIAN, BODIES, FRAMES, PARTICLES, POINTS
%---------------------------------------------------------------------%

% DECLARE THE INERTIAL REFERENCE FRAME

NEWTONIAN N

% DECLARE FOUR INTERMEDIATE FRAMES
% A: YAW FRAME
% B: ROLL FRAME
% H: REAR HEAD ANGLE FRAME
% I: FRONT HEAD ANGLE FRAME

FRAMES A,B,H,I

% DECLARE FIVE BODIES
% C: REAR FRAME
% D: REAR WHEEL
% E: MID FRAME
% F: FRONT FORK
% G: FRONT WHEEL

BODIES C,D,E,F,G

% DECLARE FOUR POINTS
% ND: REAR CONTACT POINT ON GROUND
% DN: REAR CONTACT POINT ON WHEEL
% NG: FRONT CONTACT POINT ON GROUND
% GN: FRONT CONTACT POINT ON WHEEL

POINTS ND,DN,NG,GN

%---------------------------------------------------------------------%
%         CONSTANTS AND VARIABLES
%---------------------------------------------------------------------%

% W:        WHEELBASE                          [M]
% CF:       FRONT WHEEL TRAIL                  [M]
% CR:       REAR WHEEL TRAIL                   [M]
% H:        HEIGHT OF TOP TUBE                 [H] 
% LAMBDAF:  FRONT STEER AXIS TILT              [RAD]
% LAMBDAR:  REAR STEER AXIS TILT               [RAD]
% G:        GRAVITY                            [N/KG]
% V:        FORWARD SPEED                      [M/S]
% RR:       REAR WHEEL RADIUS                  [M]
% M_R:      REAR WHEEL MASS                    [KG]
% IRXX:     REAR WHEEL MASS MOMENT OF INERTIA  [KG*M^2]
% IRYY:     REAR WHEEL MASS MOMENT OF INERTIA  [KG*M^2]
% XB:       REAR BODY CENTER OF MASS LOCATION  [M]
% ZB:       REAR BODY CENTER OF MASS LOCATION  [M]
% M_B:      REAR BODY MASS                     [KG]
% IBXX:     REAR BODY MASS MOMENT OF INERTIA   [KG*M^2]
% IBYY:     REAR BODY MASS MOMENT OF INERTIA   [KG*M^2]
% IBZZ:     REAR BODY MASS MOMENT OF INERTIA   [KG*M^2]
% IBXZ:     REAR BODY MASS PRODUCT OF INERTIA  [KG*M^2]
% XE:       MID FRAME CENTER OF MASS LOCATION  [M]
% ZE:       MID FRAME CENTER OF MASS LOCATION  [M]
% M_E:      MID FRAME MASS                     [KG]
% IEXX:     MID FRAME MASS MOMENT OF INERTIA   [KG*M^2]
% IEYY:     MID FRAME MASS MOMENT OF INERTIA   [KG*M^2]
% IEZZ:     MID FRAME MASS MOMENT OF INERTIA   [KG*M^2]
% IEXZ:     MID FRAME MASS PRODUCT OF INERTIA  [KG*M^2]
% XH:       FORK CENTER OF MASS LOCATION       [M]
% ZH:       FORK CENTER OF MASS LOCATION       [M]
% M_H:      FORK MASS                          [KG]
% IHXX:     FORK MASS MOMENT OF INERTIA        [KG*M^2]
% IHYY:     FORK MASS MOMENT OF INERTIA        [KG*M^2]
% IHZZ:     FORK MASS MOMENT OF INERTIA        [KG*M^2]
% IHXZ:     FORK MASS PRODUCT OF INERTIA       [KG*M^2]
% RF:       FRONT WHEEL RADIUS                 [M]
% M_F:      FRONT WHEEL MASS                   [KG]
% IFXX:     FRONT WHEEL MASS MOMENT OF INERTIA [KG*M^2]
% IFYY:     FRONT WHEEL MASS MOMENT OF INERTIA [KG*M^2]
% T_PHI:    LEAN TORQUE                        [N*M]
% T_DELTAF: FRONT STEER TORQUE                 [N*M]
% T_DELTAR: REAR STEER TORQUE                  [N*M]
% T_THETAR: REAR WHEEL TORQUE                  [N*M]

CONSTANTS W,CF,CR,H,LAMBDAF,LAMBDAR,G,V
CONSTANTS RR,M_R,IRXX,IRYY
CONSTANTS XB,ZB,M_B,IBXX,IBYY,IBZZ,IBXZ
CONSTANTS XE,ZE,M_E,IEXX,IEYY,IEZZ,IEXZ
CONSTANTS XH,ZH,M_H,IHXX,IHYY,IHZZ,IHXZ
CONSTANTS RF,M_F,IFXX,IFYY
SPECIFIED T_PHI,T_DELTAF,T_DELTAR,T_THETAR

% CONVERT THE BENCHMARK CONSTANTS TO SWING BIKE'S CONSTANTS
% RF: RADIUS OF FRONT WHEEL
% RR: RADIUS OF REAR WHEEL
% D1: THE REAR FORK OFFSET
% D2: THE REAR HEABTUBE LENGTH
% D3: THE HORIZONTAL TOP TUBE LENGTH
% D4: THE FRONT HEABTUBE LENGTH
% D5: THE FRONT FORK OFFSET
% L1: THE DISTANCE IN THE C1> DIRECTION FROM THE CENTER OF THE REAR    
%     WHEEL TO THE REAR FRAME CENTER OF MASS
% L2: THE DISTANCE IN THE C3> DIRECTION FROM THE CENTER OF THE REAR    
%     WHEEL TO THE REAR FRAME CENTER OF MASS
% L3: THE DISTANCE IN THE E1> DIRECTION FROM THE SEAT TUBE/TOP TUBE
%     INTERSECTION TO THE MID FRAME CENTER OF MASS
% L4: THE DISTANCE IN THE E3> DIRECTION FROM THE SEAT TUBE/TOP TUBE
%     INTERSECTION TO THE MID FRAME CENTER OF MASS
% L5: THE DISTANCE IN THE F1> DIRECTION FROM THE HEAD TUBE/TOP TUBE
%     INTERSECTION TO THE CENTER OF MASS OF THE FORK
% L6: THE DISTANCE IN THE F3> DIRECTION FROM THE HEAD TUBE/TOP TUBE
%     INTERSECTION TO THE CENTER OF MASS OF THE FORK

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST SECTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%W=1.02
%CF=0.08
%CR=0.08
%H=0.8
%LAMBDAF=3.141592653589793/10
%LAMBDAR=3.141592653589793/10
%G=9.81
%RR=0.3
%M_R=2
%IRXX=0.0603
%IRYY=0.12
%XB=0.3
%ZB=-0.9
%M_B=85
%IBXX=9.2
%IBYY=11
%IBZZ=2.8
%IBXZ=2.4
%XH=0.9
%ZH=-0.7
%M_H=4
%IHXX=0.05892
%IHYY=0.06
%IHZZ=0.00708
%IHXZ=-0.00756
%RF=0.35
%M_F=3
%IFXX=0.1405
%IFYY=0.28
%M_E=3
%XE=0.59
%ZE=-0.8
%IEXX=0.0053
%IEYY=0.1
%IEZZ=0.1
%IEXZ=0.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END TEST SECTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

D1    = CR*COS(LAMBDAR)+RR*SIN(LAMBDAR)
D2    = (H-RR+D1*SIN(LAMBDAR))/COS(LAMBDAR)  
D5    = RF*SIN(LAMBDAF)-CF*COS(LAMBDAF)
D4    = (RR-D1*SIN(LAMBDAR)+D2*COS(LAMBDAR)+D5*SIN(LAMBDAF)-RF)/&
        COS(LAMBDAF)
D3    = W-D1*COS(LAMBDAR)-D2*SIN(LAMBDAR)-D4*SIN(LAMBDAF)-D5*&
        COS(LAMBDAF)

MC    =  M_B-M_E                     
MD    =  M_R
ME    =  M_E                    
MF    =  M_H                     
MG    =  M_F
          
XC    = (XB*M_B-XE*M_E)/MC
ZC    = (ZB*M_B-ZE*M_E)/MC

ID11  =  IRxx                   
ID22  =  IRyy                   
ID33  =  IRxx

% USE THE PARALLEL AXIS THEREOM TO SHIFT THE INERTIA OF E TO THE CG OF
% B AND SUBTRACT THE INERTIA OF E FROM B
A1    = ABS(XB-XE)
B1    = ABS(0)
C1    = ABS(ZB-ZE)
PA1   = [B1^2+C1^2,-A1*B1,-A1*C1;-A1*B1,C1^2+A1^2,-B1*C1;-A1*C1,&
         -B1*C1,A1^2+B1^2]
IE    = [IEXX,0,IEXZ;0,IEYY,0;IEXZ,0,IEZZ]
IEP   = IE + ME*PA1
IB    = [IBXX,0,IBXZ;0,IBYY,0;IBXZ,0,IBZZ]
% SUBTRACT THE SHIFTED INERTIA OF E FROM B
IB2   = IB-IEP
A2    = ABS(XB-XC)
B2    = ABS(0)
C2    = ABS(ZB-ZC)
PA2   = [B2^2+C2^2,-A2*B2,-A2*C2;-A2*B2,C2^2+A2^2,-B2*C2;-A2*C2,&
         -B2*C2,A2^2+B2^2]
IB2P  = IB2 + MC*PA2

% ROTATE REAR FRAME INERTIA MATRIX THROUGH LAMBDAR
RMR   =  [COS(LAMBDAR),0,-SIN(LAMBDAR);0,1,0;SIN(LAMBDAR),0,&
COS(LAMBDAR)]
IC    =  RMR*IB2P*TRANSPOSE(RMR)
IC11  =  IC[1,1]             
IC12  =  IC[1,2]             
IC22  =  IC[2,2]             
IC23  =  IC[2,3]             
IC31  =  IC[3,1]             
IC33  =  IC[3,3]

IE11  =  IEXX             
IE12  =  0             
IE22  =  IEYY             
IE23  =  0             
IE31  =  IEXZ             
IE33  =  IEZZ
IE    =  [IE11,IE12,IE31;IE23,IE22,IE23;IE31,IE23,IE33]                

% ROTATE FORK INERTIA MATRIX THROUGH LAMBDAF
IH    =  [IHXX,0,IHXZ;0,IHYY,0;IHXZ,0,IHZZ]
RMF   =  [COS(LAMBDAF),0,-SIN(LAMBDAF);0,1,0;SIN(LAMBDAF),0,&
COS(LAMBDAF)]
IHROT =  RMF*IH*TRANSPOSE(RMF)
                   
IF11  =  IHROT[1,1]             
IF12  =  IHROT[1,2]             
IF22  =  IHROT[2,2]             
IF23  =  IHROT[2,3]             
IF31  =  IHROT[3,1]             
IF33  =  IHROT[3,3]
IF    =  [IF11,IF12,IF31;IF23,IF22,IF23;IF31,IF23,IF33]             

IG11  =  IFxx                   
IG22  =  IFyy                   
IG33  =  IFxx
       
L2    = (RR-XC*TAN(LAMBDAR)+ZC)/(TAN(LAMBDAR)*SIN(LAMBDAR)+&
        COS(LAMBDAR))
L1    = XC/COS(LAMBDAR)+L2*TAN(LAMBDAR)
L3    = XE-D1*COS(LAMBDAR)-D2*SIN(LAMBDAR)
L4    = ZE+RR-D1*SIN(LAMBDAR)+D2*COS(LAMBDAR)
L5    = COS(LAMBDAF)*XH-SIN(LAMBDAF)*ZH-CF*COS(LAMBDAF)-W*COS(LAMBDAF)
L6    = COS(LAMBDAF)*(H+ZH) + SIN(LAMBDAF)^2*D4 + SIN(LAMBDAF)*&
        COS(LAMBDAF)*D5 - SIN(LAMBDAF)*(W-XH)
                     
T4    =  T_PHI
T6    =  T_THETAR
T7    =  T_DELTAR
T8    =  T_DELTAF

%PAUSE

% DECLARE THE GENERALIZED COORDINATES
% Q1:  PERPENDICULAR DISTANCE FROM THE N2> AXIS TO THE REAR CONTACT
%      POINT IN THE GROUND PLANE
% Q2:  PERPENDICULAR DISTANCE FROM THE N1> AXIS TO THE REAR CONTACT
%      POINT IN THE GROUND PLANE
% Q3:  REAR FRAME YAW ANGLE
% Q4:  REAR FRAME ROLL ANGLE
% Q5:  REAR FRAME PITCH ANGLE
% Q6:  REAR WHEEL ROTATION ANGLE
% Q7:  MID FRAME STEER ANGLE
% Q8:  FRONT FORK STEER ANGLE
% Q9:  FRONT WHEEL ROTATION ANGLE

VARIABLES Q{9}'

%---------------------------------------------------------------------%
%         GENERALIZED SPEEDS
%---------------------------------------------------------------------%

MOTIONVARIABLES' U{9}'

%---------------------------------------------------------------------%
%         MASS AND INERTIA PROPERTIES
%---------------------------------------------------------------------%

MASS C=MC,D=MD,E=ME,F=MF,G=MG
INERTIA C,IC11,IC22,IC33,IC12,IC23,IC31
INERTIA D,ID11,ID22,ID33
INERTIA E,IE11,IE22,IE33,IE12,IE23,IE31
INERTIA F,IF11,IF22,IF33,IF12,IF23,IF31
INERTIA G,IG11,IG22,IG33

%---------------------------------------------------------------------%
%         ANGULAR RELATIONSHIPS                                       %
%---------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BEGIN TEST CODE SECTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Q1 = 0
%Q2 = 0
%Q3 = 0
%Q4 = 0
%Q5 = 0
%Q6 = 0
%Q7 = 0
%Q8 = 0
%Q9 = 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END TEST CODE SECTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FRAME YAW

SIMPROT(N,A,3,Q3)

% FRAME ROLL

SIMPROT(A,B,1,Q4)

% FRAME PITCH

SIMPROT(B,C,2,Q5)

% REAR WHEEL ROTATION

SIMPROT(C,D,2,Q6)

% MID FRAME ROTATION ANGLE

SIMPROT(C,H,3,Q7)

% MID FRAME SEAT TUBE ANGLE

SIMPROT(H,E,2,LAMBDAR)

% HEAD TUBE ANGLE

SIMPROT(E,I,2,LAMBDAF)

% STEERING ANGLE

SIMPROT(I,F,3,Q8)

% FRONT WHEEL ROTATION

SIMPROT(F,G,2,Q9)

%---------------------------------------------------------------------%
%         POSITION VECTORS
%---------------------------------------------------------------------%

% LOCATE THE CENTER OF MASS FOR EACH BODY

P_NO_DO>=Q1*N1>+Q2*N2>-RR*B3>

P_DO_CO>=L1*C1>+L2*C3>

P_DO_EO>=D1*C1>-D2*C3>+L3*E1>+L4*E3>

P_DO_FO>=D1*C1>-D2*C3>+D3*E1>+L5*F1>+L6*F3>

P_DO_GO>=D1*C1>-D2*C3>+D3*E1>+D4*F3>+D5*F1>

% LOCATE THE GROUND CONTACT POINTS

P_DO_DN>=RR*B3>

P_DN_ND>=0>

P_GO_GN>=RF*UNITVEC(N3>-DOT(F2>,N3>)*F2>)

P_GN_NG>=0>


%---------------------------------------------------------------------%
%         DEFINE THE PITCH CONFIGURATION CONSTRAINT
%---------------------------------------------------------------------%

% SET THE N3> COMPONENT OF P_ND_NG> EQUAL TO ZERO
PZERO=DOT(P_ND_NG>,N3>)

%---------------------------------------------------------------------%
%         DEFINE THE GENERALIZED SPEEDS
%---------------------------------------------------------------------%

Q1'=U1
Q2'=U2
Q3'=U3
Q4'=U4
Q5'=U5
Q6'=U6
Q7'=U7
Q8'=U8
Q9'=U9

%---------------------------------------------------------------------%
%         ANGULAR VELOCITIES
%---------------------------------------------------------------------%

ANGVEL(N,A)
ANGVEL(N,B)
ANGVEL(N,C)
ANGVEL(N,D)
ANGVEL(N,E)
ANGVEL(N,F)
ANGVEL(N,G)

%---------------------------------------------------------------------%
%         VELOCITIES
%---------------------------------------------------------------------%

V_CO_N>=DT(P_NO_CO>,N)
V_DO_N>=DT(P_NO_DO>,N)
V_EO_N>=DT(P_NO_EO>,N)
V_FO_N>=DT(P_NO_FO>,N)
V_GO_N>=DT(P_NO_GO>,N)

V2PTS(N,D,DO,DN)
V2PTS(N,G,GO,GN)

%---------------------------------------------------------------------%
%         MOTION CONSTRAINTS
%---------------------------------------------------------------------%

% DUE TO THE ASSUMPTIONS OF NO SIDE SLIP AND NO SLIP ROLLING THE
% VELOCITIES OF THE FRONT AND REAR WHEEL CONTACT POINTS, DN AND GN,
% CANNOT HAVE COMPONENTS OF VELOCITY IN THE GROUND PLANE

DEPENDENT[1]=DOT(V_DN_N>,N1>)
DEPENDENT[2]=DOT(V_DN_N>,N2>)
DEPENDENT[3]=DOT(V_GN_N>,N1>)
DEPENDENT[4]=DOT(V_GN_N>,N2>)
DEPENDENT[5]=DT(PZERO)

% THE REAR WHEEL ANGULAR SPEED, U6, THE ROLL RATE, U4, THE REAR
% STEERING, U7, AND THE FRONT STEERING RATE, U8, ARE TAKEN TO BE THE
% INDEPENDENT GENERALIZED SPEEDS

CONSTRAIN(DEPENDENT[U1,U2,U3,U5,U9])

%---------------------------------------------------------------------%
%         ANGULAR ACCELERATIONS
%---------------------------------------------------------------------%

ALF_C_N>=DT(W_C_N>,N)

ALF_D_N>=DT(W_D_N>,N)

ALF_E_N>=DT(W_E_N>,N)

ALF_F_N>=DT(W_F_N>,N)

ALF_G_N>=DT(W_G_N>,N)

%---------------------------------------------------------------------%
%         ACCELERATIONS
%---------------------------------------------------------------------%

A_CO_N>=DT(V_CO_N>,N)
A_DO_N>=DT(V_DO_N>,N)
A_EO_N>=DT(V_EO_N>,N)  
A_FO_N>=DT(V_FO_N>,N)  
A_GO_N>=DT(V_GO_N>,N)

%---------------------------------------------------------------------%
%         FORCES AND TORQUES
%---------------------------------------------------------------------%

GRAVITY(G*N3>,C,D,E,F,G)
TORQUE(A/B,T4*A1>) % LEAN TORQUE
TORQUE(C/D,T6*D2>) % REAR WHEEL TORQUE
TORQUE(C/E,T7*C3>) % REAR STEER TORQUE
TORQUE(E/F,T8*F3>) % FRONT STEER TORQUE

%---------------------------------------------------------------------%
%         EQUATIONS OF MOTION
%---------------------------------------------------------------------%

ZERO=FR()+FRSTAR()
SOLVE(ZERO,U4',U6',U7',U8')

%---------------------------------------------------------------------%
%         SIMULATION
%---------------------------------------------------------------------%

UNITSYSTEM  KG,M,S

INPUT W=1.02 M,CF=0.08 M,CR=0.08 M,H=0.8 M,G=9.81 KG

INPUT LAMBDAF=3.141592653589793/10 RAD

INPUT LAMBDAR=3.141592653589793/10 RAD
  
INPUT RR=0.3 M,M_R=2 KG,IRXX=0.0603 KG*M^2,IRYY=0.12 KG*M^2

INPUT XB=0.3 M,ZB=-0.9 M,M_B=85 KG

INPUT IBXX=9.2 KG*M^2,IBYY=11 KG*M^2,IBZZ=2.8 KG*M^2,IBXZ=2.4 KG*M^2

INPUT XH=0.9 M,ZH=-0.7 M,M_H=4 KG

INPUT IHXX=0.05892 KG*M^2,IHYY=0.06 KG*M^2,IHZZ=0.00708 KG*M^2

INPUT IHXZ=-0.00756 KG*M^2

INPUT RF=0.35 M,M_F=3 KG,IFXX=0.1405 KG*M^2,IFYY=0.28 KG*M^2

INPUT XE=0.59 M,ZE=-0.8 M,M_E=3 KG

INPUT IEXX=0.0053 KG*M^2,IEYY=0.1 KG*M^2,IEZZ=0.1 KG*M^2

INPUT IEXZ=0.0 KG*M^2

CODE DYNAMICS() SB_01.M

%---------------------------------------------------------------------%
%         SAVE OUTPUT
%---------------------------------------------------------------------%

SAVE SB_01.ALL
MODULE DyDtMod
! Pleiades problem

USE GlobalCommonMod 

IMPLICIT NONE 

CONTAINS


SUBROUTINE dev(tau,y,dy,test)
  ! Please modify this file and use your own y' = f(t,y)
  ! tau    : current time
  ! y(n) : current y (say, position of a particle in a binary system)
  ! dy(n) : dy/dtau 
  ! test: logical, set it to .True. if you would like to force the solver to stop right now
  !       you may use it to detect Nan values or bad values

  REAL(KIND=8), INTENT(IN) :: tau
  REAL(KIND=8), DIMENSION(:), INTENT(IN) :: y
  REAL(KIND=8), DIMENSION(SIZE(y)), INTENT(OUT) :: dy
  LOGICAL, INTENT(OUT) :: test

  REAL(KIND=8) ::                                                      &
      t,r,theta,phi,V0,V1,V2,V3,Lambda,DlambdaDtau,                    &
      Vdot_0,Vdot_1,Vdot_2,Vdot_3,M,Sg,Dl
  REAL(KIND=8) ::                                                   &
      G0_00,G0_01,G0_02,G0_03,G0_11,G0_12,G0_13,G0_22,G0_23,G0_33,      &
      G1_00,G1_01,G1_02,G1_03,G1_11,G1_12,G1_13,G1_22,G1_23,G1_33,      &
      G2_00,G2_01,G2_02,G2_03,G2_11,G2_12,G2_13,G2_22,G2_23,G2_33,      &
      G3_00,G3_01,G3_02,G3_03,G3_11,G3_12,G3_13,G3_22,G3_23,G3_33,      &
      g00,g01,g02,g03,g11,g12,g13,g22,g23,g33,                          &
      g10,g20,g30,g21,g31,g32,                                          &
      h00,h01,h02,h03,h11,h12,h13,h22,h23,h33,                          &
      h10,h20,h30,h21,h31,h32
  REAL(KIND=8) :: SinTheta, CosTheta
  
  ! not used 
  test = .False.

  M = 1.D0
  
  ! Read in t, r, theta, phi, P0, P1, P2, P3, S0, S1, S2, S3
  t = y(1)
  r = y(2)
  theta = y(3)
  phi = y(4)
  V0 = y(5)
  V1 = y(6)
  V2 = y(7)
  V3 = y(8)
  
  ! global 
  SinTheta = sin(theta)
  CosTheta = cos(theta)
  
  ! Metric
  
  Sg = r**2+a**2*CosTheta**2
  Dl = r**2 + a**2 - 2.D0*M*r
  h00 = -((r**2+a**2)+2.D0*M*r*a**2*SinTheta**2/Sg)/Dl 
  h01 = 0.D0
  h02 = 0.D0
  h03 = -2.D0*M*r*a/(Sg*Dl)
  h11 = Dl/Sg
  h12 = 0.D0
  h13 = 0.D0
  h22 = 1.D0/Sg
  h23 = 0.D0
  h33 = (1.D0-2.D0*M*r/Sg)/(Dl*SinTheta**2)
  
  h10 = h01
  h20 = h02
  h30 = h03
  h21 = h12
  h31 = h13
  h32 = h23
  
  g00 = -(1.D0 - 2.D0*M*r/Sg)
  g01 = 0.D0
  g02 = 0.D0
  g03 = -2.D0*M*r*a*SinTheta**2/Sg
  g11 = Sg/Dl
  g12 = 0.D0
  g13 = 0.D0
  g22 = Sg
  g23 = 0.D0
  g33 = (r**2 + a**2 + 2.D0*M*r*a**2*SinTheta**2/Sg)*SinTheta**2
  
  g10 = g01
  g20 = g02
  g30 = g03
  g21 = g12
  g31 = g13
  g32 = g23
  
  !  Levi-Civita connection.
  
  G0_00 = 0.D0
  G0_01 = M*(r**2+a**2)*(r**2-a**2*CosTheta**2)/(Sg**2*Dl)
  G0_02 = -2.D0*M*r*a**2*SinTheta*CosTheta/Sg**2
  G0_03 = 0.D0
  G0_11 = 0.D0
  G0_12 = 0.D0
  G0_13 = -M*a*((3.D0*r**2-a**2)*(r**2+a**2)-a**2*(r**2-a**2)*SinTheta**2)        &
  *SinTheta**2/(Sg**2*Dl)
  G0_22 = 0.D0
  G0_23 = 2.D0*M*r*a**3*SinTheta**3*CosTheta/Sg**2
  G0_33 = 0.D0
  
  G1_00 = M*Dl*(r**2-a**2*CosTheta**2)/Sg**3
  G1_01 = 0.D0
  G1_02 = 0.D0
  G1_03 = -M*a*Dl*(r**2-a**2*CosTheta**2)*SinTheta**2/Sg**3
  G1_11 = (-M*(r**2-a**2)+a**2*SinTheta**2*(r-M))/(Sg*Dl)
  G1_12 = -a**2*SinTheta*CosTheta/Sg
  G1_13 = 0.D0
  G1_22 = -r*Dl/Sg
  G1_23 = 0.D0
  G1_33 = -(r*(a**2+r**2)**2-a**2*SinTheta**2*(-a**2*SinTheta**2*(r-M)            &
  -M*a**2+r*(2.D0*(r**2+a**2)+M*r)))*Dl*SinTheta**2/Sg**3
  
  G2_00 = -2.D0*M*r*a**2*SinTheta*CosTheta/Sg**3
  G2_01 = 0.D0
  G2_02 = 0.D0
  G2_03 = 2.D0*M*r*a*(r**2+a**2)*SinTheta*CosTheta/Sg**3
  G2_11 = a**2*SinTheta*CosTheta/(Sg*Dl)
  G2_12 = r/Sg
  G2_13 = 0.D0
  G2_22 = -a**2*SinTheta*CosTheta/Sg
  G2_23 = 0.D0
  G2_33 = -((r**2+a**2)**3-a**2*Dl*(2.D0*(r**2+a**2)-a**2*SinTheta**2)*SinTheta**2) &
  *SinTheta*CosTheta/Sg**3
  
  G3_00 = 0.D0
  G3_01 = M*a*(r**2-a**2*CosTheta**2)/(Sg**2*Dl)
  G3_02 = -2.D0*M*r*a*(CosTheta/SinTheta)/Sg**2
  G3_03 = 0.D0
  G3_11 = 0.D0
  G3_12 = 0.D0
  G3_13 = (r*Dl*(r**2+a**2)+a**4*(r-M)*SinTheta**4-a**2*(a**2+r**2)*(2.D0*r-M)*SinTheta**2) &
  /(Sg**2*Dl)
  G3_22 = 0.D0
  G3_23 = ((r**2+a**2)**2+SinTheta**4*a**4-2.D0*a**2*(r**2-M*r+a**2)*SinTheta**2)*(CosTheta/SinTheta)/Sg**2
  G3_33 = 0.D0
  
  !  Force due to geodesic.
  
  Vdot_0 = -G0_00*V0*V0-G0_01*V0*V1-G0_02*V0*V2-G0_03*V0*V3-G0_01*V1*V0-G0_11*V1*V1-G0_12*V1 &
  *V2-G0_13*V1*V3-G0_02*V2*V0-G0_12*V2*V1-G0_22*V2*V2-G0_23*V2*V3-G0_03*V3*V0-G0_13*V3*V1     &
  -G0_23*V3*V2-G0_33*V3*V3
  
  Vdot_1 = -G1_00*V0*V0-G1_01*V0*V1-G1_02*V0*V2-G1_03*V0*V3-G1_01*V1*V0-G1_11*V1*V1-G1_12*V1 &
  *V2-G1_13*V1*V3-G1_02*V2*V0-G1_12*V2*V1-G1_22*V2*V2-G1_23*V2*V3-G1_03*V3*V0-G1_13*V3*V1     &
  -G1_23*V3*V2-G1_33*V3*V3
  
  Vdot_2 = -G2_00*V0*V0-G2_01*V0*V1-G2_02*V0*V2-G2_03*V0*V3-G2_01*V1*V0-G2_11*V1*V1-G2_12*V1 &
  *V2-G2_13*V1*V3-G2_02*V2*V0-G2_12*V2*V1-G2_22*V2*V2-G2_23*V2*V3-G2_03*V3*V0-G2_13*V3*V1     &
  -G2_23*V3*V2-G2_33*V3*V3
  
  Vdot_3 = -G3_00*V0*V0-G3_01*V0*V1-G3_02*V0*V2-G3_03*V0*V3-G3_01*V1*V0-G3_11*V1*V1-G3_12*V1 &
  *V2-G3_13*V1*V3-G3_02*V2*V0-G3_12*V2*V1-G3_22*V2*V2-G3_23*V2*V3-G3_03*V3*V0-G3_13*V3*V1     &
  -G3_23*V3*V2-G3_33*V3*V3
  
  ! assign values 
  dy(1) = V0
  dy(2) = V1
  dy(3) = V2 !velocity in theta direction
  dy(4) = V3
  dy(5) = Vdot_0
  dy(6) = Vdot_1
  dy(7) = Vdot_2
  dy(8) = Vdot_3

  Evaluated = Evaluated + 1
  RETURN 
END SUBROUTINE
 

END MODULE

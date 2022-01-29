#define  PHYSICS                        HD
#define  DIMENSIONS                     3
#define  GEOMETRY                       CARTESIAN
#define  BODY_FORCE                     POTENTIAL
#define  COOLING                        NO
#define  RECONSTRUCTION                 PARABOLIC
#define  TIME_STEPPING                  RK3
#define  NTRACER                        1
#define  PARTICLES                      NO
#define  USER_DEF_PARAMETERS            7

/* -- physics dependent declarations -- */

#define  DUST_FLUID                     NO
#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 NO
#define  THERMAL_CONDUCTION             NO
#define  VISCOSITY                      NO
#define  ROTATING_FRAME                 NO

/* -- user-defined parameters (labels) -- */

#define  T_INIT                         0
#define  RHOC                           1
#define  DN_RATIO                       2
#define  PR_RATIO                       3
#define  VELJ                           4
#define  RJET                           5
#define  JET_TIME                       6

/* [Beg] user-defined constants (do not change this line) */

#define  UNIT_DENSITY                   (0.01*CONST_mp)
#define  UNIT_LENGTH                    (1000*CONST_pc)
#define  UNIT_VELOCITY                  1.e5
#define  SHOCK_FLATTENING               YES
#define  SHOW_TIME_STEPS                YES
#define  INTERNAL_BOUNDARY              YES
#define  STOP_IF                        YES

/* [End] user-defined constants (do not change this line) */

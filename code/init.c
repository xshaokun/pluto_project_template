/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Contains basic functions for problem initialization.

  The init.c file collects most of the user-supplied functions useful
  for problem configuration.
  It is automatically searched for by the makefile.

  \author A. Mignone (mignone@ph.unito.it)
  \date   March 5, 2017
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void Init (double *v, double x1, double x2, double x3)
/*!
 * The Init() function can be used to assign initial conditions as
 * as a function of spatial position.
 *
 * \param [out] v   a pointer to a vector of primitive variables
 * \param [in] x1   coordinate point in the 1st dimension
 * \param [in] x2   coordinate point in the 2nd dimension
 * \param [in] x3   coordinate point in the 3rd dimension
 *
 * The meaning of x1, x2 and x3 depends on the geometry:
 * \f[ \begin{array}{cccl}
 *    x_1  & x_2    & x_3  & \mathrm{Geometry}    \\ \noalign{\medskip}
 *     \hline
 *    x    &   y    &  z   & \mathrm{Cartesian}   \\ \noalign{\medskip}
 *    R    &   z    &  -   & \mathrm{cylindrical} \\ \noalign{\medskip}
 *    R    & \phi   &  z   & \mathrm{polar}       \\ \noalign{\medskip}
 *    r    & \theta & \phi & \mathrm{spherical}
 *    \end{array}
 *  \f]
 *
 * Variable names are accessed by means of an index v[nv], where
 * nv = RHO is density, nv = PRS is pressure, nv = (VX1, VX2, VX3) are
 * the three components of velocity, and so forth.
 *
 *********************************************************************** */
{
  double mu,vinit[NVAR];
  mu = MeanMolecularWeight(vinit);

  v[RHO] = 0.0;
  v[VX1] = 0.0;
  v[VX2] = 0.0;
  v[VX3] = 0.0;
  #if HAVE_ENERGY
  v[PRS] = 0.0;
  #endif
  v[TRC] = 0.0;

  #if PHYSICS == MHD || PHYSICS == RMHD
  v[BX1] = 0.0;
  v[BX2] = 0.0;
  v[BX3] = 0.0;

  v[AX1] = 0.0;
  v[AX2] = 0.0;
  v[AX3] = 0.0;
  #endif
}

/* ********************************************************************* */
void InitDomain (Data *d, Grid *grid)
/*!
 * Assign initial condition by looping over the computational domain.
 * Called after the usual Init() function to assign initial conditions
 * on primitive variables.
 * Value assigned here will overwrite those prescribed during Init().
 *
 * - Details:
 *   o Isothermal halo gas in hydrostatic equilibrium.
 *   o Given the value of central cell, integral along k direction first, then
 *     along j direction (if included), and finally i direction.
 *
 *********************************************************************** */
{
  int i,j,k,dim,nv;

/* ------------------------------------------------------------------
   Locate the procossor having origin
   ------------------------------------------------------------------ */
  int count=0;
  int rc_orig[3]; // rank coordinate of origin block
  int *rc = grid->rank_coord;  // parallel coordinate in a Cartesian topology
  int orig[3]; // local index of the origin

  DIM_LOOP(dim)
    if (grid->xbeg[dim] * grid->xend[dim] <= 0.0)
      count++;

  nv = 0;
  if (count == 3) { // locate the procossor having origin
    nv = prank;
    DIM_LOOP(dim) {
      rc_orig[dim] = rc[dim];
      for (i=grid->lbeg[dim]; i<=grid->lend[dim]; i++)
        if (grid->x[dim][i-1]*grid->x[dim][i] < 0.0) {
          orig[dim] = i;
          break;
        }
    }
  }

  int rank_orig=0; // store the id of procossor having the origin
  MPI_Allreduce (&nv, &rank_orig, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  MPI_Bcast (rc_orig, 3, MPI_INT, rank_orig, MPI_COMM_WORLD);
  MPI_Bcast (orig, 3, MPI_INT, rank_orig, MPI_COMM_WORLD);

/* ------------------------------------------------------------------
   Define start and end indices for loop
   ------------------------------------------------------------------ */

  int ip, times, beg[3], pdir[3]; // index for processor loop
  double mu,vinit[NVAR];
  double phi,phi_old,dphi,alpha,rho_old;
  int par_dim[3] = {0, 0, 0};
  int *np = grid->nproc;   // number of processors in each dim
  double *x1 = grid->xgc[IDIR];
  double *x2 = grid->xgc[JDIR];
  double *x3 = grid->xgc[KDIR];

  mu = MeanMolecularWeight(vinit);

  DIM_EXPAND(par_dim[0] = grid->nproc[IDIR] > 1;  ,
             par_dim[1] = grid->nproc[JDIR] > 1;  ,
             par_dim[2] = grid->nproc[KDIR] > 1;)

  if (prank == rank_orig) {
    k = orig[KDIR];
    j = orig[JDIR];
    i = orig[IDIR];
    d->Vc[RHO][k][j][i] = g_inputParam[RHOC];
    d->Vc[PRS][k][j][i] = d->Vc[RHO][k][j][i] * g_inputParam[T_INIT]/(KELVIN*mu);
  }

  DIM_LOOP(dim)
  if (rc[dim]>rc_orig[dim]) {
    beg[dim] = grid->lbeg[dim];
    pdir[dim] = 1;
  } else if (rc[dim]<rc_orig[dim]) {
    beg[dim] = grid->lend[dim];
    pdir[dim] = -1;
  } else {
    beg[dim] = orig[dim];
    pdir[dim] = 1;
  }

  #if DEBUG == TRUE
  print("prank=%d, rank_orig=%d, rc_orig=", prank, rank_orig);
  DIM_LOOP(dim)
    print("%d ", rc_orig[dim]);
  print("rc=");
  DIM_LOOP(dim)
    print("%d ", rc[dim]);
  print("beg=");
  DIM_LOOP(dim)
    print("%d ", beg[dim]);
  print("orig=");
  DIM_LOOP(dim)
    print("%d ", orig[dim]);
  print("\n");
  #endif

  for (dim=2; dim>=0; dim--) {

/* ------------------------------------------------------------------
   Then integrate from that procossor
   ------------------------------------------------------------------ */

    if (rc[IDIR] == rc_orig[IDIR] || dim < 0)
    if (rc[JDIR] == rc_orig[JDIR] || dim < 1)
    if (rc[KDIR] == rc_orig[KDIR] || dim < 2) {
      switch (dim) {
        case KDIR:
          // integrate in positive direction
          j = orig[JDIR];
          i = orig[IDIR];
          for (k=orig[KDIR]+1,phi_old=BodyForcePotential(x1[i],x2[j],x3[k-1]); k<=KEND; k++) {
            phi = BodyForcePotential(x1[i],x2[j],x3[k]);
            dphi = phi - phi_old;
            alpha = dphi/(2.*g_inputParam[T_INIT]/KELVIN/mu);
            d->Vc[RHO][k][j][i] = d->Vc[RHO][k-1][j][i] * (1.-alpha)/(1.+alpha);
            d->Vc[PRS][k][j][i] = d->Vc[RHO][k][j][i] * g_inputParam[T_INIT]/(KELVIN*mu);
            phi_old = phi;
          }
          // integrate in negative direction
          for (k=orig[KDIR]-1,phi_old=BodyForcePotential(x1[i],x2[j],x3[k+1]); k>=KBEG; k--) {
            phi = BodyForcePotential(x1[i],x2[j],x3[k]);
            dphi = phi - phi_old;
            alpha = dphi/(2.*g_inputParam[T_INIT]/KELVIN/mu);
            d->Vc[RHO][k][j][i] = d->Vc[RHO][k+1][j][i] * (1.-alpha)/(1.+alpha);
            d->Vc[PRS][k][j][i] = d->Vc[RHO][k][j][i] * g_inputParam[T_INIT]/(KELVIN*mu);
            phi_old = phi;
          }
          break;

        case JDIR:
          i = orig[IDIR];
          KDOM_LOOP(k) {
            // integrate in positive direction
            for (j=orig[JDIR]+1,phi_old=BodyForcePotential(x1[i],x2[j-1],x3[k]); j<=JEND; j++) {
              phi = BodyForcePotential(x1[i],x2[j],x3[k]);
              dphi = phi - phi_old;
              alpha = dphi/(2.*g_inputParam[T_INIT]/KELVIN/mu);
              d->Vc[RHO][k][j][i] = d->Vc[RHO][k][j-1][i] * (1.-alpha)/(1.+alpha);
              d->Vc[PRS][k][j][i] = d->Vc[RHO][k][j][i] * g_inputParam[T_INIT]/(KELVIN*mu);
              phi_old = phi;
            }
            // integrate in negative direction
            for (j=orig[JDIR]-1,phi_old=BodyForcePotential(x1[i],x2[j+1],x3[k]); j>=JBEG; j--) {
              phi = BodyForcePotential(x1[i],x2[j],x3[k]);
              dphi = phi - phi_old;
              alpha = dphi/(2.*g_inputParam[T_INIT]/KELVIN/mu);
              d->Vc[RHO][k][j][i] = d->Vc[RHO][k][j+1][i] * (1.-alpha)/(1.+alpha);
              d->Vc[PRS][k][j][i] = d->Vc[RHO][k][j][i] * g_inputParam[T_INIT]/(KELVIN*mu);
              phi_old = phi;
            }
          }
          break;

        case IDIR:
          KDOM_LOOP(k)
          JDOM_LOOP(j) {
            // integrate in positive direction
            for (i=orig[IDIR]+1,phi_old=BodyForcePotential(x1[i-1],x2[j],x3[k]); i<=IEND; i++) {
              phi = BodyForcePotential(x1[i],x2[j],x3[k]);
              dphi = phi - phi_old;
              alpha = dphi/(2.*g_inputParam[T_INIT]/KELVIN/mu);
              d->Vc[RHO][k][j][i] = d->Vc[RHO][k][j][i-1] * (1.-alpha)/(1.+alpha);
              d->Vc[PRS][k][j][i] = d->Vc[RHO][k][j][i] * g_inputParam[T_INIT]/(KELVIN*mu);
              phi_old = phi;
            }
            // integrate in negative direction
            for (i=orig[IDIR]-1,phi_old=BodyForcePotential(x1[i+1],x2[j],x3[k]); i>=IBEG; i--) {
              phi = BodyForcePotential(x1[i],x2[j],x3[k]);
              dphi = phi - phi_old;
              alpha = dphi/(2.*g_inputParam[T_INIT]/KELVIN/mu);
              d->Vc[RHO][k][j][i] = d->Vc[RHO][k][j][i+1] * (1.-alpha)/(1.+alpha);
              d->Vc[PRS][k][j][i] = d->Vc[RHO][k][j][i] * g_inputParam[T_INIT]/(KELVIN*mu);
              phi_old = phi;
            }
          }
          break;
      }
    }
    NVAR_LOOP(nv) AL_Exchange_dim ((char *)d->Vc[nv][0][0], par_dim, SZ);

/* ------------------------------------------------------------------
   Integrate along outer region in certain direction
   ------------------------------------------------------------------ */

    ip = rc_orig[dim]+pdir[dim];
    times = MAX(rc_orig[dim], grid->nproc[dim]-rc_orig[dim]-1);
    for (; times>0; ip+=pdir[dim],times--) {

      if (rc[dim] == ip)
      if (rc[IDIR] == rc_orig[IDIR] || dim < 1)
      if (rc[JDIR] == rc_orig[JDIR] || dim < 2) {
        switch (dim) {
          case KDIR:
            j = orig[JDIR];
            i = orig[IDIR];
            for (k = beg[KDIR],phi_old=BodyForcePotential(x1[i],x2[j],x3[k-pdir[KDIR]]); k>=KBEG && k<=KEND; k+=pdir[KDIR]) {
              phi = BodyForcePotential(x1[i],x2[j],x3[k]);
              dphi = phi - phi_old;
              alpha = dphi/(2*g_inputParam[T_INIT]/KELVIN/mu);
              d->Vc[RHO][k][j][i] = d->Vc[RHO][k-pdir[KDIR]][j][i] * (1.-alpha)/(1.+alpha);
              d->Vc[PRS][k][j][i] = d->Vc[RHO][k][j][i] * g_inputParam[T_INIT]/(KELVIN*mu);
              phi_old = phi;
            }
            break;

          case JDIR:
            i=orig[IDIR];
            KDOM_LOOP(k) {
              for (j = beg[JDIR],phi_old=BodyForcePotential(x1[i],x2[j-pdir[JDIR]],x3[k]); j>=JBEG && j<=JEND; j+=pdir[JDIR]) {
                phi = BodyForcePotential(x1[i],x2[j],x3[k]);
                dphi = phi - phi_old;
                alpha = dphi/(2*g_inputParam[T_INIT]/KELVIN/mu);
                d->Vc[RHO][k][j][i] = d->Vc[RHO][k][j-pdir[JDIR]][i] * (1.-alpha)/(1.+alpha);
                d->Vc[PRS][k][j][i] = d->Vc[RHO][k][j][i] * g_inputParam[T_INIT]/(KELVIN*mu);
                phi_old = phi;
              }
            }
            break;

          case IDIR:
            KDOM_LOOP(k)
            JDOM_LOOP(j)
            for (i = beg[IDIR],phi_old=BodyForcePotential(x1[i-pdir[IDIR]],x2[j],x3[k]); i>=IBEG && i<=IEND; i+=pdir[IDIR]) {
              phi = BodyForcePotential(x1[i],x2[j],x3[k]);
              dphi = phi - phi_old;
              alpha = dphi/(2*g_inputParam[T_INIT]/KELVIN/mu);
              d->Vc[RHO][k][j][i] = d->Vc[RHO][k][j][i-pdir[IDIR]] * (1.-alpha)/(1.+alpha);
              d->Vc[PRS][k][j][i] = d->Vc[RHO][k][j][i] * g_inputParam[T_INIT]/(KELVIN*mu);
              phi_old = phi;
            }
            break;
        }
      }
      NVAR_LOOP(nv) AL_Exchange_dim ((char *)d->Vc[nv][0][0], par_dim, SZ);
    }
  }
}

/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/*!
 *  Perform runtime data analysis.
 *
 * \param [in] d the PLUTO Data structure
 * \param [in] grid   pointer to array of Grid structures
 *
 *********************************************************************** */
{
  int    i, j, k;
  double ***dV = grid->dV;
  double vol, scrh, vr2, vz2;
  double Mtot, Ekin, Eth, Epot, Etot, Pr, Pz;
  double *rl = grid->xl[IDIR];
  double *rr = grid->xr[IDIR];
  double *dz = grid->dx[KDIR];
  double *x1 = grid->xgc[IDIR];
  double *x2 = grid->xgc[JDIR];
  double *x3 = grid->xgc[KDIR];

  Mtot = Ekin = Eth = Epot = Etot = Pz= Pr = 0.0;
  DOM_LOOP(k,j,i) {
    if (INCLUDE_JDIR==YES)
      vol = dV[k][j][i];
    else
      vol = 2*CONST_PI*dV[k][j][i];

    vr2 = d->Vc[VX1][k][j][i]*d->Vc[VX1][k][j][i];
    vz2 = d->Vc[VX3][k][j][i]*d->Vc[VX3][k][j][i];

    Mtot += d->Vc[RHO][k][j][i]*vol;

    scrh = 0.5*d->Vc[RHO][k][j][i]*(vr2 + vz2);
    Ekin += scrh*vol;

    scrh = d->Vc[PRS][k][j][i]/(g_gamma - 1.0);
    Eth += scrh*vol;

    scrh = d->Vc[RHO][k][j][i]*BodyForcePotential(x1[i],x2[j],x3[k]);
    Epot += scrh*vol;
  }

  Etot = Ekin + Eth + Epot;

  MPI_Allreduce(&Mtot, &vol, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  Mtot = vol;

  MPI_Allreduce(&Ekin, &vol, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  Ekin = vol;

  MPI_Allreduce(&Eth, &vol, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  Eth = vol;

  MPI_Allreduce(&Epot, &vol, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  Epot = vol;

  MPI_Allreduce(&Etot, &vol, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  Etot = vol;

  MPI_Barrier(MPI_COMM_WORLD);

/* -- save the values computed above to hist.out -- */
  if (prank == 0){
    char fname[512];
    static double tpos = -1.0;
    FILE *fp;

    sprintf(fname, "%s/hist.out", RuntimeGet()->output_dir);
    if (g_stepNumber == 0) {
      fp = fopen(fname,"w");
      fprintf(fp,"%7s  %12s  %12s  %12s  %12s  %12s  %12s  %12s  %12s\n", "t", "dt", "Mtot", "Ekin", "Eth", "Epot", "Etot", "Pr", "Pz");
    } else {
      if (tpos < 0.0) {     // for restarting the simulation
        char sline[512];
        fp = fopen(fname,"r");
        while (fgets(sline,512,fp)) {}
        sscanf(sline, "%lf\n",&tpos);
        fclose(fp);
      }
      fp = fopen(fname,"a");
    }
    if (g_time > tpos) {
      fprintf(fp, "%12.6e  %12.6e  %12.6e  %12.6e  %12.6e  %12.6e  %12.6e  %12.6e  %12.6e \n", g_time, g_dt, Mtot, Ekin, Eth, Epot, Etot, Pr, Pz);
    }
    fclose(fp);
  }
}

#if PHYSICS == MHD
/* ********************************************************************* */
void BackgroundField (double x1, double x2, double x3, double *B0)
/*!
 * Define the component of a static, curl-free background
 * magnetic field.
 *
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 * \param [out] B0 array containing the vector componens of the background
 *                 magnetic field
 *********************************************************************** */
{
   B0[0] = 0.0;
   B0[1] = 0.0;
   B0[2] = 0.0;
}
#endif

/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid)
/*!
 *  Assign user-defined boundary conditions.
 *
 * \param [in,out] d  pointer to the PLUTO data structure containing
 *                    cell-centered primitive quantities (d->Vc) and
 *                    staggered magnetic fields (d->Vs, when used) to
 *                    be filled.
 * \param [in] box    pointer to a RBox structure containing the lower
 *                    and upper indices of the ghost zone-centers/nodes
 *                    or edges at which data values should be assigned.
 * \param [in] side   specifies the boundary side where ghost zones need
 *                    to be filled. It can assume the following
 *                    pre-definite values: X1_BEG, X1_END,
 *                                         X2_BEG, X2_END,
 *                                         X3_BEG, X3_END.
 *                    The special value side == 0 is used to control
 *                    a region inside the computational domain.
 * \param [in] grid  pointer to an array of Grid structures.
 *
 * - Details:
 *   o Jet nozzle is implemented in ghost cells, within jet radius.
 *   o Jet is parameterized by density raio(of jet density to halo densty at center)
 *     and pressure/energy ratio(of jet pressure/energy to halo pressure/energy at center)
 *********************************************************************** */
{
  int   i, j, k, nv;
  static int first_call = 1;
  double *x1 = grid->x[IDIR];
  double *x2 = grid->x[JDIR];
  double *x3 = grid->x[KDIR];
  double R, dx3;
  double mu, vj[NVAR], momtj, energj;
  double phi, phi_old, dphi, alpha;

  mu = MeanMolecularWeight(vj);

  vj[RHO] = g_inputParam[DN_RATIO]*g_inputParam[RHOC];
  vj[VX3] = g_inputParam[VELJ];
  vj[PRS] = g_inputParam[PR_RATIO] * g_inputParam[RHOC] * g_inputParam[T_INIT]/(KELVIN*mu);
  vj[TRC] = 1.0;

  if (first_call) {
    /* ----------------------------------------------------
     Calculate jet power, and corresponding accretion rate
     ------------------------------------------------------ */
    double flux;
    double P_kin, P_th, P_jet, P0, Mdot, Mdot0;
    double cs_a, mach;

    momtj = vj[RHO] * vj[VX3];
    energj = vj[RHO] * (vj[PRS]/(g_gamma-1) + 0.5*vj[VX3]*vj[VX3]);
    P0 = UNIT_DENSITY * UNIT_LENGTH * UNIT_LENGTH * UNIT_VELOCITY * UNIT_VELOCITY * UNIT_VELOCITY;  // in erg/s
    Mdot0 = UNIT_DENSITY * UNIT_LENGTH * UNIT_LENGTH * UNIT_VELOCITY / CONST_Msun * 3.16e7;   // in Msun/yr
    flux = CONST_PI * g_inputParam[RJET] * g_inputParam[RJET] * vj[VX3];
    Mdot = vj[RHO] * flux;
    P_kin = 0.5*Mdot*vj[VX3]*vj[VX3];
    P_th = vj[PRS]/(g_gamma-1.0) * flux;
    P_jet = P_kin + P_th;

    Mdot *= Mdot0;
    P_kin *= P0;
    P_th *= P0;
    P_jet *= P0;

    cs_a = sqrt(g_gamma*d->Vc[PRS][KBEG][JBEG][IBEG]/d->Vc[RHO][KBEG][JBEG][IBEG]);  // only if isothermal
    mach = vj[VX3] / cs_a * sqrt(g_inputParam[DN_RATIO]/g_inputParam[PR_RATIO]);

    if (prank==0) {
      print ("> Jet Parameters:\n\n");
      print ("  [Density]:          %8.3e (g/cm**3)\n",vj[RHO]*UNIT_DENSITY);
      print ("  [Energy Density]:   %8.3e (erg/cm**3)\n",vj[PRS]*UNIT_DENSITY*UNIT_VELOCITY*UNIT_VELOCITY);
      print ("  [Mass Flux]:        %8.3e (Msun/yr)\n",Mdot*Mdot0);
      print ("  [Luminisity]:       %8.3e (erg/s)\n",P_jet);
      print ("  [Kinetic Power]:    %8.3e (erg/s)\n",P_kin);
      print ("  [Thermal Power]:    %8.3e (erg/s)\n",P_th);
      print ("  [Kinetic Fraction]: %8.3f \% \n",P_kin/P_jet*100);
      print ("  [Mach Number]:      %8.3e\n",mach);
      print ("  [Sound Speed]:      %8.3e (km/s)  [ambiant for reference]\n",cs_a);
      print ("  \n");

    first_call = 0;
    }
  }
  /* --------------------------------------------------------- */

  if (side == X1_BEG){
    /* -- Hydrostatic equilibrium outflow boundary condition -- */
    BOX_LOOP(box, k,j,i) {
      phi_old = BodyForcePotential(x1[i+1],x2[j],x3[k]);
      phi = BodyForcePotential(x1[i],x2[j],x3[k]);
      dphi = phi - phi_old;
      alpha = dphi/(2.*g_inputParam[T_INIT]/KELVIN/mu);
      NVAR_LOOP(nv) d->Vc[nv][k][j][i] = d->Vc[nv][k][j][i+1];
      d->Vc[RHO][k][j][i] = d->Vc[RHO][k][j][i+1] * (1.-alpha)/(1.+alpha);
      d->Vc[PRS][k][j][i] = d->Vc[RHO][k][j][i] * g_inputParam[T_INIT]/(KELVIN*mu);
    }
  }

  if (side == X1_END){
    /* -- Hydrostatic equilibrium outflow boundary condition -- */
    BOX_LOOP(box, k,j,i) {
      phi_old = BodyForcePotential(x1[i-1],x2[j],x3[k]);
      phi = BodyForcePotential(x1[i],x2[j],x3[k]);
      dphi = phi - phi_old;
      alpha = dphi/(2.*g_inputParam[T_INIT]/KELVIN/mu);
      NVAR_LOOP(nv) d->Vc[nv][k][j][i] = d->Vc[nv][k][j][i-1];
      d->Vc[RHO][k][j][i] = d->Vc[RHO][k][j][i-1] * (1.-alpha)/(1.+alpha);
      d->Vc[PRS][k][j][i] = d->Vc[RHO][k][j][i] * g_inputParam[T_INIT]/(KELVIN*mu);
    }
  }

  if (side == X2_BEG){
    /* -- Hydrostatic equilibrium outflow boundary condition -- */
    BOX_LOOP(box, k,j,i) {
      phi_old = BodyForcePotential(x1[i],x2[j+1],x3[k]);
      phi = BodyForcePotential(x1[i],x2[j],x3[k]);
      dphi = phi - phi_old;
      alpha = dphi/(2.*g_inputParam[T_INIT]/KELVIN/mu);
      NVAR_LOOP(nv) d->Vc[nv][k][j][i] = d->Vc[nv][k][j+1][i];
      d->Vc[RHO][k][j][i] = d->Vc[RHO][k][j+1][i] * (1.-alpha)/(1.+alpha);
      d->Vc[PRS][k][j][i] = d->Vc[RHO][k][j][i] * g_inputParam[T_INIT]/(KELVIN*mu);
    }
  }

  if (side == X2_END){
    /* -- Hydrostatic equilibrium outflow boundary condition -- */
    BOX_LOOP(box, k,j,i) {
      phi_old = BodyForcePotential(x1[i],x2[j-1],x3[k]);
      phi = BodyForcePotential(x1[i],x2[j],x3[k]);
      dphi = phi - phi_old;
      alpha = dphi/(2.*g_inputParam[T_INIT]/KELVIN/mu);
      NVAR_LOOP(nv) d->Vc[nv][k][j][i] = d->Vc[nv][k][j-1][i];
      d->Vc[RHO][k][j][i] = d->Vc[RHO][k][j-1][i] * (1.-alpha)/(1.+alpha);
      d->Vc[PRS][k][j][i] = d->Vc[RHO][k][j][i] * g_inputParam[T_INIT]/(KELVIN*mu);
    }
  }

  if (side == X3_BEG){  /* -- X3_BEG boundary -- */
    /* -- Hydrostatic equilibrium reflective boundary condition -- */
    BOX_LOOP(box, k,j,i) {
      phi_old = BodyForcePotential(x1[i],x2[j],x3[k+1]);
      phi = BodyForcePotential(x1[i],x2[j],x3[k]);
      dphi = phi - phi_old;
      alpha = dphi/(2.*g_inputParam[T_INIT]/KELVIN/mu);
      NVAR_LOOP(nv) d->Vc[nv][k][j][i] = d->Vc[nv][k+1][j][i];
      d->Vc[RHO][k][j][i] = d->Vc[RHO][k+1][j][i] * (1.-alpha)/(1.+alpha);
      d->Vc[PRS][k][j][i] = d->Vc[RHO][k][j][i] * g_inputParam[T_INIT]/(KELVIN*mu);
    }
  }

  if (side == X3_END){
    /* -- Hydrostatic equilibrium outflow boundary condition -- */
    BOX_LOOP(box, k,j,i) {
      phi_old = BodyForcePotential(x1[i],x2[j],x3[k-1]);
      phi = BodyForcePotential(x1[i],x2[j],x3[k]);
      dphi = phi - phi_old;
      alpha = dphi/(2.*g_inputParam[T_INIT]/KELVIN/mu);
      NVAR_LOOP(nv) d->Vc[nv][k][j][i] = d->Vc[nv][k-1][j][i];
      d->Vc[RHO][k][j][i] = d->Vc[RHO][k-1][j][i] * (1.-alpha)/(1.+alpha);
      d->Vc[PRS][k][j][i] = d->Vc[RHO][k][j][i] * g_inputParam[T_INIT]/(KELVIN*mu);
    }
  }

  if (side == 0) {
    DOM_LOOP(k,j,i) {
      if (g_time <= g_inputParam[JET_TIME]) {
        R = sqrt(x1[i]*x1[i] + x2[j]*x2[j]);
        dx3 = grid->dx[KDIR][KBEG];
        if (R <= g_inputParam[RJET])
        if (fabs(x3[k]) <= 3*dx3 && x3[k]>=0.0) {  /* -- 3 cells thickness -- */
          d->Vc[RHO][k][j][i] = vj[RHO]/cosh(pow(R/g_inputParam[RJET],6));
          d->Vc[VX1][k][j][i] = 0.0;
          d->Vc[VX2][k][j][i] = 0.0;
          d->Vc[VX3][k][j][i] = vj[VX3];
          d->Vc[PRS][k][j][i] = vj[PRS]/cosh(pow(R/g_inputParam[RJET],6));
          d->Vc[TRC][k][j][i] = vj[TRC];
          d->flag[k][j][i] |= FLAG_INTERNAL_BOUNDARY;
        } else if (fabs(x3[k]) <= 3*dx3 && x3[k]<=0.0) {
          d->Vc[RHO][k][j][i] = vj[RHO]/cosh(pow(R/g_inputParam[RJET],6));
          d->Vc[VX1][k][j][i] = 0.0;
          d->Vc[VX2][k][j][i] = 0.0;
          d->Vc[VX3][k][j][i] = -vj[VX3];
          d->Vc[PRS][k][j][i] = vj[PRS]/cosh(pow(R/g_inputParam[RJET],6));
          d->Vc[TRC][k][j][i] = vj[TRC];
          d->flag[k][j][i] |= FLAG_INTERNAL_BOUNDARY;
        }
      }
    }
  }
}

#if BODY_FORCE != NO
/* ********************************************************************* */
void BodyForceVector(double *v, double *g, double x1, double x2, double x3)
/*!
 * Prescribe the acceleration vector as a function of the coordinates
 * and the vector of primitive variables *v.
 *
 * \param [in] v  pointer to a cell-centered vector of primitive
 *                variables
 * \param [out] g acceleration vector
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 *
 *********************************************************************** */
{
  g[IDIR] = 0.0;
  g[JDIR] = 0.0;
  g[KDIR] = 0.0;
}
/* ********************************************************************* */
double BodyForcePotential(double x1, double x2, double x3)
/*!
 * Return the gravitational potential as function of the coordinates.
 *
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 *
 * \return The body force potential \f$ \Phi(x_1,x_2,x_3) \f$.
 *
 * - Details:
 *   o Dark matter halo: NFW model
 *   o Stellar disk + stellar bulge: model used in Guo & Mathews (2012)
 *********************************************************************** */
{
  double phi0,phi,phi_nfw,phi_bulge,phi_disc;
  const double rho_crit=8.6e-30, Delta=200, M200=5.71e12, h=0.6766;
  const double Mbulge=7.9e8, d=0.7;
  const double Mdisc=9.42e10, a=6.5, b=0.26;
  double R, z, r, rs, rd, c200;

  phi0 = UNIT_VELOCITY*UNIT_VELOCITY;
  c200 = 5.26*pow(M200/1.e14/h,-0.1);  //concentration, refer to Neto+2007
  rs = 1/c200*pow(3*M200/(4*CONST_PI*rho_crit*Delta),1./3.);  // characteristic radius

  #if GEOMETRY == CARTESIAN
    R = sqrt(x1*x1+x2*x2);
    z = fabs(x3);
    r = sqrt(x1*x1+x2*x2+x3*x3);
  #elif GEOMETRY == POLAR
    R = x1;
    z = x3;
    r = sqrt(x1*x1+x3*x3);
  #elif GEOMETRY == SPHERICAL
    r = x1;
    R = r*sin(x2);
    z = r*cos(x2);
  #endif

  phi_nfw = -3*M200/c200 * log(r/rs+1)/r / phi0;
  phi_bulge = -Mbulge/(r+d) / phi0;

  rd = sqrt(R*R+pow(a+sqrt(z*z+b*b),2.));
  phi_disc = -Mdisc/rd / phi0;

  phi = CONST_G*(phi_nfw+phi_bulge+phi_disc)*CONST_Msun/(1.e3*CONST_pc);

  return phi;
}
#endif

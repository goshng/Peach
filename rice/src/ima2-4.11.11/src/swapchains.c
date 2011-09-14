/* IMa2  2009-2011  Jody Hey, Rasmus Nielsen and Sang Chul Choi*/

#undef GLOBVARS
#include "imamp.h"

/*********** LOCAL STUFF **********/
static int swap01, swap01d, swap12, swap12d;
static unsigned long swapcount[MAXCHAINS][MAXCHAINS];

static double swapweight (int ci, int cj);


double
swapweight (int ci, int cj)
{
  int li;
  double sumi = 0, sumj = 0, w;
  for (li = 0; li < nloci; li++)
  {
    sumi += C[ci]->G[li].pdg;
    sumj += C[cj]->G[li].pdg;
  }
  sumi += C[ci]->allpcalc.probg;
  sumj += C[cj]->allpcalc.probg;
  w = exp ((beta[ci] - beta[cj]) * (sumj - sumi));
  return (w);
}


/************ GLOBAL FUNCTIONS ******************/

void
setheat (double hval1, double hval2, int heatmode)
{
  int ci;
  beta[0] = 1.0;
  for (ci = 1; ci < numchains; ci++)
  {
    switch (heatmode)
    {
    case HLINEAR:
      beta[ci] = 1.0 / (1.0 + hval1 * ci);
      break;
    case HTWOSTEP:
    case HGEOMETRIC:
      /* sometimes, with hval 1 values > 1 get identical adjacent beta values */ 
      beta[ci] =
        1 - (1 - hval2) * ci * pow (hval1, (double) (numchains - 1 - ci)) / (double) (numchains - 1);
      break;
    }
    //JH 6/11/2010  added this, mostly in case users use h1>1 (which may be ok but can also easily give negative values for beta)
    if (beta[ci] <= 0.0 || beta[ci] > 1.0)
      IM_err (IMERR_COMMANDLINEHEATINGTERMS, "command line heating terms have caused a heating value out of range. chain %d beta %lf", ci, beta[ci]);
  }
}                               /* setheat */

/* swaps chains,  adjust heating levels if adaptive heating is invoked.
will do multiple swaps if swaptries > 1.  If a swap was done involving chain0 then swap0ok =1. Sometimes, with 
swaptries > 1,  chain 0 will swap with another chain, and then swap back, restoring the current parameter 
values and genealogies. This is not detected, so stuff later that checks the return from this function
will think that parameters have changed.  This should not matter */
int
swapchains (int swaptries, int heatmode, double hval1, double hval2)
{
  int ci, cj, i, swap0ok;
  double metropolishastingsterm;
  void *swapptr;
  static int sw01 = 0, sw01d = 0, sw12 = 0, sw12d = 0;
  int cjmin, cjrange;

#define  MINSWAP  0.1
#define  BETADJUST  1.414
#define  INCADJUST  1.414
#define  MINHEAT  0.0001
#define  MAXHEAT  0.2
#define  MAXINC  100
#define  MININC  0.1
#define  PAUSESWAP 1000
#define SWAPDIST 7
// 5/27/2010  removed HADAPT stuff

  for (i = 0, swap0ok = 0; i < swaptries; i++)
  {

    do
    {
      ci = (int) (uniform () * numchains);
    } while (ci < 0 || ci >= numchains);

    if (numchains < 2*SWAPDIST + 3)
    {
      cjmin = 0;
      cjrange = numchains;
    }
    else
    {
      cjmin = IMAX(0,ci-SWAPDIST);
      cjrange = IMIN(numchains, ci+SWAPDIST) -cjmin;
    }
    do
    {
      cj = cjmin + (int) (uniform () * cjrange);
    } while (cj == ci || cj < 0 || cj >= numchains);
    if ((ci == 0 && cj == 1) || (ci == 1 && cj == 0))
      sw01d++;
    if ((ci == 1 && cj == 2) || (ci == 2 && cj == 1))
      sw12d++;
    if (ci < cj)
    {
      swapcount[cj][ci]++;
    }
    else
    {
      swapcount[ci][cj]++;
    }
    metropolishastingsterm = swapweight (ci, cj);
    if (metropolishastingsterm >= 1.0
        || metropolishastingsterm > uniform ())

    {
      swapptr = C[ci];
      C[ci] = C[cj];
      C[cj] = swapptr;
      if (ci < cj)

      {
        swapcount[ci][cj]++;
      }
      else
      {
        swapcount[cj][ci]++;
      }
      if ((ci == 0 && cj == 1) || (ci == 1 && cj == 0))
        sw01++;
      if ((ci == 1 && cj == 2) || (ci == 2 && cj == 1))
        sw12++;
      if (ci == 0 || cj == 0)
        swap0ok |= 1;
    }
  }
  return swap0ok;
}                               /* swapchains */


void
printchaininfo (FILE * outto, int heatmode, double hval1,
                double hval2)
{
  int i;
  fprintf (outto, "\nCHAIN SWAPPING BETWEEN SUCCESSIVE CHAINS: ");
  switch (heatmode)

  {
  case HLINEAR:
    fprintf (outto, " Linear Increment  term: %.4f\n", hval1);
    break;
  case HTWOSTEP:
    fprintf (outto, " Twostep Increment  term1: %.4f term2: %.4f\n", hval1,
             hval2);
    break;
  case HGEOMETRIC:
    fprintf (outto, " Geometric Increment  term1: %.4f term2: %.4f\n",
             hval1, hval2);
    break;
  }
  fprintf (outto,
           "-----------------------------------------------------------------------------\n");
  if (outto == stdout)

  {
    fprintf (outto, "beta terms :");
    for (i = 0; i < numchains; i++)
      fprintf (outto, "|%2d %5.3f", i, beta[i]);
    fprintf (outto, "\n");
    fprintf (outto, "Swap rates :");
    for (i = 0; i < numchains - 1; i++)
      if (swapcount[i + 1][i])
        fprintf (outto, "|%2d %5.3f", i,
                 swapcount[i][i + 1] / (float) swapcount[i + 1][i]);
    fprintf (outto, "\n");
    fprintf (outto, "Swap counts:");
    for (i = 0; i < numchains - 1; i++)
      fprintf (outto, "|%2d %5ld", i, swapcount[i][i + 1]);
    fprintf (outto, "\n\n");
  }
  else
  {
    fprintf (outto, "Chain   beta   #Swaps  Rate \n");
    for (i = 0; i < numchains - 1; i++)

    {
      if (swapcount[i + 1][i])
        fprintf (outto, " %3d  %7.4f %5ld  %7.4f\n", i, beta[i],
                 swapcount[i][i + 1],
                 swapcount[i][i + 1] / (float) swapcount[i + 1][i]);

      else
        fprintf (outto, " %3d  %7.4f %5ld  \n", i, beta[i],
                 swapcount[i][i + 1]);
    }
    fprintf (outto, " %3d  %7.4f   na      na\n", i, beta[i]);
    fprintf (outto, "\n");
  }
}                               /* printchaininfo */

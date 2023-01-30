TITLE skm95.mod  
 
COMMENT
----------------------------------------------------------------
Stochastic version of the K channel mechanism kd3h5.mod by
Z. Mainen in Mainen & Sejnowski 95.

This represents a potassium channel, with Hodgkin-Huxley like kinetics,
based on the gates model, assuming stochastic opening and closing.

Kinetic rates based roughly on Sah et al. and Hamill et al. (1991)
The main kinetic difference from the standard H-H model (shh.mod) is 
that the K+ kinetic is different, not n^4, but just n, 
and the activation curves are different.

The rate functions are adapted directly from the Kd3h5.mod file
by Zach Mainen.

The stochastic model is as following:

Potassium

       = alpha_n =>      
   [N0]             [N1]
      <= beta_n =      


The model keeps track on the number of channels in each state, and 
uses a binomial distribution to update these number.

Jan 1999, Mickey London, Hebrew University, mikilon@lobster.ls.huji.ac.il
        Peter N. Steinmetz, Caltech, peter@klab.caltech.edu
14 Sep 99 PNS. Added deterministic flag.
19 May 2002 Kamran Diba.  Changed gamma and deterministic from GLOBAL to RANGE.
23 Nov 2011 Werner Van Geit @ BBP. Changed the file so that it can use the neuron random number generator. Tuned voltage dependence
----------------------------------------------------------------
ENDCOMMENT

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
    SUFFIX StochKv
    THREADSAFE
    USEION k READ ek WRITE ik
    RANGE N,eta, gk, gamma, deterministic, gkbar, ik
    RANGE ninf, ntau,a,b,P_a,P_b
    GLOBAL Ra, Rb
    GLOBAL vmin, vmax, q10, temp
    RANGE tadj
    BBCOREPOINTER rng
}

UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
    (pS) = (picosiemens)
    (S) = (siemens)
    (um) = (micron)
} 

PARAMETER {
    v           (mV)
    dt      (ms)
    area    (um2)
    
    gamma  =  30          (pS)
    eta              (1/um2)
    gkbar = .75      (S/cm2)
    
    tha  = -40   (mV)        : v 1/2 for inf
    qa   = 9            : inf slope     
    Ra   = 0.02 (/ms)       : max act rate
    Rb   = 0.002    (/ms)       : max deact rate
    
    celsius (degC)
    temp = 23 (degC)   : original temperature for kinetic set
    q10 = 2.3               : temperature sensitivity
    
    deterministic = 0   : if non-zero, will use deterministic version
    vmin = -120 (mV)    : range to construct tables for
    vmax = 100  (mV)
} 

ASSIGNED {
    a       (/ms)
    b       (/ms)
    ik      (mA/cm2)
    gk      (S/cm2)
    ek      (mV)
    ninf        : steady-state value
    ntau (ms)   : time constant for relaxation
    tadj

    N 
    scale_dens (pS/um2) 
    P_a     : probability of one channel making alpha transition
    P_b     : probability of one channel making beta transition

    rng
    usingR123

    n0_n1_new

}


STATE {
    n         : state variable of deterministic description
    N0 N1       : N states populations
    n0_n1 n1_n0 : number of channels moving from one state to the other 
}

COMMENT
The Verbatim block is needed to generate random nos. from a uniform distribution between 0 and 1 
for comparison with Pr to decide whether to activate the synapse or not
ENDCOMMENT
   
VERBATIM
#include "nrnran123.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#ifndef CORENEURON_BUILD
double nrn_random_pick(void* r);
void* nrn_random_arg(int argpos);
#endif

ENDVERBATIM
: ----------------------------------------------------------------
: initialization
INITIAL {
    VERBATIM
    if( usingR123 ) {
        nrnran123_setseq((nrnran123_State*)_p_rng, 0, 0);
    }
    ENDVERBATIM
  
    eta = gkbar / gamma
    trates(v)
    n = ninf
    scale_dens = gamma/area
    N = floor(eta*area + 0.5)
    
    N1 = floor(n * N + 0.5)
    N0 = N-N1       : any round off into non-conducting state
    
    n0_n1 = 0
    n1_n0 = 0
}

: ----------------------------------------------------------------
: Breakpoint for each integration step
BREAKPOINT {
  SOLVE states
  
  gk =  (strap(N1) * scale_dens * tadj)
  
  ik = 1e-4 * gk * (v - ek)
} 


: ----------------------------------------------------------------
: states - updates number of channels in each state
PROCEDURE states() {

    trates(v)
    
    P_a = strap(a*dt)
    P_b = strap(b*dt)

    : check that will represent probabilities when used
    ChkProb( P_a)
    ChkProb( P_b)
    
    : transitions
    n0_n1 = BnlDev(P_a, N0)
    n1_n0 = BnlDev(P_b, N1)

    : move the channels
    N0    = strap(N0 - n0_n1 + n1_n0)
    N1    = N - N0
}

: ----------------------------------------------------------------
: trates - compute rates, using table if possible
PROCEDURE trates(v (mV)) {     
    TABLE ntau, ninf, a, b, tadj
    DEPEND dt, Ra, Rb, tha, qa, q10, temp, celsius
    FROM vmin TO vmax WITH 199
    
    tadj = q10 ^ ((celsius - temp)/(10 (K)))
    a = SigmoidRate(v, tha, Ra, qa)
    a = a * tadj
    b = SigmoidRate(-v, -tha, Rb, qa)
    b = b * tadj
    ntau = 1/(a+b)
    ninf = a*ntau
}


: ----------------------------------------------------------------
: SigmoidRate - Compute a sigmoid rate function given the 
: 50% point th, the slope q, and the amplitude a.
FUNCTION SigmoidRate(v (mV),th (mV),a (1/ms),q) (1/ms){
    UNITSOFF
    if (fabs(v-th) > 1e-6 ) {
        SigmoidRate = a * (v - th) / (1 - exp(-(v - th)/q))
    UNITSON

    } else {
        SigmoidRate = a * q
    }
}   


: ----------------------------------------------------------------
: sign trap - trap for negative values and replace with zero
FUNCTION strap(x) {
    if (x < 0) {
        strap = 0
VERBATIM
        fprintf (stderr,"skv.mod:strap: negative state");
ENDVERBATIM
    } else {
        strap = x
    }
}

: ----------------------------------------------------------------
: ChkProb - Check that number represents a probability
PROCEDURE ChkProb(p) {

  if (p < 0.0 || p > 1.0) {
    VERBATIM
// ToDo: should be disabled during ForwardSkip and enabled right after
//    fprintf(stderr, "StochKv.mod:ChkProb: argument not a probability.\n");
    ENDVERBATIM
  }

}

PROCEDURE setRNG() {

VERBATIM
    // For compatibility, allow for either MCellRan4 or Random123.  Distinguish by the arg types
    // Object => MCellRan4, seeds (double) => Random123
#ifndef CORENEURON_BUILD
    usingR123 = 0;
    if( ifarg(1) && hoc_is_double_arg(1) ) {
        nrnran123_State** pv = (nrnran123_State**)(&_p_rng);
        uint32_t a2 = 0;
        uint32_t a3 = 0;

        if (*pv) {
            nrnran123_deletestream(*pv);
            *pv = (nrnran123_State*)0;
        }
        if (ifarg(2)) {
            a2 = (uint32_t)*getarg(2);
        }
        if (ifarg(3)) {
            a3 = (uint32_t)*getarg(3);
        }
        *pv = nrnran123_newstream3((uint32_t)*getarg(1), a2, a3);
        usingR123 = 1;
    } else if( ifarg(1) ) {
        void** pv = (void**)(&_p_rng);
        *pv = nrn_random_arg(1);
    } else {
        void** pv = (void**)(&_p_rng);
        *pv = (void*)0;
    }
#endif
ENDVERBATIM
}

FUNCTION urand() {

VERBATIM
    double value;
    if( usingR123 ) {
        value = nrnran123_dblpick((nrnran123_State*)_p_rng);
    } else if (_p_rng) {
#ifndef CORENEURON_BUILD
        value = nrn_random_pick(_p_rng);
#endif
    } else {
        value = 0.5;
    }
    _lurand = value;
ENDVERBATIM
}

: Returns random numbers drawn from a binomial distribution
FUNCTION brand(P, N) {

VERBATIM
        /*
        :Supports separate independent but reproducible streams for
        : each instance. However, the corresponding hoc Random
        : distribution MUST be set to Random.uniform(0,1)
        */

        // Should probably be optimized
        double value = 0.0;
        int i;
        for (i = 0; i < _lN; i++) {
           if (urand(_threadargs_) < _lP) {
              value = value + 1;
           }
        }
        return(value);

ENDVERBATIM

        brand = value
}

VERBATIM
#define        PI 3.141592654
#define        r_ia     16807
#define        r_im     2147483647
#define        r_am     (1.0/r_im)
#define        r_iq     127773
#define        r_ir     2836
#define        r_ntab   32
#define        r_ndiv   (1+(r_im-1)/r_ntab)
#define        r_eps    1.2e-7
#define        r_rnmx   (1.0-r_eps)
ENDVERBATIM

VERBATIM
/* ---------------------------------------------------------------- */
/* gammln - compute natural log of gamma function of xx */
static double
gammln(double xx)
{
    double x,tmp,ser;
    static double cof[6]={76.18009173,-86.50532033,24.01409822,
        -1.231739516,0.120858003e-2,-0.536382e-5};
    int j;
    x=xx-1.0;
    tmp=x+5.5;
    tmp -= (x+0.5)*log(tmp);
    ser=1.0;
    for (j=0;j<=5;j++) {
        x += 1.0;
        ser += cof[j]/x;
    }
    return -tmp+log(2.50662827465*ser);
}
ENDVERBATIM


: ----------------------------------------------------------------
: BnlDev - draw a uniform deviate from the generator
FUNCTION BnlDev (ppr, nnr) {

VERBATIM
        int j;
        static int nold=(-1);
        double am,em,g,angle,p,bnl,sq,bt,y;
        static double pold=(-1.0),pc,plog,pclog,en,oldg;
        
        /* prepare to always ignore errors within this routine */
         
        
        p=(_lppr <= 0.5 ? _lppr : 1.0-_lppr);
        am=_lnnr*p;
        if (_lnnr < 25) {
            bnl=0.0;
            for (j=1;j<=_lnnr;j++)
                if (urand(_threadargs_) < p) bnl += 1.0;
        }
        else if (am < 1.0) {
            g=exp(-am);
            bt=1.0;
            for (j=0;j<=_lnnr;j++) {
                bt *= urand(_threadargs_);
                if (bt < g) break;
            }
            bnl=(j <= _lnnr ? j : _lnnr);
        }
        else {
            if (_lnnr != nold) {
                en=_lnnr;
                oldg=gammln(en+1.0);
                nold=_lnnr;
            }
            if (p != pold) {
                pc=1.0-p;
                 plog=log(p);
                pclog=log(pc);
                pold=p;
            }
            sq=sqrt(2.0*am*pc);
            do {
                do {
                    angle=PI*urand(_threadargs_);
                        angle=PI*urand(_threadargs_);
                    y=tan(angle);
                    em=sq*y+am;
                } while (em < 0.0 || em >= (en+1.0));
                em=floor(em);
                    bt=1.2*sq*(1.0+y*y)*exp(oldg-gammln(em+1.0) - 
                    gammln(en-em+1.0)+em*plog+(en-em)*pclog);
            } while (urand(_threadargs_) > bt);
            bnl=em;
        }
        if (p != _lppr) bnl=_lnnr-bnl;
        
        /* recover error if changed during this routine, thus ignoring
            any errors during this routine */
       
        
        return bnl;
        
    ENDVERBATIM
    BnlDev = bnl
}  

VERBATIM
static void bbcore_write(double* x, int* d, int* xx, int* offset, _threadargsproto_) {
    if (d) {
        uint32_t* di = ((uint32_t*)d) + *offset;
      // temporary just enough to see how much space is being used
      if (!_p_rng) {
        di[0] = 0; di[1] = 0, di[2] = 0;
      }else{
        nrnran123_State** pv = (nrnran123_State**)(&_p_rng);
        nrnran123_getids3(*pv, di, di+1, di+2);
        // write stream sequence
        char which;
        nrnran123_getseq(*pv, di+3, &which);
        di[4] = (int)which;
      }
      //printf("StochKv3.mod %p: bbcore_write offset=%d %d %d\n", _p, *offset, d?di[0]:-1, d?di[1]:-1);
    }
    *offset += 5;
}
static void bbcore_read(double* x, int* d, int* xx, int* offset, _threadargsproto_) {
    assert(!_p_rng);
    uint32_t* di = ((uint32_t*)d) + *offset;
        if (di[0] != 0 || di[1] != 0|| di[2] != 0)
        {
      nrnran123_State** pv = (nrnran123_State**)(&_p_rng);
      *pv = nrnran123_newstream3(di[0], di[1], di[2]);
      nrnran123_setseq(*pv, di[3], (char)di[4]);
        }
      //printf("StochKv3.mod %p: bbcore_read offset=%d %d %d\n", _p, *offset, di[0], di[1]);
    *offset += 5;
}
ENDVERBATIM

FUNCTION bbsavestate() {
        bbsavestate = 0
VERBATIM
 #ifndef CORENEURON_BUILD
        // TODO: since N0,N1 are no longer state variables, they will need to be written using this callback
        //  provided that it is the version that supports multivalue writing
        /* first arg is direction (-1 get info, 0 save, 1 restore), second is value*/
        double *xdir, *xval, *hoc_pgetarg();
        long nrn_get_random_sequence(void* r);
        void nrn_set_random_sequence(void* r, int val);
        xdir = hoc_pgetarg(1);
        xval = hoc_pgetarg(2);
        if (_p_rng) {
                // tell how many items need saving
                if (*xdir == -1.) {
                    if( usingR123 ) {
                        *xdir = 2.0;
                    } else {
                        *xdir = 1.0;
                    }
                    return 0.0;
                }
                else if (*xdir == 0.) {
                    if( usingR123 ) {
                        uint32_t seq;
                        char which;
                        nrnran123_getseq( (nrnran123_State*)_p_rng, &seq, &which );
                        xval[0] = (double) seq;
                        xval[1] = (double) which;
                    } else {
                        xval[0] = (double)nrn_get_random_sequence(_p_rng);
                    }
                } else{
                    if( usingR123 ) {
                        nrnran123_setseq( (nrnran123_State*)_p_rng, (uint32_t)xval[0], (char)xval[1] );
                    } else {
                        nrn_set_random_sequence(_p_rng, (long)(xval[0]));
                    }
                }
        }

        // TODO: check for random123 and get the seq values
#endif
ENDVERBATIM
}


/* Produced by f2c and Martin's  f2c-clean,v 1.11 2012/05/04 --
   and some pretty editing
*/
#include <math.h>

void cldaisy(int *nn, int *jpp, double *x,
	     double *valmd, double *weights, int *jtmd, int *jdat,
	     int *vtype, int *ndyst, int *mdata, double *disv)
{
/* Calculating dissimilarities between objects or variables

            nn  = number of objects
            jpp = number of variables used for the calculations

   The following vectors and matrices must be dimensioned in the main program :

       vtype was character originally
       vtype(j) is the type of variable j:
              = 1 (A) for an Asymmetric binary variable
              = 2 (S) for a  Symmetric  binary variable
              = 3 (N) for a  Nominal  variable
              = 4 (O) for an Ordinal  variable
              = 5 (I) for an Interval variable (additive)
              = 6 (T) for a  raTio    variable (log transformed)
     ndyst is the "kind of dissimilarity/distance"  aka 'diss_type'
              = 0  "mixed" / gower == here treated _identically_ to "manhattan" (L1) !
              = 1  "euclidean" (= L2 )
              = 2  "manhattan" (= L1 )
       vector jtmd is only read if there are missing values : if(mdata)
       jtmd(j) =  0 if variable j is binary
               = -1 if variable j is not binary and has missing values
               = +1 if variable j is not binary and has no missing values
*/

    /* int nbad;  */

/* VARs Parameter adjustments --- all this just so we can use 1-indexing  (FIXME?) */
    --disv;
    --vtype;
    --jtmd;
    --weights;
    --valmd;
    int
	x_dim1 = *nn,
	x_offset = 1 + x_dim1;
    x -= x_offset;

    /* Function Body */
    /*logical*/int hasna = *mdata != 0;

/*         calculation of the dissimilarities */
    int nlk = 0;
    if (*jdat == 1) {
	/* nbad=0; */

/* Case I: `mixed' type variables  <==> ndyst == 0 */
	for (int l = 2; l <= *nn; ++l) {
	    int la = l - 1;
	    for (int k = 1; k <= la; ++k) {
		++nlk;
		/*  disv(nlk) :=  d(l,k) := Dissimilarity between obs. l & k := dlk/ppa */
		double ppa = 0., dlk = 0.;
		for (int j = 1; j <= *jpp; ++j) {
		    if (vtype[j] >= 3) {
			if (hasna) {
			    if (jtmd[j] < 0) {
				if ((x[l + j * x_dim1] == valmd[j]) ||
				    (x[k + j * x_dim1] == valmd[j])) {
				    continue; // goto L420;
				}
			    }
			}
			ppa += weights[j];
			if (vtype[j] == 3) { /* type = "N"ominal */
			    if (x[l + j * x_dim1] != x[k + j * x_dim1]) {
				dlk += weights[j];
			    }
			} else { /* type in {"O", "I", "T"} */
			    dlk += weights[j] * fabs(x[l + j * x_dim1] - x[k + j * x_dim1]);
			}
		    } else { /* vtype(j) in {1, 2} <---> binary variable x(*,j) */
			if((x[l + j * x_dim1] != 0.f && x[l + j * x_dim1] != 1.f) ||
			   (x[k + j * x_dim1] != 0.f && x[k + j * x_dim1] != 1.f)) {
			    continue; // goto L420;
			}
			if (vtype[j] == 2 || x[l + j * x_dim1] != 0. || x[k + j * x_dim1] != 0.) {
			    ppa += weights[j];
			}
			if (x[l + j * x_dim1] != x[k + j * x_dim1]) {
			    dlk += weights[j];
			}
		    }
		    /* L420: -- continue -- */
		}
		if (ppa == 0.f) {
/* was          if(ppa <= 0.5) then
                 nbad=nbad+1 */
		    disv[nlk] = -1.;
		} else {
		    disv[nlk] = dlk / ppa;
		}
	    }
	}
    } else
    {
/* Case II : jdat != 1:  all variables are interval scaled
   -------   ~~~~~~~~~ { basically === dysta() in ./dysta.f
                         FIXME: common code! } */
	double pp = (double) (*jpp);
	for (int l = 2; l <= *nn; ++l) {
	    int lsubt = l - 1;
	    for (int k = 1; k <= lsubt; ++k) {
		double clk = 0.;
		++nlk;
		int npres = 0;
		/*  disv(nlk) :=  d(l,k) := Dissimilarity between obs.  l & k */
		for (int j = 1; j <= *jpp; ++j) {
		    if (hasna) {
			if (jtmd[j] < 0) {
			    if ((x[l + j * x_dim1] == valmd[j]) ||
				(x[k + j * x_dim1] == valmd[j])) {
				continue; // goto L530;
			    }
			}
		    }
		    ++npres;
		    if (*ndyst == 1) { /* euclidean */
			clk += (x[l + j * x_dim1] - x[k + j * x_dim1]) *
			       (x[l + j * x_dim1] - x[k + j * x_dim1]);
		    } else { /* manhattan */
			clk += fabs(x[l + j * x_dim1] - x[k + j * x_dim1]);
		    }
		    /* L530: ; */
		}
		double rpres = (double) npres;
		if (npres == 0) { /* all missing */
		    disv[nlk] = -1.f;
		} else if (*ndyst == 1) {
		    disv[nlk] = sqrt(clk * (pp / rpres));
		} else {
		    disv[nlk] = clk * (pp / rpres);
		}
	    }
	}
    }
} /* cldaisy */


/* -*- mode: c; kept-new-versions: 25; kept-old-versions: 20 -*- */

/*   Clustering LARge Applications
     ~		~~~   ~
     Clustering program based upon the k-medoid approach,
     and suitable for data sets of at least 100 objects.
     (for smaller data sets, please use program pam.)
 */

/* $Id$
 * original Id: clara.f,v 1.10 2002/08/27 15:43:58 maechler translated by
 * f2c (version 20010821) and run through f2c-clean,v 1.10 2002/03/28
 */

#include <math.h>

#include <R_ext/Print.h>/* for diagnostics */

#include "cluster.h"

void clara(int *n,  /* = number of objects */
	   int *jpp,/* = number of variables */
	   int *kk, /* = number of clusters, 1 <= kk <= n-1 */
	   double *x,	/* Input:  the data x[n, jpp] _rowwise_ (transposed)
			 * Output: the first `n' values are the `clustering'
			 *	   (integers in 1,2,..,kk) */
	   int *nran,	/* = #{random samples} drawn	   (= `samples' in R)*/
	   int *nsam,	/* = #{objects} drawn from data set (`sampsize' in R) */
	   double *dys,/* [1:(1 + (nsam * (nsam - 1))/2)]
			* Output: to contain the distances */
	   int *mdata,	/*= {0,1}; 1: min(x) is missing value (NA);  0: no NA */
	   double *valmd,/*[j]= missing value code (instead of NA) for x[,j]*/
	   int *jtmd,	/* [j]= {-1,1};	 -1: x[,j] has NA; 1: no NAs in x[,j] */
	   int *diss_kind,/* = {1,2};  1 : euclidean;  2 : manhattan*/
	   int *nrepr,
	   int *nsel,
	   int *nbest,/* x[nbest[j]] will be the j-th obs in the final sample */
	   int *nr, int *nrx,
	   double *radus, double *ttd, double *ratt,
	   double *ttbes, double *rdbes, double *rabes,
	   int *mtt, double *obj,
	   double *avsyl, double *ttsyl, double *sylinf,
	   int *jstop, int *trace_lev,
	   double *tmp, /* = double [ 3 * nsam ] */
	   int *itmp	/* = integer[ 6 * nsam ] */
    )
{

#define tmp1 tmp
#define tmp2 &tmp[*nsam]

#define ntmp1 itmp
#define ntmp2 &itmp[*nsam]
#define ntmp3 &itmp[nsamb]
#define ntmp4 &itmp[nsamb+ *nsam]
#define ntmp5 &itmp[2*nsamb]
#define ntmp6 &itmp[2*nsamb+ *nsam]

    /* Local variables */

    Rboolean nafs, kall, full_sample, lrg_sam, has_NA = *mdata;
    int j, jk, jkk, js, jsm, jhalt, jran, l, n_sam;
    int nsm, ntt, nad, nadv, kran, kans, nrun, nexap, nexbp;
    int n_dys, nsamb, nunfs;
    double rnn, sky, zb, s, sx = -1., zba = -1.;/* Wall */

    /* Parameter adjustments */
    --nsel;

    *jstop = 0;
    rnn = (double) (*n);
    /* n_dys := size of distance array dys[] */
    n_dys = *nsam * (*nsam - 1) / 2 + 1;/* >= 1 */

    full_sample = (*n == *nsam);/* only one sub sample == full data */
    nsamb = *nsam * 2;
    lrg_sam = (*n < nsamb);/* sample more than *n/2 */
    if (lrg_sam)/* generate indices for the other, smaller half */
	n_sam = *n - *nsam;
    else
	n_sam = *nsam;
    nunfs = 0;
    kall = FALSE;

    nrun = 0; /* << initialize `random seed' of the very simple randm() below */

/* __LOOP__ :  random subsamples are drawn and partitioned into kk clusters */

    for (jran = 1; jran <= *nran; ++jran) {
	jhalt = 0;
	if (!full_sample) {/* `real' case: sample size < n */
	    ntt = 0;
	    if (kall/*was jran != 1 */ && nunfs != jran && !lrg_sam) {
		/* nsel[] := sort(nrx[])   for the first j=1:k	?? */
		for (jk = 0; jk < *kk; ++jk)
		    nsel[jk+1] = nrx[jk];
		for (jk = 1; jk < *kk; ++jk) {
		    nsm = nsel[jk];
		    jsm = jk;
		    for (jkk = jk + 1; jkk <= *kk; ++jkk) {
			if (nsel[jkk] < nsm) {
			    nsm = nsel[jkk];
			    jsm = jkk;
			}
		    }
		    nsel[jsm] = nsel[jk];
		    nsel[jk] = nsm;
		}
		ntt = *kk;
	    }
	    else {
		/* Loop finding random index `kran' not yet in nrx[] : */
	    L180:
		kran = (int) (rnn * randm(&nrun) + 1.);
		if(*trace_lev >= 3)
		    Rprintf("... {180} nrun=%d -> k{ran}=%d\n", nrun,kran);
		if (kran > *n) kran = *n;

		if (kall/*jran != 1*/) {
		    for (jk = 0; jk < *kk; ++jk)
			if (kran == nrx[jk])
			    goto L180;
		}
		/* end Loop */
		nsel[++ntt] = kran;
		if (ntt == n_sam)
		    goto L295;
	    }

	    if(*trace_lev >= 2) {
		Rprintf(".. kall(T/F)=%d , nsel[ntt=%d] = %d\n",
			kall, ntt, nsel[ntt]);
		if(*trace_lev >= 3) {
		    Rprintf("... nrx[0:%d]= ",*kk-1);
		    for (jk = 0; jk < *kk; jk++)
			Rprintf("%d ",nrx[jk]); Rprintf("\n");
		    Rprintf("... nsel[1:%d]= ",ntt);
		    for (jk = 1; jk <= ntt; jk++)
			Rprintf("%d ",nsel[jk]); Rprintf("\n");
		}
	    }

	    do {
	    L210:
		/* find `kran', a random `k' in {1:n},
		 * not in nrx[0:(k-1)] nor nsel[1:ntt] : */
		kran = (int) (rnn * randm(&nrun) + 1.);
		if(*trace_lev >= 3)
		    Rprintf("... {210} nrun=%d -> k{ran}=%d\n", nrun,kran);
		if (kran > *n) kran = *n;

		if (kall/*jran != 1*/ && lrg_sam) {
		    for (jk = 0; jk < *kk; ++jk) {
			if (kran == nrx[jk])
			    goto L210;
		    }
		}
		/* insert kran into nsel[1:ntt] or after  and increase ntt : */
		for (kans = 1; kans <= ntt; ++kans)
		    if (nsel[kans] >= kran) {
			if (nsel[kans] == kran)
			    goto L210;
			else {
			    for (nad = kans; nad <= ntt; ++nad) {
				nadv = ntt - nad + kans;
				nsel[nadv + 1] = nsel[nadv];
			    }
			    nsel[kans] = kran;
			    /* continue _outer_ loop */ goto L290;
			}
		    }
		nsel[ntt+1] = kran;

	    L290:
		++ntt;
	    } while (ntt < n_sam);

	L295:
	    if (lrg_sam) {
		/* have indices for smaller _nonsampled_ half; revert this: */
		for (j = 1, nexap = 1, nexbp = 0; j <= *n; j++) {
		    if (nsel[nexap] == j)
			++nexap;
		    else
			nrepr[nexbp++] = j;
		}
		for (j = 0; j < *nsam; ++j)
		    nsel[j+1] = nrepr[j];
	    }
	    if(*trace_lev)
		Rprintf("C clara(): sample %d [ntt=%d, nunfs=%d] -> dysta2()\n",
			jran, ntt, nunfs);
	}
	else { /* full_sample : *n = *nsam -- one sample is enough ! */
	    for (j = 1; j <= *nsam; ++j)
		nsel[j] = j;
	}

	dysta2(*nsam, *jpp, &nsel[1], x, *n, dys, *diss_kind,
	       jtmd, valmd, &jhalt);
	if (jhalt == 1) {
	    if(*trace_lev)
		Rprintf("  dysta2() gave jhalt = 1 --> new sample\n");
	    continue;/* random sample*/
	}

	s = 0.;
	l = 0;/* dys[0] is not used here */
	do {
	    ++l;
	    if (s < dys[l])
		s = dys[l];
	} while (l+1 < n_dys);

	if(*trace_lev >= 2)
	    Rprintf(". clara(): s:= max dys[l=%d] = %g\n", l,s);

	bswap2(*kk, *nsam, nrepr, dys, &sky, s,
	       /* dysma */tmp1, /*dysmb*/tmp2,
	       /* beter[], only used here */&tmp[nsamb]);

	selec(*kk, *n, *jpp, *diss_kind, &zb, *nsam, has_NA, jtmd, valmd,
	      nrepr, &nsel[1], dys, x, nr, &nafs, ttd, radus, ratt,
	      ntmp1, ntmp2, ntmp3, ntmp4, ntmp5, ntmp6, tmp1, tmp2);

	if (nafs) {
	    ++nunfs;
	}
	else if(!kall /*jran == 1*/ || zba > zb) {
	    /* 1st proper sample  or  new best */
	    zba = zb;
	    for (jk = 0; jk < *kk; ++jk) {
		ttbes[jk] = ttd	 [jk];
		rdbes[jk] = radus[jk];
		rabes[jk] = ratt [jk];
		nrx  [jk] = nr	 [jk];
	    }
	    for (js = 0; js < *nsam; ++js)
		nbest[js] = nsel[js+1];
	    sx = s;
	}

	kall = TRUE;

	if(full_sample) break; /* out of resampling */
    }
/* --- end random sampling loop */

    if (nunfs >= *nran) { *jstop = 1; return; }


/*     for the best subsample, the objects of the entire data set
     are assigned to their clusters */

    if (!kall) { *jstop = 2; return; }

    if(*trace_lev) {
	Rprintf("C clara(): sample _found_ --> dysta2(nbest), resul(), end\n");
    }
    *obj = zba / rnn;
    dysta2(*nsam, *jpp, nbest, x, *n, dys, *diss_kind, jtmd, valmd, &jhalt);

    resul(*kk, *n, *jpp, *diss_kind, has_NA, jtmd, valmd, x, nrx, mtt);

    if (*kk > 1)
	black(*kk, *jpp, *nsam, nbest, dys, sx, x,
	      /* compute --> */
	      avsyl, ttsyl, sylinf,
	      ntmp1, ntmp2, ntmp3, ntmp4, tmp1, tmp2);
    return;
} /* End clara() ---------------------------------------------------*/
#undef tmp1
#undef tmp2

#undef ntmp1
#undef ntmp2
#undef ntmp3
#undef ntmp4
#undef ntmp5
#undef ntmp6


void dysta2(int nsam, int jpp, int *nsel,
	    double *x, int n, double *dys, int diss_kind,
	    int *jtmd, double *valmd, int *jhalt)
{
/* Compute Dissimilarities for the selected sub-sample	---> dys[,] */

    /* Local variables */
    int j, k, kj, l, lj, ksel, lsel, nlk, npres;
    double clk, d1;

    nlk = 0;
    dys[0] = 0.;/* very first index; *is* used in ?? (it is, not in R!) */
    for (l = 1; l < nsam; ++l) {
	lsel = nsel[l];
	if(lsel <= 0 || lsel > n)
	    REprintf(" ** dysta2(): nsel[l= %d] = %d is OUT\n", l, lsel);
	for (k = 0; k < l; ++k) {
	    ksel = nsel[k];
	    if(ksel <= 0 || ksel > n)
		REprintf(" ** dysta2(): nsel[k= %d] = %d is OUT\n", k, ksel);
	    clk = 0.;
	    ++nlk;
	    npres = 0;
	    for (j = 0; j < jpp; ++j) {
		lj = (lsel - 1) * jpp + j;
		kj = (ksel - 1) * jpp + j;
		if (jtmd[j] < 0) {
		    /* in the following line (Fortran!), x[-2] ==> seg.fault
		       {BDR to R-core, Sat, 3 Aug 2002} */
		    if (x[lj] == valmd[j]) {
			continue /* next j */;
		    }
		    if (x[kj] == valmd[j]) {
			continue /* next j */;
		    }
		}
		++npres;
		if (diss_kind == 1)
		    clk += (x[lj] - x[kj]) * (x[lj] - x[kj]);
		else
		    clk += fabs(x[lj] - x[kj]);
	    }
	    if (npres == 0) {
		*jhalt = 1;
		dys[nlk] = -1.;
	    } else {
		d1 = clk * (((double) (jpp)) / (double) npres);
		dys[nlk] = (diss_kind == 1) ? sqrt(d1) : d1 ;
	    }
	}
    }
    return;
} /* End dysta2() -----------------------------------------------------------*/

double randm(int *nrun)
{
/* we programmed this generator ourselves because we wanted it
   to be machine independent. it should run on most computers
   because the largest int used is less than 2^30 . the period
   is 2^16=65536, which is good enough for our purposes. */
    /* MM: improved the original speed-wise only: */
    *nrun = (*nrun * 5761 + 999) & 0177777;
    /* Masking off all but the last 16 bits is equivalent to  % 65536 */
    return ((double) (*nrun) / 65536.);
} /* randm() */

/* bswap2() : called once [per random sample] from clara() : */
void bswap2(int kk, int nsam, int *nrepr,
	    double *dys, double *sky, double s,
	    double *dysma, double *dysmb, double *beter)
{
    int j, ja, k, kbest = -1, nbest = -1;/* init for -Wall */
    int nkj, njn, njaj, nny, nmax;

    double ammax, small, asky, dzsky, dz, cmd;

    /* Parameter adjustments */
    --nrepr;
    --beter;

    --dys;/*index via meet() */
    --dysma;	--dysmb;


/* ====== first algorithm: BUILD. ====== */

    for (j = 1; j <= nsam; ++j) {
	nrepr[j] = 0;
	dysma[j] = s * 1.1 + 1.;
    }

    for(nny = 0; nny < kk; nny++) {
	for (ja = 1; ja <= nsam; ++ja) {
	    if (nrepr[ja] == 0) {
		beter[ja] = 0.;
		for (j = 1; j <= nsam; ++j) {
		    njaj = F77_CALL(meet)(&ja, &j);
		    cmd = dysma[j] - dys[njaj];
		    if (cmd > 0.)
			beter[ja] += cmd;
		}
	    }
	}
	ammax = 0.;
	for (ja = 1; ja <= nsam; ++ja) {
	    if (nrepr[ja] == 0 && ammax <= beter[ja]) {
		ammax = beter[ja];
		nmax = ja;
	    }
	}
	nrepr[nmax] = 1;
	for (j = 1; j <= nsam; ++j) {
	    njn = F77_CALL(meet)(&nmax, &j);
	    if (dysma[j] > dys[njn])
		dysma[j] = dys[njn];
	}
    }

    *sky = 0.;
    for (j = 1; j <= nsam; ++j)
	*sky += dysma[j];

    if (kk == 1)
	return;

    asky = *sky / ((double) nsam);

/* ====== second algorithm: SWAP. ====== */

/* Big LOOP : */
L60:

    for (j = 1; j <= nsam; ++j) {
	dysma[j] = s * 1.1 + 1.;
	dysmb[j] = s * 1.1 + 1.;
	for (ja = 1; ja <= nsam; ++ja) {
	    if (nrepr[ja] != 0) {
		njaj = F77_CALL(meet)(&ja, &j);
		if (dys[njaj] < dysma[j]) {
		    dysmb[j] = dysma[j];
		    dysma[j] = dys[njaj];
		} else if (dys[njaj] < dysmb[j])
		    dysmb[j] = dys[njaj];
	    }
	}
    }

    dzsky = 1.;
    for (k = 1; k <= nsam; ++k) {
	if (nrepr[k] != 1) {
	    for (ja = 1; ja <= nsam; ++ja) {
		if (nrepr[ja] != 0) {
		    dz = 0.;
		    for (j = 1; j <= nsam; ++j) {
			njaj = F77_CALL(meet)(&ja, &j);
			nkj  = F77_CALL(meet)(&k, &j);
			if (dys[njaj] == dysma[j]) {
			    small = dysmb[j];
			    if (small > dys[njaj])
				small = dys[nkj];
			    dz = dz - dysma[j] + small;
			}
			else if (dys[nkj] < dysma[j])
			    dz = dz - dysma[j] + dys[nkj];
		    }
		    if (dz < dzsky) {
			dzsky = dz;
			kbest = k;
			nbest = ja;
		    }
		}
	    }
	}
    }
    if (dzsky >= 0.)
	return;

    nrepr[kbest] = 1;
    nrepr[nbest] = 0;
    *sky += dzsky;
    goto L60;

} /* End of bswap2() -------------------------------------------------- */

/* selec() : called once [per random sample] from clara() */
void selec(int kk, int n, int jpp, int diss_kind,
	   double *zb, int nsam, Rboolean has_NA, int *jtmd, double *valmd,
	   int *nrepr, int *nsel, double *dys, double *x, int *nr,
	   Rboolean *nafs, /* := TRUE if a distance cannot be calculated */
	   double *ttd, double *radus, double *ratt,
	   int *nrnew, int *nsnew, int *npnew, int *ns, int *np, int *new,
	   double *ttnew, double *rdnew)
{

    /* Local variables */
    int j, jk, jj, jp, jnew, ka, kb, jkabc = -1/* -Wall */;
    int newf, nrjk,  npab, nstrt, na, nb, npa, npb, njk, nobs;

    double dsum, pp, tra, rns, dnull = -9./* -Wall */;

/* Parameter adjustments */
    --nsel;    --nrepr;

    --ratt;
    --radus; --ttd;    --np;	--nr;	 --ns;

    --rdnew; --ttnew; --npnew; --nrnew; --nsnew;
    --new;

    --dys;


    /* nafs := TRUE  if a distance cannot be calculated*/
    *nafs = FALSE;

/* identification of representative objects, and initializations */

    jk = 0;
    for (j = 1; j <= nsam; ++j) {
	if (nrepr[j] != 0) {
	    ++jk;
	    nr	 [jk] = nsel[j];
	    ns	 [jk] = 0;
	    ttd	 [jk] = 0.;
	    radus[jk] = -1.;
	    np	 [jk] = j;
	}
    }

/* - assignment of the objects of the entire data set to a cluster,
 * - computation of some statistics,
 * - determination of the new ordering of the clusters */

    *zb = 0.;
    pp = (double) (jpp);
    newf = 0;

    for(jj = 1; jj <= n; jj++) {
	if (!has_NA) {
	    for (jk = 1; jk <= kk; ++jk) {
		dsum = 0.;
		nrjk = nr[jk];
		if (nrjk != jj) {
		    for (jp = 0; jp < jpp; ++jp) {
			na = (nrjk - 1) * jpp + jp;
			nb = (jj   - 1) * jpp + jp;
			tra = fabs(x[na] - x[nb]);
			if (diss_kind == 1)
			    tra *= tra;
			dsum += tra;
		    }
		    if (jk != 1 && dsum >= dnull)
			continue /* next jk */;
		}
		dnull = dsum;
		jkabc = jk;
	    }
	}
	else { /* _has_ missing data */
	    Rboolean pres;
	    pres = FALSE;
	    for (jk = 1; jk <= kk; ++jk) {
		dsum = 0.;
		nrjk = nr[jk];
		if (nrjk != jj) {
		    nobs = 0;
		    for (jp = 0; jp < jpp; ++jp) {
			na = (nrjk - 1) * jpp + jp;
			nb = (jj   - 1) * jpp + jp;
			if (jtmd[jp] < 0) {
			    if (x[na] == valmd[jp] || x[nb] == valmd[jp])
				continue /* next jp */;
			}
			nobs++;
			tra = fabs(x[na] - x[nb]);
			if (diss_kind == 1)
			    tra *= tra;
			dsum += tra;
		    }
		    if (nobs == 0) /* all pairs partially missing */
			continue /* next jk */;
		    dsum *= (nobs / pp);
		}
		if (!pres)
		    pres = TRUE;
		else if (dnull <= dsum)
		    continue /* next jk */;
		/* here : pres was FALSE {i.e. 1st time} or
		 *	  dnull > dsum	 {i.e. new best} */
		dnull = dsum;
		jkabc = jk;
	    }/* for(jk ..) */

	    if (!pres) { /* found nothing */
		*nafs = TRUE; return;
	    }
	} /* if (has_NA..) else */

	if (diss_kind == 1)
	    dnull = sqrt(dnull);

	*zb += dnull;
	ttd[jkabc] += dnull;
	if (radus[jkabc] < dnull)
	    radus[jkabc] = dnull;

	++ns[jkabc];
	if (newf < kk) {
	    if (newf != 0) {
		for (jnew = 1; jnew <= newf; ++jnew) {
		    if (jkabc == new[jnew])
			goto L90;/* next jj */
		}
	    }
	    ++newf;
	    new[newf] = jkabc;
	}
    L90:
	;
    } /* for( jj = 1..n ) */


/*     a permutation is carried out on vectors nr,ns,np,ttd,radus
     using the information in vector new. */

    for (jk = 1; jk <= kk; ++jk) {
	njk = new[jk];
	nrnew[jk] = nr[njk];
	nsnew[jk] = ns[njk];
	npnew[jk] = np[njk];
	ttnew[jk] = ttd[njk];
	rdnew[jk] = radus[njk];
    }
    for (jk = 1; jk <= kk; ++jk) {
	nr[jk] = nrnew[jk];
	ns[jk] = nsnew[jk];
	np[jk] = npnew[jk];
	ttd[jk] = ttnew[jk];
	radus[jk] = rdnew[jk];
    }
    for (j = 1; j <= kk; ++j) {
	rns = (double) ns[j];
	ttd[j] /= rns;
    }

    if (kk > 1) {

	/* computation of minimal distance of medoid ka to any
	   other medoid for comparison with the radius of cluster ka. */

	for (ka = 1; ka <= kk; ++ka) {
	    nstrt = 0;
	    npa = np[ka];
	    for (kb = 1; kb <= kk; ++kb) {
		if (kb == ka)
		    continue /* next kb */;

		npb = np[kb];
		npab = F77_CALL(meet)(&npa, &npb);
		if (nstrt == 0)
		    nstrt = 1;
		else if (dys[npab] >= ratt[ka])
		    continue /* next kb */;

		ratt[ka] = dys[npab];
		if (ratt[ka] != 0.)
		    continue /* next kb */;

		ratt[ka] = -1.;
	    }
	    if (ratt[ka] > -0.5)
		ratt[ka] = radus[ka] / ratt[ka];
	}
    }
    return;
} /* End selec() -----------------------------------------------------------*/

void resul(int kk, int n, int jpp, int diss_kind, Rboolean has_NA,
	   int *jtmd, double *valmd, double *x, int *nrx, int *mtt)
{
    /* Local variables */
    int j, jk, jj, ka, na, nb, njnb, nrjk, jksky = -1/* Wall */;
    double pp, abc, dsum, tra, dnull = -9./* Wall */;

/* clustering vector is incorporated into x, and ``printed''. */

    pp = (double) (jpp);

    for(jj = 0; jj < n; jj++) {

	for (jk = 0; jk < kk; ++jk) {
	    if (nrx[jk] == jj + 1)/* 1-indexing */
		goto L220; /* continue next jj (i.e., outer loop) */
	}
	njnb = jj * jpp;

	if (!has_NA) {
	    for (jk = 0; jk < kk; ++jk) {
		dsum = 0.;
		nrjk = (nrx[jk] - 1) * jpp;
		for (j = 0; j < jpp; ++j) {
		    na = nrjk + j;
		    nb = njnb + j;
		    tra = fabs(x[na] - x[nb]);
		    if (diss_kind == 1)
			tra *= tra;
		    dsum += tra;
		}
		if (diss_kind == 1)
		    dsum = sqrt(dsum);
		if (jk == 0)
		    dnull = dsum + .1f;
		if (dnull > dsum) {
		    dnull = dsum;
		    jksky = jk;
		}
	    }
	}
	else { /* _has_ missing data */
	    for (jk = 0; jk < kk; ++jk) {
		dsum = 0.;
		nrjk = (nrx[jk] - 1) * jpp;
		abc = 0.;
		for (j = 0; j < jpp; ++j) {
		    na = nrjk + j;
		    nb = njnb + j;
		    if (jtmd[j] < 0) {
			if (x[na] == valmd[j] || x[nb] == valmd[j])
			    continue /* next j */;
		    }
		    abc += 1.;
		    tra = fabs(x[na] - x[nb]);
		    if (diss_kind == 1)
			tra *= tra;
		    dsum += tra;
		}
		if (diss_kind == 1)
		    dsum = sqrt(dsum);
		dsum *= (abc / pp);
		if (jk == 0)
		    dnull = dsum + .1f;

		if (dnull > dsum) {
		    dnull = dsum;
		    jksky = jk;
		}
	    }
	}
	x[njnb] = (double) jksky + 1;/* 1-indexing */

    L220:
	;
    } /* for(jj)  while (jj < n);*/

    for (jk = 0; jk < kk; ++jk) {
	nrjk = nrx[jk];
	x[(nrjk - 1) * jpp] = (double) jk + 1;/* 1-indexing */
    }
    for (ka = 0; ka < kk; ++ka) {
	mtt[ka] = 0;
	for(j = 0; j < n; j++) {
	    if (((int) (x[j * jpp] + .1f)) == ka + 1)/* 1-indexing */
		++mtt[ka];
	}
    }
    return;
} /* end resul() -----------------------------------------------------------*/


void black(int kk, int jpp, int nsam, int *nbest,
	   double *dys, double s, double *x,
	   /* --> Output : */
	   double *avsyl, double *ttsyl, double *sylinf,
	   int *ncluv, int *nsend, int *nelem, int *negbr,
	   double *syl, double *srank)
{
/* Silhouettes computation and "drawing"  --> syl[] and sylinf[] */

    /* System generated locals */
    int sylinf_dim1, sylinf_offset;

    /* Local variables */

    double att, btt, db, dysa, dysb, symax;
    int lang = -1/* -Wall */;
    int j, l, lplac, nj, nl, nbb, ncase, nclu, numcl, nsylr, ntt;

/* Parameter adjustments */
    --avsyl;

    --srank;
    --syl;
    --negbr;
    --nelem;
    --nsend;
    --ncluv;	--nbest;
    --dys;

    sylinf_dim1 = nsam;
    sylinf_offset = 1 + sylinf_dim1 * 1;
    sylinf -= sylinf_offset;

/*
     construction of clustering vector (ncluv)
     of selected sample (nbest).
*/

    /* Function Body */
    for (l = 1; l <= nsam; ++l) {
	ncase = nbest[l];
	ncluv[l] = (int) (x[(ncase - 1) * jpp] + .1f);
    }

/*     drawing of the silhouettes */

    nsylr = 0;
    *ttsyl = 0.;
    for (numcl = 1; numcl <= kk; ++numcl) {
	ntt = 0;
	for (j = 1; j <= nsam; ++j) {
	    if (ncluv[j] == numcl) {
		++ntt;
		nelem[ntt] = j;
	    }
	}
	for (j = 1; j <= ntt; ++j) {
	    nj = nelem[j];
	    dysb = s * 1.1 + 1.;
	    negbr[j] = -1;

	    for (nclu = 1; nclu <= kk; ++nclu) {
		if (nclu != numcl) {
		    nbb = 0;
		    db = 0.;
		    for (l = 1; l <= nsam; ++l) {
			if (ncluv[l] == nclu) {
			    ++nbb;
			    db += dys[F77_CALL(meet)(&nj, &l)];
			}
		    }
		    btt = (double) nbb;
		    db /= btt;
		    if (db < dysb) {
			dysb = db;
			negbr[j] = nclu;
		    }
		}
	    }

	    if (ntt == 1) {
		syl[j] = 0.;	    continue /* j */;
	    }
	    dysa = 0.;
	    for (l = 1; l <= ntt; ++l) {
		nl = nelem[l];
		dysa += dys[F77_CALL(meet)(&nj, &nl)];
	    }
	    att = (double) (ntt - 1);
	    dysa /= att;
	    if (dysa <= 0.) {
		if (dysb > 0.)
		    syl[j] = 1.;
		else
		    syl[j] = 0.;

		continue /* j */;
	    }

	    if (dysb > 0.) {
		if (dysb > dysa)
		    syl[j] = 1. - dysa / dysb;
		else if (dysb < dysa)
		    syl[j] = dysb / dysa - 1.;
		else /* (dysb == dysa) */
		    syl[j] = 0.;
	    }
	    else {
		syl[j] = -1.;
	    }

	    if (syl[j] < -1.)
		syl[j] = -1.;
	    else if (syl[j] > 1.)
		syl[j] = 1.;

	} /* for(j ..) */

	avsyl[numcl] = 0.;
	for (j = 1; j <= ntt; ++j) {
	    symax = -2.;
	    for (l = 1; l <= ntt; ++l) {
		if (syl[l] > symax) {
		    symax = syl[l];
		    lang = l;
		}
	    }
	    nsend[j] = lang;
	    srank[j] = syl[lang];
	    avsyl[numcl] += srank[j];
	    syl[lang] = -3.;
	}
	*ttsyl += avsyl[numcl];
	avsyl[numcl] /= ntt;

	if (ntt >= 2) {
	    for (l = 1; l <= ntt; ++l) {
		lplac = nsend[l];
		ncase = nelem[lplac];
		++nsylr;
		sylinf[nsylr + sylinf_dim1] = (double) numcl;
		sylinf[nsylr + (sylinf_dim1 << 1)] = (double) negbr[lplac];
		sylinf[nsylr + sylinf_dim1 * 3] = srank[l];
		sylinf[nsylr + (sylinf_dim1 << 2)] = (double) nbest[ncase];
	    }
	}
	else {
	    ncase = nelem[1];
	    ++nsylr;
	    sylinf[nsylr + sylinf_dim1] = (double) numcl;
	    sylinf[nsylr + (sylinf_dim1 << 1)] = (double) negbr[1];
	    sylinf[nsylr + sylinf_dim1 * 3] = 0.;
	    sylinf[nsylr + (sylinf_dim1 << 2)] = (double) nbest[ncase];
	}

    }
    *ttsyl /= (double) (nsam);
    return;
} /* black */

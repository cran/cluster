/*
 * PAM := Partitioning Around Medoids
 *
 * original Id: pam.f,v 1.16 2003/06/03 13:40:56 maechler translated by
 * f2c (version 20031025) and run through f2c-clean,v 1.10 2002/03/28
 */

#include <float.h>

#include <R_ext/Print.h>/* for diagnostics */

#include "cluster.h"
#include "ind_2.h"

/*     carries out a clustering using the k-medoid approach.
 */
void pam(int *nn, int *p, int *kk, double *x, double *dys,
	 int *jdyss, /* jdyss = 0 : compute distances from x
		      *	      = 1 : distances provided	in x */
	 double *valmd, int *jtmd,
	 int *ndyst, int *nsend, int/*logical*/ *nrepr, int *nelem,
	 double *radus, double *damer, double *ttd, double *separ,
	 double *ttsyl, double *obj, int *med, int *ncluv,
	 double *clusinf, double *sylinf, int *nisol)
{
    int clusinf_dim1 = *kk;

    /* Local variables */
    Rboolean all_stats = (obj[0] == 0.),/* if false, only return 'ncluv[]' */
	med_given = (med[0] != 0);/* if true, med[] contain initial medoids */
    int k, i, nhalf, jhalt, trace_lev = (int) obj[1];
    double s, sky;

    /* Function Body */
    if (*jdyss != 1) {
	jhalt = 0;
	F77_CALL(dysta)(nn, p, x, dys, ndyst, jtmd, valmd, &jhalt);
	if (jhalt != 0) {
	    *jdyss = -1; return;
	}
    }
    nhalf = *nn * (*nn - 1) / 2 + 1; /* nhalf := #{distances}+1 = length(dys) */

    /* s := max( dys[.] ), the largest distance */
    for (i = 1, s = 0.; i < nhalf; ++i) /* dys[0] == 0. not used here */
	if (s < dys[i])
	    s = dys[i];

    /* FIXME: work with med[] = (i_1, i_2, ..., i_k)
     * ----- instead nrepr[] = (b_1, ... b_n)   b_i in {0,1} */
    for (i = 0; i < *nn; ++i)
	nrepr[i] = 0;
    if(med_given) { /* if true, med[] contain initial medoids */

	/* for the moment, translate these to nrepr[] 0/1 :
	 * not assuming that the med[] indices are sorted */
	for (k = 0; k < *kk; k++)
	    nrepr[med[k] - 1] = 1;
    }

/*     Build + Swap : */
    bswap(*kk, *nn, nrepr, med_given, trace_lev,
	  radus, damer, ttd, dys, &sky, s, obj);

/*     Compute Clustering & STATs if(all_stats): */
    cstat(kk, nn, nsend, nrepr, all_stats,
	  radus, damer, ttd, separ, &s, dys, ncluv, nelem, med, nisol);
    if(all_stats) {
	for (k = 0; k < *kk; ++k) {
	    clusinf[k]=		(double)       nrepr[k];
	    clusinf[k + clusinf_dim1]	     = radus[k];
	    clusinf[k + (clusinf_dim1 << 1)] = ttd  [k];
	    clusinf[k + clusinf_dim1 * 3]    = damer[k];
	    clusinf[k + (clusinf_dim1 << 2)] = separ[k];
	}
	if (1 < *kk && *kk < *nn) {
	    /* Compute Silhouette info : */
	    dark(*kk, *nn, ncluv, nsend, nelem, nrepr,
		 radus, damer, ttd, ttsyl, dys, &s, sylinf);
	}
    }
} /* pam */

/* NOTE: dysta() is *kept* in Fortran for now --> ./dysta.f */

/* -----------------------------------------------------------

     bswap(): the clustering algorithm in 2 parts:  I. build,	II. swap
*/
void bswap(int kk, int nsam, int *nrepr, Rboolean med_given, int trace_lev,
	   /* nrepr[]: here is boolean (0/1): 1 = "is representative object"  */
	   double *dysma, double *dysmb, double *beter,
	   double *dys, double *sky, double s, double *obj)
{
    int i, j, ij, k;

     /* Parameter adjustments */
    --nrepr;
    --beter;

    --dysma; --dysmb;

    s = s * 1.1 + 1.;/* larger than all dys[];
			replacing by DBL_MAX  changes result - why ? */

/* IDEA: when n (= nsam) is large compared to k (= kk),
 * ----  rather use a "sparse" representation:
 * instead of boolean vector nrepr[] , use  ind_repr <- which(nrepr) !!
 */
    for (i = 1; i <= nsam; ++i)
	dysma[i] = s;

    if(med_given) {
	if(trace_lev)
	    Rprintf("pam()'s bswap(): medoids given\n");

	/* compute dysma[] : dysma[j] = D(j, nearest_representative) */
	for (i = 1; i <= nsam; ++i) {
	    if (nrepr[i] == 1)
		for (j = 1; j <= nsam; ++j) {
		    ij = ind_2(i, j);
		    if (dysma[j] > dys[ij])
			dysma[j] = dys[ij];
		}
	}
    }
    else {

/*  ====== first algorithm: BUILD. ====== */

	if(trace_lev)
	    Rprintf("pam()'s bswap(): build %d medoids:\n", kk);

	/* find  kk  representatives  aka medoids :  */

	for (k = 1; k <= kk; ++k) {

	    /* compute beter[i] for all non-representatives:
	     * also find ammax := max_{..} and nmax := argmax_i{beter[i]} ... */
	    int nmax = -1; /* -Wall */
	    double ammax, cmd;
	    ammax = 0.;
	    for (i = 1; i <= nsam; ++i) {
		if (nrepr[i] == 0) {
		    beter[i] = 0.;
		    for (j = 1; j <= nsam; ++j) {
			cmd = dysma[j] - dys[ind_2(i, j)];
			if (cmd > 0.)
			    beter[i] += cmd;
		    }
		    if (ammax <= beter[i]) {
			/*  does < (instead of <= ) work too? -- NO! */
			ammax = beter[i];
			nmax = i;
		    }
		}
	    }

	    nrepr[nmax] = 1;/* = .true. : found new representative */
	    if (trace_lev >= 2)
		Rprintf("    new repr. %d\n", nmax);

	    /* update dysma[] : dysma[j] = D(j, nearest_representative) */
	    for (j = 1; j <= nsam; ++j) {
		ij = ind_2(nmax, j);
		if (dysma[j] > dys[ij])
		    dysma[j] = dys[ij];
	    }
	}
	/* output of the above loop:  nrepr[], dysma[], ... */
    }

    if(trace_lev) /* >= 2 (?) */ {
	Rprintf("  after build: medoids are");
	for (i = 1; i <= nsam; ++i)
	    if(nrepr[i] == 1) Rprintf(" %2d", i);
	if(trace_lev >= 2) {
	    Rprintf("\n  and min.dist dysma[1:n] are\n");
	    for (i = 1; i <= nsam; ++i) {
		Rprintf(" %6.3g", dysma[i]);
		if(i % 10 == 0) Rprintf("\n");
	    }
	    if(nsam % 10 != 0) Rprintf("\n");
	} else Rprintf("\n");
    }

    *sky = 0.;
    for (j = 1; j <= nsam; ++j)
	*sky += dysma[j];
    obj[0] = *sky / nsam;

    if (kk > 1 || med_given) {

	double small, dz, dzsky;
	int kj, kbest = -1, nbest = -1;/* init: -Wall*/

/* ====== second algorithm: SWAP. ====== */

	/* Hmm: In the following, we RE-compute dysma[];
	 *      don't need it first time; then only need *update* after swap */

/*--   Loop : */
L60:
	for (j = 1; j <= nsam; ++j) {
	    /*  dysma[j] := D_j  d(j, <closest medi>)  [KR p.102, 104]
	     *  dysmb[j] := E_j  d(j, <2-nd cl.medi>)  [p.103] */
	    dysma[j] = s;
	    dysmb[j] = s;
	    for (i = 1; i <= nsam; ++i) {
		if (nrepr[i] == 1) {
		    ij = ind_2(i, j);
		    if (dysma[j] > dys[ij]) {
			dysmb[j] = dysma[j];
			dysma[j] = dys[ij];
		    } else if (dysmb[j] > dys[ij]) {
			dysmb[j] = dys[ij];
		    }
		}
	    }
	}

	dzsky = 1.;
	for (k = 1; k <= nsam; ++k) if (nrepr[k] == 0) {
	    for (i = 1; i <= nsam; ++i) if (nrepr[i] != 0) {
		dz = 0.;
		/* dz := T_{ih} := sum_j C_{jih}  [p.104] : */
		for (j = 1; j <= nsam; ++j) {
		    ij = ind_2(i, j);
		    kj = ind_2(k, j);
		    if (dys[ij] == dysma[j]) {
			small = dysmb[j] > dys[kj]? dys[kj] : dysmb[j];
			dz += (- dysma[j] + small);
		    } else if (dys[kj] < dysma[j])
			dz += (- dysma[j] + dys[kj]);
		}
		if (dzsky > dz) {
		    dzsky = dz;
		    kbest = k;
		    nbest = i;
		}
	    }
	}
	if (dzsky < 0.) { /* found an improving swap */
	    if(trace_lev >= 2)
		Rprintf( "   swp new %d <-> %d old; decreasing diss. by %g\n",
			kbest, nbest, dzsky);
	    nrepr[kbest] = 1;
	    nrepr[nbest] = 0;
	    *sky += dzsky;
	    goto L60;
	}
    }
    obj[1] = *sky / nsam;
} /* bswap */

/* -----------------------------------------------------------
 cstat(): Compute STATistics (numerical output) concerning each partition
*/
void cstat(int *kk, int *nn, int *nsend, int *nrepr, Rboolean all_stats,
	   double *radus, double *damer, double *ttd, double *separ, double *s,
	   double *dys, int *ncluv, int *nelem, int *med, int *nisol)
{
    Rboolean kand;
    int j, k, m = -1, ja, jb, jk, jndz, ksmal = -1/* -Wall */;
    int mevj, njaj, nel, nvn, ntt, nvna, numl, nplac;
    double aja, ajb, dam, djm, dsmal, sep, ttt;

    /* Parameter adjustments */
    --nisol;
    --med;
    --nelem;
    --ncluv;
    --separ;
    --ttd;
    --damer;
    --radus;
    --nrepr;
    --nsend;

    /* nsend[j] := i,  where x[i,] is the medoid to which x[j,] belongs */
    for (j = 1; j <= *nn; ++j) {
	if (nrepr[j] == 0) {
	    dsmal = *s * 1.1f + 1.;
	    for (k = 1; k <= *nn; ++k) {
		if (nrepr[k] == 1) {
		    njaj = ind_2(k, j);
		    if (dsmal > dys[njaj]) {
			dsmal = dys[njaj];
			ksmal = k;
		    }
		}
	    }
	    nsend[j] = ksmal;
	} else {
	    nsend[j] = j;
	}
    }
    /* ncluv[j] := k , the cluster number  (k = 1..*kk) */
    jk = 1;
    nplac = nsend[1];
    for (j = 1; j <= *nn; ++j) {
	ncluv[j] = 0;
	if (nsend[j] == nplac)
	    ncluv[j] = 1;
    }
    for (ja = 2; ja <= *nn; ++ja) {
	nplac = nsend[ja];
	if (ncluv[nplac] == 0) {
	    ++jk;
	    for (j = 2; j <= *nn; ++j) {
		if (nsend[j] == nplac)
		    ncluv[j] = jk;
	    }
	    if (jk == *kk)
		break;
	}
    }

    if(all_stats) { /*	   analysis of the clustering. */

	for (k = 1; k <= *kk; ++k) {
	    ntt = 0;
	    radus[k] = -1.;
	    ttt = 0.;
	    for (j = 1; j <= *nn; ++j) {
		if (ncluv[j] == k) {
		    ++ntt;
		    m = nsend[j];
		    nelem[ntt] = j;
		    djm = dys[ind_2(j, m)];
		    ttt += djm;
		    if (radus[k] < djm)
			radus[k] = djm;
		}
	    }
	    ttd[k] = ttt / ntt;
	    med[k] = m;
	}
	if (*kk == 1) {
	    damer[1] = *s;
	    nrepr[1] = *nn;
	    return;
	}
	/*  ELSE	  kk > 1 : */

	/* numl = number of L-clusters. */
	numl = 0;
	for (k = 1; k <= *kk; ++k) {
	    /*
	      identification of cluster k:
	      nel  = number of objects
	      nelem= vector of object indices */

	    nel = 0;
	    for (j = 1; j <= *nn; ++j) {
		if (ncluv[j] == k) {
		    ++nel;
		    nelem[nel] = j;
		}
	    }
	    nrepr[k] = nel;
	    if (nel == 1) {
		nvn = nelem[1];
		damer[k] = 0.;
		separ[k] = *s * 1.1f + 1.;
		for (j = 1; j <= *nn; ++j) {
		    if (j != nvn) {
			mevj = ind_2(nvn, j);
			if (separ[k] > dys[mevj])
			    separ[k] = dys[mevj];
		    }
		}

		/* Is cluster k
		   1) an L-cluster	 or
		   2) an L*-cluster ? */
		if (separ[k] == 0.)
		    ++numl;

	    }
	    else { /*	       nel != 1 : */
		dam = -1.;
		sep = *s * 1.1f + 1.;
		kand = TRUE;
		for (ja = 1; ja <= nel; ++ja) {
		    nvna = nelem[ja];
		    aja = -1.;
		    ajb = *s * 1.1f + 1.;
		    for (jb = 1; jb <= *nn; ++jb) {
			jndz = ind_2(nvna, jb);
			if (ncluv[jb] == k) {
			    if (aja < dys[jndz])
				aja = dys[jndz];
			} else {
			    if (ajb > dys[jndz])
				ajb = dys[jndz];
			}
		    }
		    if (kand && aja >= ajb)
			kand = FALSE;
		    if (dam < aja)
			dam = aja;
		    if (sep > ajb)
			sep = ajb;
		}
		separ[k] = sep;
		damer[k] = dam;
		if (kand) {
		    ++numl;
		    if (dam >= sep) /*	L-cluster */
			nisol[k] = 1;
		    else/*		L*-cluster */
			nisol[k] = 2;
		    continue /* k */;
		}
	    }
	    /* nel = 1 or (!kand) : */
	    nisol[k] = 0;

	}/* for(k) */

    } /* all_stats */

} /* cstat */

/* -----------------------------------------------------------
     Compute Silhouette Information :
 */
void dark(int kk, int nn, int *ncluv,
	  int *nsend, int *nelem, int *negbr,
	  double *syl, double *srank, double *avsyl, double *ttsyl,
	  double *dys, double *s, double *sylinf)
{
    int sylinf_d = nn; /* sylinf[nn, 4] */
    int j, k, l, lang=-1 /*Wall*/, lplac, k_, nj, nl, nbb, ntt, nsylr;
    double db, dysa, dysb, symax;
    /* pointers to sylinf[] columns:*/
    double *sylinf_2, *sylinf_3, *sylinf_4;
    sylinf_2 = sylinf	+ sylinf_d;
    sylinf_3 = sylinf_2 + sylinf_d;
    sylinf_4 = sylinf_3 + sylinf_d;

    /* Parameter adjustments */
    --avsyl;
    --ncluv;

    nsylr = 0;
    *ttsyl = 0.;
    for (k = 1; k <= kk; ++k) {

	/* nelem[0:(ntt-1)] := indices (1-based) of obs. in cluster k : */
	ntt = 0;
	for (j = 1; j <= nn; ++j) {
	    if (ncluv[j] == k) {
		nelem[ntt] = j;
		++ntt;
	    }
	}

	for (j = 0; j < ntt; ++j) {/* (j+1)-th obs. in cluster k */
	    nj = nelem[j];
	    dysb = *s * 1.1 + 1.;
	    negbr[j] = -1;
	    /* for all clusters  k_ != k : */
	    for (k_ = 1; k_ <= kk; ++k_) if (k_ != k) {
		nbb = 0;
		db = 0.;
		for (l = 1; l <= nn; ++l) if (ncluv[l] == k_) {
		    ++nbb;
		    if (l != nj)
			db += dys[ind_2(nj, l)];
		}
		db /= nbb; /* now  db(k_) := mean( d[j, l]; l in C_{k_} ) */
		if (dysb > db) {
		    dysb = db;
		    negbr[j] = k_;
		}
	    }/* negbr[j] := arg max_{k_} db(k_) */
	    if (ntt > 1) {
		dysa = 0.;
		for (l = 0; l < ntt; ++l) {
		    nl = nelem[l];
		    if (nj != nl)
			dysa += dys[ind_2(nj, nl)];
		}
		dysa /= ntt - 1;
		if (dysa > 0.) {
		    if (dysb > 0.) {
			if (dysb > dysa)
			    syl[j] = 1. - dysa / dysb;
			else if (dysb < dysa)
			    syl[j] = dysb / dysa - 1.;
			else /* dysb == dysa: */
			    syl[j] = 0.;

			if (syl[j] < -1.)
			    syl[j] = -1.;
			else if (syl[j] > 1.)
			    syl[j] = 1.;

		    } else {
			syl[j] = -1.;
		    }
		}
		else /* dysa == 0 */ if (dysb > 0.)
		    syl[j] = 1.;
		else
		    syl[j] = 0.;
	    }
	    else { /*	  ntt == 1: */
		syl[j] = 0.;
	    }
	} /* for( j ) */
	avsyl[k] = 0.;
	if (ntt == 0) /* this can happen when medoids are user-specified !*/
	    continue; /* next k */

	for (j = 0; j < ntt; ++j) {
	    symax = -2.;
	    for (l = 0; l < ntt; ++l) {
		if (symax < syl[l]) {
		    symax = syl[l];
		    lang = l;
		}
	    }
	    nsend[j] = lang;
	    srank[j] = symax; /* = syl[lang] */
	    avsyl[k] += srank[j];
	    syl[lang] = -3.;
	}
	*ttsyl += avsyl[k];
	avsyl[k] /= ntt;
	if (ntt == 1) {
	    sylinf  [nsylr] = (double) k;
	    sylinf_2[nsylr] = (double) negbr[0];
	    sylinf_3[nsylr] = 0.;
	    sylinf_4[nsylr] = (double) nelem[0];
	    ++nsylr;
	} else {
	    for (j = 0; j < ntt; ++j) {
		lplac = nsend[j];
		sylinf	[nsylr] = (double) k;
		sylinf_2[nsylr] = (double) negbr[lplac];
		sylinf_3[nsylr] = srank[j];
		sylinf_4[nsylr] = (double) nelem[lplac];
		++nsylr;
	    }
	}
    } /* for (k) */
    *ttsyl /= nn;
} /* dark */

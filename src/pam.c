/* -*- mode: c; kept-new-versions: 25; kept-old-versions: 20 -*- */

/*
 * PAM := Partitioning Around Medoids
 *
 * $Id$
 * original Id: pam.f,v 1.16 2003/06/03 13:40:56 maechler translated by
 * f2c (version 20031025) and run through f2c-clean,v 1.10 2002/03/28
 */

#include <R_ext/Print.h>/* for diagnostics */

#include "cluster.h"

/*     carries out a clustering using the k-medoid approach.
 */
void pam(int *nn, int *jpp, int *kk, double *x, double *dys,
	 int *jdyss, /* jdyss = 0 : compute distances from x
		      *	      = 1 : distances provided	in x */
	 double *valmd, int *jtmd,
	 int *ndyst, int *nsend, int *nrepr, int *nelem,
	 double *radus, double *damer, double *ttd, double *separ,
	 double *ttsyl, double *obj, int *med, int *ncluv,
	 double *clusinf, double *sylinf, int *nisol)
{
    int clusinf_dim1 = *kk;

    /* Local variables */
    Rboolean all_stats = (obj[0] == 0.);/* if false, only return 'ncluv[]' */
    int k, l, nhalf, jhalt;
    double s, sky;

    /* Function Body */
    if (*jdyss != 1) {
	jhalt = 0;
	F77_CALL(dysta)(nn, jpp, x, dys, ndyst,
			jtmd, valmd, &jhalt);
	if (jhalt != 0) {
	    *jdyss = -1; return;
	}
    }
    nhalf = *nn * (*nn - 1) / 2 + 1; /* nhalf := #{distances}+1 = length(dys) */

    /* s := max( dys[.] ), the largest distance */
    for (l = 1, s = 0.; l < nhalf; ++l) /* dys[0] == 0. is NOT used !*/
	if (s < dys[l])
	    s = dys[l];

/*     Build + Swap : */
    bswap(kk, nn, nrepr, radus, damer, ttd, dys, &sky, &s, obj);

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
	    dark(kk, nn, &nhalf, ncluv, nsend, nelem, nrepr,
		 radus, damer, ttd, ttsyl, dys, &s, sylinf);
	}
    }
} /* pam */

/* NOTE: dysta() is *kept* in Fortran for now --> ./dysta.f */

/* -----------------------------------------------------------

     bswap(): the clustering algorithm in 2 parts:  I. build,	II. swap
*/
void bswap(int *kk, int *nn, int *nrepr,
	   /* nrepr[]: here is boolean (0/1): 1 = "is representative object"  */
	   double *dysma, double *dysmb, double *beter,
	   double *dys, double *sky, double *s, double *obj)
{
    int i, j, ij, k, kj, kbest = -1, nbest = -1, njn, nmax = -1;/* init: -Wall*/
    double ammax, small, cmd, dz, dzsky;

     /* Parameter adjustments */
    --beter;
    --dysmb;
    --dysma;
    --nrepr;


/*     first algorithm: build. */

    for (i = 1; i <= *nn; ++i) {
	nrepr[i] = 0;
	dysma[i] = *s * 1.1f + 1.;
    }
    for (k = 1; k <= *kk; ++k) {
	for (i = 1; i <= *nn; ++i) {
	    if (nrepr[i] == 0) {
		beter[i] = 0.;
		for (j = 1; j <= *nn; ++j) {
		    cmd = dysma[j] - dys[ind_2(i, j)];
		    if (cmd > 0.)
			beter[i] += cmd;
		}
	    }
	}
	ammax = 0.;
	for (i = 1; i <= *nn; ++i) {
	    if (nrepr[i] == 0 && ammax <= beter[i]) {
		/*  does < (instead of <= ) work too? -- NO! */
		ammax = beter[i];
		nmax = i;
	    }
	}
	nrepr[nmax] = 1;/* = .true. : *is* a representative */
	for (j = 1; j <= *nn; ++j) {
	    njn = ind_2(nmax, j);
	    if (dysma[j] > dys[njn])
		dysma[j] = dys[njn];
	}
    }
    *sky = 0.;
    for (j = 1; j <= *nn; ++j)
	*sky += dysma[j];
    obj[0] = *sky / *nn;

    if (*kk > 1) {

	/*     second algorithm: swap. */

/*--   Loop : */
L60:
	for (j = 1; j <= *nn; ++j) {
	    dysma[j] = *s * 1.1f + 1.;
	    dysmb[j] = *s * 1.1f + 1.;
	    for (i = 1; i <= *nn; ++i) {
		if (nrepr[i] == 1) {
		    ij = ind_2(i, j);
		    if (dys[ij] < dysma[j]) {
			dysmb[j] = dysma[j];
			dysma[j] = dys[ij];
		    } else {
			if (dysmb[j] > dys[ij])
			    dysmb[j] = dys[ij];
		    }
		}
	    }
	}
	dzsky = 1.;
	for (k = 1; k <= *nn; ++k) if (nrepr[k] == 0) {
	    for (i = 1; i <= *nn; ++i) if (nrepr[i] == 1) {
		dz = 0.;
		for (j = 1; j <= *nn; ++j) {
		    ij = ind_2(i, j);
		    kj = ind_2(k, j);
		    if (dys[ij] == dysma[j]) {
			small = dysmb[j];
			if (small > dys[kj])
			    small = dys[kj];
			dz += (- dysma[j] + small);
		    } else {
			if (dys[kj] < dysma[j])
			    dz += (- dysma[j] + dys[kj]);
		    }
		}
		if (dzsky > dz) {
		    dzsky = dz;
		    kbest = k;
		    nbest = i;
		}
	    }
	}
	if (dzsky < 0.) {
	    nrepr[kbest] = 1;
	    nrepr[nbest] = 0;
	    *sky += dzsky;
	    goto L60;
	}
    }
    obj[1] = *sky / *nn;
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
void dark(int *kk, int *nn, int *hh, int *ncluv,
	  int *nsend, int *nelem, int *negbr,
	  double *syl, double *srank, double *avsyl, double *ttsyl,
	  double *dys, double *s, double *sylinf)
{
    int sylinf_dim1 = *nn; /* sylinf[*nn, 4] */
    int j, k, l, lang=-1 /*Wall*/, lplac, k_, nj, nl, nbb, ntt, nsylr;
    double db, dysa, dysb, symax;
    /* pointers to sylinf[] columns:*/
    double *sylinf_2, *sylinf_3, *sylinf_4;
    sylinf_2 = sylinf	+ sylinf_dim1;
    sylinf_3 = sylinf_2 + sylinf_dim1;
    sylinf_4 = sylinf_3 + sylinf_dim1;
    /* Parameter adjustments */
    --avsyl;
    --ncluv;

    nsylr = 0;
    *ttsyl = 0.;
    for (k = 1; k <= *kk; ++k) {
	ntt = 0;
	for (j = 1; j <= *nn; ++j) {
	    if (ncluv[j] == k) {
		nelem[ntt] = j;
		++ntt;
	    }
	}
	for (j = 0; j < ntt; ++j) {/* (j+1)-th obs. in cluster k */
	    nj = nelem[j];
	    dysb = *s * 1.1f + 1.;
	    negbr[j] = -1;
	    for (k_ = 1; k_ <= *kk; ++k_) if (k_ != k) {
		nbb = 0;
		db = 0.;
		for (l = 1; l <= *nn; ++l) if (ncluv[l] == k_) {
		    ++nbb;
		    if (l != nj)
			db += dys[ind_2(nj, l)];
		}
		db /= nbb;
		if (dysb > db) {
		    dysb = db;
		    negbr[j] = k_;
		}
	    }
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
	if (ntt < 2) {
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
    *ttsyl /= *nn;
} /* dark */

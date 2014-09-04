/* Produced by
 * $Id: f2c-clean,v 1.10 2002/03/28 16:37:27 maechler Exp $
 *
 * twins.f -- translated by f2c (version 20031025).
 */

#include <math.h>

#include <Rmath.h>
#include <R_ext/Print.h>/* for diagnostics */
#include <R_ext/Utils.h>/* for interrupting */

#include "cluster.h"
#include "ind_2.h"

// the auxiliary routines
static void agnes(int nn, int *kwan, int *ner, double *ban, double dys[],
		  int method, double *alpha, int *merge, int trace_lev);
static void splyt(int nn, int *kwan, int *ner, double *ban, double dys[],
		  int *merge, int trace_lev);
static double min_dis(double dys[], int ka, int kb, int ner[]);

/*     This program performs agglomerative nesting (AGNES) using the */
/*     group average method (_or others_) of Sokal and Michener (1958), */
/*     as well as divisive analysis (DIANA) using the method of */
/*     Mcnaughton-Smith, Williams, Dale, and Mockett (1964). */

/*     Extended by Martin Maechler to allow the (flexible) */
/*     Lance-Williams clustering method (with parameters alpha[1:4]) */

void twins(int *nn, // = maximal number of objects
	   int *jpp,// = maximal number of variables used in the analysis
	   double *x,
	   double *dys,
	   double *dys2,// dys2(.) can have length 1, if(!keep.diss)
	   int *jdyss, /* jdyss (in/out): initially, jdyss mod 10 = 1 : <==> diss = TRUE
			* jdyss < 10 : don't save dissimilarities (C.keep.diss = FALSE) */
	   double *valmd,
	   int *jtmd, int *ndyst, int *jalg,
	   int *method,  // for agnes()  only
	   int *kwan,
	   int *ner,     // = order []   (in R)
	   double *ban,  // = height[]
	   double *coef,
	   double *alpha, int *merge, int *trace_lev)
{
    if (*jdyss % 10 == 1) {
	*jpp = 1;
    } else { // compute distances
	int jhalt = 0;
	F77_CALL(dysta)(nn, jpp, x, dys, ndyst, jtmd, valmd, &jhalt);
	/*       ------ in ./dysta.f */
	if (jhalt != 0) { *jdyss = -1; return; }
    }
    if (*jdyss >= 10) { /*        save distances for S */
	Memcpy(dys2, dys, (*nn * (*nn - 1) / 2 + 1));
    }
    if (*jalg != 2) {
	// AGNES
	agnes(*nn, kwan, ner, ban, dys, *method, alpha, merge, trace_lev[0]);
    } else {
	// DIANA
	splyt(*nn, kwan, ner, ban, dys,                 merge, trace_lev[0]);
    }
    // Compute agglomerative/divisive coefficient from banner:
    *coef = bncoef(*nn, ban);
    return;
} /* twins */

/*     ----------------------------------------------------------- */
/*     AGNES agglomeration */
static void
agnes(int nn, int *kwan, int *ner, double *ban,
      double dys[], int method, double *alpha, int *merge, int trace_lev)
{

/* VARs */
    double akb, d_min;
    int j, k, l1, l2, lq,
	nclu, lnum, lput, lenda, lendb, lnext, n_1 = nn - 1,
	la = -1, lb = -1, llast = -1, lfyrs = -1, // <- against (unnecessary) warnings [-Wall]
	_d, nmerge;
    Rboolean has_a3 = FALSE, has_a4 = FALSE;// is alpha[3] or [4] == 0 -- for Lance-Williams

    /* System generated locals */
    int merge_dim1 = n_1;
    int merge_offset = 1 + merge_dim1;

    /* Parameter adjustments */
    merge -= merge_offset;
    --ban;
    --ner;
    --kwan;
    --alpha;

    if(trace_lev) {
	_d = (nn >= 100) ? 3 : (nn >= 10) ? 2 : 1;
	Rprintf("C agnes(n=%*d, method = %d, ..): ", _d,nn, method);
    } else _d = -1;// -Wall

    if(method == 6 || method == 7) {
	// 6: "flexible": "Flexible Strategy" (K+R p.236 f) extended to 'Lance-Williams'
	// 7: "gaverage" aka Flexible UPGMA (Belbin et al., 1992)
	has_a3 = (alpha[3] != 0.);
	has_a4 = (alpha[4] != 0.);
	if(trace_lev) {
	    if(has_a4)
		Rprintf("|par| = 4, alpha[1:4] = (%g,%g,%g,%g); ",
			alpha[1],alpha[2],alpha[3],alpha[4]);
	    else if(has_a3)
		Rprintf("|par| = 3, alpha[1:3] = (%g,%g,%g); ",
			alpha[1],alpha[2],alpha[3]);
	}
    }

//  Starting with nn clusters, kwan[j] := #{obj} in cluster j
    for (j = 1; j <= nn; ++j) {
	kwan[j] = 1;
	ner[j] = j;
    }

// ----------------------------------------------------------------------------
/*     find closest clusters */
    nmerge = 1;
    if(trace_lev) Rprintf("%d merging steps\n", n_1);
    for (nclu = n_1; nclu >= 1; --nclu) {
	// j := min_j { kwan[j] > 0} = first non-empty cluster
	j = 1; do { j++; } while(kwan[j] == 0);
	if(trace_lev >= 2) Rprintf(" nclu=%*d, j=%*d, ", _d,nclu, _d,j);

	d_min = dys[ind_2(1, j)] * 1.1f + 1.;
	for (k = 1; k <= n_1; ++k) if (kwan[k] > 0) {
		for (j = k + 1; j <= nn; ++j) if (kwan[j] > 0) {
			int k_j = ind_2(k, j);
			if (d_min >= dys[k_j]) { // Note: also when "==" !
			    d_min =  dys[k_j];
			    la = k;
			    lb = j;
			}
		    }
	    }
	// --> closest clusters  {la < lb}  are at distance  d_min
	if(trace_lev >= 2) Rprintf("d_min = D(%*d,%*d) = %#g; ", _d,la, _d,lb, d_min);

/*     merge-structure for plotting tree in S */

	l1 = -la;
	l2 = -lb;
	for (j = 1; j <= (nmerge - 1); ++j) {
	    if (merge[j + merge_dim1] == l1 || merge[j + (merge_dim1 << 1)] == l1)  l1 = j;
	    if (merge[j + merge_dim1] == l2 || merge[j + (merge_dim1 << 1)] == l2)  l2 = j;
	}
	merge[nmerge + merge_dim1] = l1;
	merge[nmerge + (merge_dim1 << 1)] = l2;
	++nmerge;

/*     determine lfyrs and llast */

	for (k = 1; k <= nn; ++k) {
	    if (ner[k] == la) lfyrs = k;
	    if (ner[k] == lb) llast = k;
	}
	ban[llast] = d_min;

	if(trace_lev >= 2) Rprintf("nmerge = %*d, last=%*d;", _d,nmerge, _d,llast);

/*     if the two clusters are next to each other, ner must not be changed */

	lnext = lfyrs + kwan[la];
	if (lnext != llast) { /*     updating ner and ban */
	    if(trace_lev >= 2) Rprintf(" upd(n,b);");
	    lput = lfyrs + kwan[la];
	    lnum = llast - lput;
	    for (k = 1; k <= lnum; ++k) {
		int lka = ner[lput];
		akb = ban[lput];
		lenda = llast + kwan[lb] - 2;
		lendb = lenda + 1;
		for (j = lput; j <= lenda; ++j) {
		    ner[j] = ner[j + 1];
		    ban[j] = ban[j + 1];
		}
		ner[lendb] = lka;
		ban[lendb] = akb;
	    }
	}
	if(trace_lev >= 3) Rprintf("\n");

/*     We will merge A & B into  A_{new} */

	// Calculate new dissimilarities d(q, A_{new})
	for (lq = 1; lq <= nn; ++lq) { //  for each other cluster 'q'

	    if (lq == la || lq == lb || kwan[lq] == 0)
		continue;

	    int naq = ind_2(la, lq);
	    int nbq = ind_2(lb, lq);
	    if(trace_lev >= 3)
		Rprintf(" old D(A, j), D(B, j), j=%*d  = (%g,%g); ",
			_d,lq, dys[naq], dys[nbq]);

	    switch(method) {
	    case 1: { //   1: unweighted pair-]group average method, UPGMA
		double
		    ta = (double) kwan[la],
		    tb = (double) kwan[lb],
		    fa = ta / (ta + tb),
		    fb = tb / (ta + tb);
		dys[naq] = fa * dys[naq] + fb * dys[nbq];
		break;
	    }
	    case 2: /*     2: single linkage */
		dys[naq] = fmin2(dys[naq], dys[nbq]);
		break;
	    case 3: /*     3: complete linkage */
		dys[naq] = fmax2(dys[naq], dys[nbq]);
		break;
            case 4: { //   4: ward's method
		double
		    ta = (double) kwan[la],
		    tb = (double) kwan[lb],
		    tq = (double) kwan[lq],
		    fa = (ta + tq) / (ta + tb + tq),
		    fb = (tb + tq) / (ta + tb + tq),
		    fc = -tq / (ta + tb + tq);
		int nab = ind_2(la, lb);
		dys[naq] = sqrt(fa * dys[naq] * dys[naq] +
				fb * dys[nbq] * dys[nbq] +
				fc * dys[nab] * dys[nab]);
		break;
	    }
	    case 5: /*     5: weighted average linkage */
		dys[naq] = (dys[naq] + dys[nbq]) / 2.;
		break;
	    case 6: { //   6: "Flexible Strategy" (K+R p.236 f) extended to 'Lance-Williams'
		double dNew = alpha[1] * dys[naq] + alpha[2] * dys[nbq];
		if(has_a3) dNew += alpha[3] * dys[ind_2(la, lb)];
		if(has_a4) dNew += alpha[4] * fabs(dys[naq] - dys[nbq]);
		dys[naq] = dNew;
		/* Lance-Williams would allow alpha(1:2) to *depend* on |cluster|
		 * could also include the extensions of Jambu(1978) --
		 * See Gordon A.D. (1999) "Classification" (2nd ed.) p.78 ff */
		break;
	    }
	    case 7: {/*    7: generalized "average" = Flexible UPGMA (Belbin et al., 1992)
		      * Applies the flexible Lance-Williams formula to the UPGMA, aka
		      * "average" case 1 above, i.e., alpha_{1,2} depend on cluster sizes: */
		double
		    ta = (double) kwan[la],
		    tb = (double) kwan[lb],
		    fa = alpha[1] * ta / (ta + tb),
		    fb = alpha[2] * tb / (ta + tb),
		    dNew = fa * dys[naq] + fb * dys[nbq];
		if(has_a3) dNew += alpha[3] * dys[ind_2(la, lb)];
		if(has_a4) dNew += alpha[4] * fabs(dys[naq] - dys[nbq]);
		dys[naq] = dNew;
		break;
	    }
	    default:
		error(_("invalid method (code %d)"), method);
	    }
	    if(trace_lev >= 3)
		Rprintf(" new D(A', %*d) = %g\n", _d,lq, dys[naq]);
	} // for (lq ..)
	kwan[la] += kwan[lb];
	kwan[lb] = 0;
	if(trace_lev >= 2)
	    Rprintf("%s new size = %d\n", (trace_lev >= 3)? " --> " : "", kwan[la]);

    }// for(nclu ..)
    return;
} /* agnes */

/*     ----------------------------------------------------------- */

/*       cf = ac := "Agglomerative Coefficient" from AGNES banner */
/*  or   cf = dc := "Divisive Coefficient"      from DIANA banner */

void R_bncoef(int *nn, double *ban, double *cf) {
    *cf = bncoef(*nn, ban);
}

double bncoef(int n, double *ban)
{
    int k, n_1 = n-1;

    double sup = 0.;// sup := max_k ban[k]
    for(k = 1; k < n; ++k)
	if (sup < ban[k])
	    sup = ban[k];

    double cf = 0.;
    for (k = 0; k < n; ) {
	int kearl = (k > 0) ? k : 1,
	    kafte = (++k < n) ? k : n_1;
	double syze = fmin2(ban[kearl], ban[kafte]);
	cf += (1. - syze / sup);
    }
    return cf / n;
} /* bncoef */

/*     ----------------------------------------------------------- */
/*     DIANA "splitting" */

static void
splyt(int nn, int *kwan, int *ner, double *ban,
      double dys[], int *merge, int trace_lev)
{
    /* Local variables */
    int j, ja, jb, k, l;
    int jma, jmb, lmm, llq, lmz,
	lxx, lxy, lmma, lmmb, lner, nclu;
    int lchan, nhalf, nmerge, n_1 = nn - 1, splyn;
    /* against (unnecessary) warnings [-Wall]: */
    int jaway = -1, lndsd = -1;

    double da, db, cs, sd, dyff;

    /* System generated locals */
    int merge_dim1 = n_1;
    int merge_offset = 1 + merge_dim1;

    /* Parameter adjustments */
    merge -= merge_offset;
    --ban;
    --ner;
    --kwan;


    /*     initialization */
    nclu = 1;
    nhalf = nn * n_1 / 2 + 1;
    for (l = 1; l <= nn; ++l) {
	kwan[l] = 0;
	ban[l] = 0.;
	ner[l] = l;
    }
    kwan[1] = nn;
    ja = 1;

/*     cs :=  diameter of data set */

    cs = 0.;
    k = 0;
    do {
	if (cs < dys[k])
	    cs = dys[k];
	++k;
    } while (k < nhalf);
    if(trace_lev)
	Rprintf("C diana(): ndist= %d, diameter = %g\n", nhalf, cs);

/*     prepare for splitting */

//____________ Big Loop _________________________________________________
L30:
    jb = ja + kwan[ja] - 1;
    jma = jb;

/*     special case of a pair of objects */

    if (kwan[ja] == 2) {
	kwan[ja] = 1;
	kwan[jb] = 1;
	ban [jb] = dys[ind_2(ner[ja], ner[jb])];
    }
    else {
	/*     finding first object to be shifted */
	double bygsd = -1.;
	for (l = ja; l <= jb; ++l) {
	    lner = ner[l];
	    sd = 0.;
	    for (j = ja; j <= jb; ++j)
		sd += dys[ind_2(lner, ner[j])];
	    if (bygsd < sd) {
		bygsd = sd;
		lndsd = l;
	    }
	}

/*     shifting the first object */
	--kwan[ja];
	kwan[jb] = 1;
	if (jb != lndsd) {
	    lchan = ner[lndsd];
	    lmm = jb - 1;
	    for (lmma = lndsd; lmma <= lmm; ++lmma) {
		lmmb = lmma + 1;
		ner[lmma] = ner[lmmb];
	    }
	    ner[jb] = lchan;
	}
	splyn = 0;
	jma = jb - 1;

/*     finding the next object to be shifted */

	do {
	    splyn++;
	    int rest = (jma - ja);
	    double bdyff = -1.;
	    for (l = ja; l <= jma; ++l) {
		lner = ner[l];
		da = 0.;
		for (j = ja; j <= jma; ++j)
		    da += dys[ind_2(lner, ner[j])];
		da /= rest;
		db = 0.;
		for (j = jma + 1; j <= jb; ++j)
		    db += dys[ind_2(lner, ner[j])];
		db /= splyn;
		dyff = da - db;
		if (bdyff < dyff) {
		    bdyff = dyff;
		    jaway = l;
		}
	    } /* end for(l ..) */
	    jmb = jma + 1;

/*     shifting the next object when necessary */

	    if (bdyff <= 0.)
		break; // out of  "object shifting"  while(.) loop

	    if (jma != jaway) {
		lchan = ner[jaway];
		lmz = jma - 1;
		for (lxx = jaway; lxx <= lmz; ++lxx)
		    ner[lxx] = ner[lxx + 1];
		ner[jma] = lchan;
	    }
	    for (lxx = jmb; lxx <= jb; ++lxx) {
		lxy = lxx - 1;
		if (ner[lxy] < ner[lxx])
		    break;
		lchan = ner[lxy]; ner[lxy] = ner[lxx]; ner[lxx] = lchan;
	    }

	    --kwan[ja];
	    kwan[jma] = kwan[jmb] + 1;
	    kwan[jmb] = 0;
	    --jma;
	    jmb = jma + 1;

	} while (jma != ja);


// 200:     switch the two parts when necessary

	if (ner[ja] >= ner[jmb]) {
	    int lxxa = ja;
	    for (int lgrb = jmb; lgrb <= jb; ++lgrb) {
		++lxxa;
		lchan = ner[lgrb];
		int lxg = -1;
		for (lxy = lxxa; lxy <= lgrb; ++lxy) {
		    int lxf = lgrb - lxy + lxxa;
		    lxg = lxf - 1;
		    ner[lxf] = ner[lxg];
		}
		ner[lxg] = lchan;
	    }
	    llq = kwan[jmb];
	    kwan[jmb] = 0;
	    jma = ja + jb - jma - 1;
	    jmb = jma + 1;
	    kwan[jmb] = kwan[ja];
	    kwan[ja] = llq;
	}

/* 300 :    compute level for banner */

	if (nclu == 1) {
	    ban[jmb] = cs;
	} else {
	    ban[jmb] = min_dis(dys, ja, jb, &ner[1]);
	}

    }

    if (++nclu < nn) { /* continue splitting until all objects are separated */
	if (jb != nn) {
 L420:
	    ja += kwan[ja];
	    if (ja <= nn) {
		if (kwan[ja] <= 1)
		    goto L420;
		else
		    goto L30;
	    }
	}
	ja = 1;
	if (kwan[ja] == 1)
	    goto L420;
	else
	    goto L30;
    }
//____________ End Big Loop _________________________________________________

/* 500 :  merge-structure for plotting tree in S */

    for (nmerge = 1; nmerge <= n_1; ++nmerge) {
	int nj = -1, l1, l2;
	double dmin = cs;
	for (j = 2; j <= nn; ++j) {
	    if (kwan[j] >= 0 && dmin >= ban[j]) {
		dmin = ban[j];
		nj = j;
	    }
	}
	kwan[nj] = -1;
	l1 = -ner[nj - 1];
	l2 = -ner[nj];
	for (j = 1; j <= (nmerge - 1); ++j) {
	    if (merge[j + merge_dim1] == l1 || merge[j + (merge_dim1 << 1)] == l1)  l1 = j;
	    if (merge[j + merge_dim1] == l2 || merge[j + (merge_dim1 << 1)] == l2)  l2 = j;
	}
	merge[nmerge + merge_dim1] = l1;
	merge[nmerge + (merge_dim1 << 1)] = l2;
    }
    return;
} /* splyt */

/*     ----------------------------------------------------------- */
/*     used in splyt() above */
static double min_dis(double dys[], int ka, int kb, int ner[])
{
    double dm = 0.;
    for(int k = ka -1; k < kb -1; ++k) {
	int ner_k = ner[k];
	for (int j = k+1; j < kb; ++j) {
	    int k_j = ind_2(ner_k, ner[j]);
	    if (dm < dys[k_j])
		dm = dys[k_j];
	}
    }
    return dm;
} /* min_dis */


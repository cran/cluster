/* inlined, to be included in pam.c and clara.c */

static
#ifdef __GNUC__
__inline__
#endif
int ind_2(int l, int j)
{
/* Utility, originally FORTRAN,	 called "meet"; called from CLARA, PAM & TWINS.
 * Original code had MEET(), MEET2(), and MEET3() in the 3 source files.

 * ind_2(l,j) returns the *index* of dys() where diss. d(l,j) is stored:
 *	  d(l,j) == dys[ind_2(l,j)]
 *
 * MM: changed to work with 0-origin matrices  dys[], but l,j are >= 1
 */
#ifdef was_orig
    if(l > j)
	return (l-2)*(l-1)/2 + j;
    else if(l == j)
	return 0;/* and the first element, dys[0] is := 0.  permanently! */
    else /* l < j */
	return (j-2)*(j-1)/2 + l;
#else
    /* from Li Long -- optimizes particularly well on Itanium2 */
    int result = 0;
    int m = l>j ? l : j;
    int n = l>j ? j : l;
    if (l != j)
	result = (m-2)*(m-1)/2 + n;
    return result;
#endif
}


#include <R.h>

void sildist(double * d, int *n, int * clustering, int *k,
  double * diC, double *bi, int * counts, double * ai,
  double * si, int * neighbor, int * ismat){
// d          : distance matrix : in dist format, length : n*(n-1)/2
// n          : number of Subjets (attr(d,'Size'))
// clustering : clustering vector
// k          : number of clusters
// diC        : diC
// bi         : min_C diC
// ai         : ai
// si         : (ai-bi)/max(ai,bi)
// neighbor   : neighbor
// ismat      : boolean : is d a matrix (1) or a dist vector (0)

    int i,j,l ;
    int nm1 = *n-1 ;
    int ci,cj;

    int indix;

    for(i=0, l=0; i<nm1; i++){
        ci = clustering[i] - 1 ;

        if(*ismat) l = n*i+i+1 ;

        for(j=i+1 ; j<*n; j++,l++){

            cj = clustering[j] - 1 ;


                indix = (*k)*i+cj ;
                diC[indix]+=d[l] ;
                counts[indix]++;

                indix = (*k)*j+ci ;
                diC[indix]+=d[l] ;
                counts[indix]++;

        }
    }

    for(i=0; i<*n; i++){
        for(j=0; j<*k; j++){
            indix = *k*i+j;
            diC[indix]/=counts[indix];
        }

        // bi = min_C diC

        ci = clustering[i]-1 ;

        ai[i] = diC[*k*i+ci]  ;


        if(ci==0){
            bi[i] = diC[*k*i+1] ;
            neighbor[i] = 2;
        }
        else{
            bi[i] = diC[*k*i] ;
            neighbor[i] = 1 ;
        }

        for(j=1; j<*k;j++){
            if(j!=ci){
                indix = *k*i+j ;
              if(bi[i]>diC[indix]){
                  bi[i] = diC[indix] ;
                  neighbor[i] = j+1 ;
              }
            }
        }

        si[i] = ai[i];
        if(bi[i]>ai[i]) si[i] = bi[i] ;
        si[i] = (bi[i]-ai[i])/si[i] ;

    }


}

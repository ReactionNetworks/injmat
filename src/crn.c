/* Copyright (C) 2010-2012, Murad Banaji
 *
 * This file is part of QUALMAT
 *
 * QUALMAT is free software; you can redistribute it and/or 
 * modify it under the terms of the GNU General Public License 
 * as published by the Free Software Foundation; either version 2, 
 * or (at your option) any later version.
 *
 * QUALMAT is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with QUALMAT: see the file COPYING.  If not, write to 
 * the Free Software Foundation, Inc., 59 Temple Place - Suite 330, 
 * Boston, MA 02111-1307, USA. 

 * For comments, queries, or bug-reports about QUALMAT, please 
 * email qualmat@gmail.com

 */

#include "crn.h"

//#include <iostream>

/* note that these are not meant to be the most efficient */
/* algorithms. Also there is very little error checking  */
/* (to maximise speed), and if the routines are used carelessly,  */
/* it is possible to leak memory.  */

/* Note that in the case of large sparse matrices, the algorithms */
/* presented here serve to highlight the need for graph-theoretic */
/* algorithms */

/* A permutation of [1, ..., n] is just an n-vector */

/* calculate the unsigned term in an n X n matrix corresponding */
/* to permutation tm */

int unsterm(int **imat, int *tm, int n){
  int i,tot=1;
  for(i=0;i<n;i++)
    tot*=imat[i][tm[i]];
  return tot;
}

/* the unsigned term in the determinant of a k X k submatrix  */
/* of n X m matrix imat, corresponding to permutation tm */

int unsterm1(int **imat, int n, int m, int *vec1, int *tm, int k){
  int i,j,tot=1;
  for(i=0;i<k;i++){
    if((j=imat[vec1[i]][tm[i]])!=0)
      tot*=j;
    else
      return 0;
  }
  return tot;
}

/* same as unsterm1, except that all positive entries in the  */
/* submatrix are replaced with zeros */

int unsterm2(int **imat, int n, int m, int *vec1, int *tm, int k){
  int i,tot=1;
  for(i=0;i<k;i++)
    tot*=min(imat[vec1[i]][tm[i]], 0);
  return tot;
}



/* check if the k X k minor of n X m matrix imat indexed by the pair */
/* (vec1, vec2) is sign nonsingular or sign singular */
/* pms is a matrix storing all permutations of vec2 */
/* generated with allperms1. The final digit of each */
/* vector in pms is the parity of the permutation */

int minorisSNSSS(int **imat, int n, int m, int *vec1, int *vec2, int k, unsigned long fk, int **pms){
  int tmp, tm1=0;
  unsigned long i;

  /* worth checking for a row or column of zeros first */

  if(minorhas0rc(imat, n, m, vec1, vec2, k))
    return 1;

  for(i=0;i<fk;i++){
    if((tmp=unsterm1(imat, n, m, vec1, pms[i],k))){ /* nonzero term */
      if(!tm1){tm1=pms[i][k]*tmp;} /* store first nonzero term */
      else if(tm1*pms[i][k]*tmp<0) /* oppositely signed terms */
	return 0;
    }
  }
  return 1;
}

/* Exactly like the previous routine, but */
/* returns 1/-1 for SNS, 2 for SS and 0 for neither. */

int minorisSNSSS1(int **imat, int n, int m, int *vec1, int *vec2, int k, unsigned long fk, int **pms){
  int tmp, tm1=0;
  unsigned long i;

  /* worth checking for a row or column of zeros first */

  if(minorhas0rc(imat, n, m, vec1, vec2, k))
    return 2;

  for(i=0;i<fk;i++){
    if((tmp=unsterm1(imat, n, m, vec1, pms[i],k))){ /* nonzero term */
      if(!tm1){tm1=pms[i][k]*tmp;} /* store first nonzero term */
      else if(tm1*pms[i][k]*tmp<0) /* oppositely signed terms */
	return 0;
    }
  }
  if(!tm1)
    return 2;
  if(tm1<0)
    return -1;
  return 1;
}

/* like minorisSNSSS except doesn't first compute */
/* all the permutations, so should avoid memory failure */

int minorisSNSSS2(int **imat, int n, int m, int *vec1, int *vec2, int k){
  int tmp, tm1=0; 
  int *vec;
  int *veclr;
  int i, par=1, flag=1;
  unsigned long t=0;
  unsigned long fk=factorial(k);

  /* worth checking for a row or column of zeros first */

  if(minorhas0rc(imat, n, m, vec1, vec2, k))
    return 1;

  vec=(int *)malloc((size_t) ((n)*sizeof(int)));
  veclr=(int *)malloc((size_t) ((n)*sizeof(int)));

  for(i=0;i<k;i++){
    vec[i]=vec2[i];
    veclr[i]=-1;
  }

  while(flag){

    /* for(i=0;i<n;i++){ */
    /*   if(veclr[i]==-1) */
    /* 	fprintf(stderr, "<%d  ", vec[i]+1); */
    /*   else */
    /* 	fprintf(stderr, "%d>  ", vec[i]+1); */
    /* } */

    if((tmp=unsterm1(imat, n, m, vec1, vec,k))){ /* nonzero term */
      if(!tm1){tm1=par*tmp;} /* store first nonzero term */
      else if(tm1*par*tmp<0){ /* oppositely signed terms */
	free((char *) vec);free((char *) veclr);
	return 0;
      }
    }
    flag=nextperm(&vec, &veclr, &par, k);
    if(t%1000000==0){
      fprintf(stderr, "%d percent complete\n", (int)(((double)t)/((double)fk)*100));
    }
    t++;
  }
  free((char *) vec);free((char *) veclr);
  return 1;
}

/* Exactly like the previous routine, but */
/* returns 1/-1 for SNS, 2 for SS and 0 for neither. */

int minorisSNSSS3(int **imat, int n, int m, int *vec1, int *vec2, int k){
  int tmp, tm1=0; 
  int *vec;
  int *veclr;
  int i, par=1, flag=1;
  unsigned long t=0;
  unsigned long fk=factorial(k);

  /* worth checking for a row or column of zeros first */

  if(minorhas0rc(imat, n, m, vec1, vec2, k))
    return 2;

  vec=(int *)malloc((size_t) ((n)*sizeof(int)));
  veclr=(int *)malloc((size_t) ((n)*sizeof(int)));

  for(i=0;i<k;i++){
    vec[i]=vec2[i];
    veclr[i]=-1;
  }

  while(flag){

    /* for(i=0;i<k;i++){ */
    /*   if(veclr[i]==-1) */
    /* 	fprintf(stderr, "<%d  ", vec[i]+1); */
    /*   else */
    /* 	fprintf(stderr, "%d>  ", vec[i]+1); */
    /* } */
    /* fprintf(stderr, "\n"); */

    if((tmp=unsterm1(imat, n, m, vec1, vec,k))){ /* nonzero term */
      if(!tm1){tm1=par*tmp;} /* store first nonzero term */
      else if(tm1*par*tmp<0){ /* oppositely signed terms */
	free((char *) vec);free((char *) veclr);
	return 0;
      }
    }
    flag=nextperm(&vec, &veclr, &par, k);
    if(t%1000000==0){
      fprintf(stderr, "%d percent complete\n", (int)(((double)t)/((double)fk)*100));
    }
    t++;
  }
  free((char *) vec);free((char *) veclr);

  if(!tm1)
    return 2;
  if(tm1<0)
    return -1;
  return 1;

}

/* Check if the n times n matrix imat is sign nonsingular */
/* or sign singular */

int matrixisSNSSS(int **imat, int n){
  int vec1[n];
  unsigned long fn;
  int **pms;
  int i, retval=0;

  for(i=0;i<n;i++)
    vec1[i]=i;

  if(n>9){
    retval = minorisSNSSS3(imat, n, n, vec1, vec1, n);
  }
  else{
    fn=factorial(n);
    pms=allperms(n);

    retval = minorisSNSSS1(imat, n, n, vec1, vec1, n, fn, pms);
    free_imatrix(pms,0, fn-1, 0, n);
  }
  return retval;

}



/* Does the k X k minor of n X m matrix imat indexed by */
/* vec1 and vec2 have a row or column of zeros? */
/* Note that a failure does *not* imply that the minor is */
/* not sign singular - to really check this you'd need to check  */
/* if all terms in the determinant expansion are zero, a much  */
/* lengthier task */

int minorhas0rc(int **imat, int n, int m, int *vec1, int *vec2, int k){
  int flag, i=0, j=0;  

  while(j<k){
    flag=0;i=0;
    while(i<k && !flag){
      if(imat[vec1[i]][vec2[j]])
	flag=1; /* not a column of zeros */
      i++;
    }
    if(!flag)
      return 1; /* found a column of zeros */
    j++;
  }

  i=0;j=0;

  while(i<k){
    flag=0;j=0;
    while(j<k && !flag){
      if(imat[vec1[i]][vec2[j]])
	flag=1; /* not a row of zeros */
      j++;
    }
    if(!flag)
      return 1; /* found a row of zeros */
    i++;
  }

  return 0;
}

/* exactly as the previous routine, except checking imat-  */
/* (imat with all positive entries replaced with zeros) */

int minorhas0rc1(int **imat, int n, int m, int *vec1, int *vec2, int k){
  int flag, i=0, j=0;  

  while(j<k){
    flag=0;i=0;
    while(i<k && !flag){
      if(imat[vec1[i]][vec2[j]] <0)
	flag=1; /* not a column of zeros */
      i++;
    }
    if(!flag)
      return 1; /* found a column of zeros */
    j++;
  }

  i=0;j=0;

  while(i<k){
    flag=0;j=0;
    while(j<k && !flag){
      if(imat[vec1[i]][vec2[j]] <0)
	flag=1; /* not a row of zeros */
      j++;
    }
    if(!flag)
      return 1; /* found a row of zeros */
    i++;
  }

  return 0;
}


/* check if the k X k minor of n X m matrix imat indexed by  */
/* (vec1, vec2) is sign nonsingular or singular */
/* 1) calculate all permutations of vec2 */
/* 2) calculate the signed terms corresponding to vec1 and */
/*    each permutation of vec2 */
/* 3) Sum as we calculate to get the determinant */
/* 4) If determinant is zero then we are done */
/* 5) Otherwise multiply each term by the determinant to check if any */
/*    product is less than zero */

int minorisSNSsing(int **imat, int n, int m, int *vec1, int *vec2, int k, unsigned long fk, int **pms){
  int dt=0, tmp, tm1=0;
  unsigned long i;
  int flag=0;

  /* worth checking for a column or row of zeros first */

  if(minorhas0rc(imat, n, m, vec1, vec2, k))
    return 1;

  for(i=0;i<fk;i++){ // only store nonzero terms
    tmp=unsterm1(imat, n, m, vec1, pms[i],k);
    if (tmp){/* nonzero term */
      if(!tm1) /* first nonzero term */
	tm1=pms[i][k]*tmp;
      else if(pms[i][k]*tmp*tm1 <0)
	flag=1; /* not SNS */

      dt+=pms[i][k]*tmp;
    }
  }
  if(!flag || !dt) /* flag=0 means SNS; dt=0 means singular */
    return 1;

  return 0;
}

/* like the previous routine, except  */
/* does not require large memory storage */
/* for all permutations */

int minorisSNSsing1(int **imat, int n, int m, int *vec1, int *vec2, int k){
  int tmp, tm1=0; 
  int *vec;
  int *veclr;
  int dt=0, i, par=1, flag=0,nxt=1;
  unsigned long t=0;
  //  unsigned long fk=factorial(k);

  /* worth checking for a row or column of zeros first */

  if(minorhas0rc(imat, n, m, vec1, vec2, k)) /* sign singular */
    return 2;

  vec=(int *)malloc((size_t) ((n)*sizeof(int)));
  veclr=(int *)malloc((size_t) ((n)*sizeof(int)));

  for(i=0;i<k;i++){
    vec[i]=vec2[i];
    veclr[i]=-1;
  }

  while(nxt){

    /* for(i=0;i<k;i++){ */
    /*   if(veclr[i]==-1) */
    /* 	fprintf(stderr, "<%d  ", vec[i]+1); */
    /*   else */
    /* 	fprintf(stderr, "%d>  ", vec[i]+1); */
    /* } */
    /* fprintf(stderr, "\n"); */

    if ((tmp=unsterm1(imat, n, m, vec1, vec,k))){/* nonzero term */
      if(!tm1){tm1=par*tmp;} /* first nonzero term */
      else if(par*tmp*tm1 <0){
	flag=1; /* not SNS */
	//fprintf(stderr, "here\n");
      }

      dt+=par*tmp;
    }

    nxt=nextperm(&vec, &veclr, &par, k);
    if(t%1000000==0){
      //      fprintf(stderr, "%d percent complete\n", (int)(((double)t)/((double)fk)*100));
    }
    //fprintf(stderr, "t=%d\n", t);
    t++;
  }
  free((char *) vec);free((char *) veclr);

  if(!flag || !dt) /* flag=0 means SNS; dt=0 means singular */
    return 1;

  if(!tm1) /* SS */
    return 2;
  if(!flag){ /* SNS */
    if(tm1<0)
      return -1;
    else
      return 1;
  }
  if(!dt) /* singular but not SS */
    return 3;

  return 0; /* neither SNS nor singular */

}


int minorisSNSsing2(int **imat, int n, int m, int *vec1, int *vec2, int k){

  /* worth checking for a column or row of zeros first */

  if(minorhas0rc(imat, n, m, vec1, vec2, k))
    return 1;
  if(detsubmat(imat, n, m, vec1, vec2, k)==0)
    return 1;
  if(qualdetsubmat(imat, n, m, vec1, vec2, k)!=2)
    return 1;
  return 0;

}

// return the determinant

int minornotsing(int **imat, int n, int m, int *vec1, int *vec2, int k){

  /* worth checking for a column or row of zeros first */
  if(minorhas0rc(imat, n, m, vec1, vec2, k))
    return 0;
  return detsubmat(imat, n, m, vec1, vec2, k);

}

// return 1 for not SNS

int minornotSNS(int **imat, int n, int m, int *vec1, int *vec2, int k){

  /* worth checking for a column or row of zeros first */
  if(qualdetsubmat(imat, n, m, vec1, vec2, k)==2)
    return 1;
  return 0;

}



/* compatible minors */
/* take two matrices imat1, and imat2 of equal dimensions */
/* assume that imat1 is a matrix, while imat2 is a sign pattern */
/* 1) check if either of imat1[vec1|vec2] or imat2[vec1|vec2] */
/*    is SS. If so, then we are done. */
/* 2) check if imat1[vec1|vec2]=0. If so done. Else: */
/* 3) check if imat2(vec1|vec2) is SNS. If not, done. Else */
/* 4) check imat1[vec1|vec2]imat2[vec1|vec2] > 0 */

int submat_signpat_compat(int **imat1, int **imat2, int n, int m, int *vec1, int *vec2, int k){
  int **pms;
  int dt=0, tmp;
  unsigned long i, fk;
  int tm1=0;

  /* worth checking if either is SS first */

  if(minorhas0rc(imat2, n, m, vec1, vec2, k) || minorhas0rc(imat1, n, m, vec1, vec2, k))
    return 1;

  if((dt=detsubmat(imat1, n, m, vec1, vec2, k))==0)
    return 1;

  fk=factorial(k);
  pms=allperms1(vec2, k);
  for(i=0;i<fk;i++){
    tmp=unsterm1(imat2, n, m, vec1, pms[i],k);
    if (tmp!=0){
      if(abs(tmp)>1){ // unsigned term
	free_imatrix(pms,0, fk-1, 0, k);
	return 0;
      }
      // fails to be SNS
      if(tm1==0)
	tm1=pms[i][k]*tmp;
      else if(pms[i][k]*tmp*tm1 < 0){
	free_imatrix(pms,0, fk-1, 0, k);
	return 0;
      }
    }
  }
  if(dt*tm1<0){ // opp signs 
    free_imatrix(pms,0, fk-1, 0, k);
    return 0;
  }

  free_imatrix(pms,0, fk-1, 0, k);
  return 1;
}

/* imat1 is an n X m matrix, while imat2 is an n X m sign-pattern */
/* check if they are compatible. */

int mat_signpat_compat(int **imat1, int **imat2, int n, int m, int q){
  int k;
  int r;
  long r1, r2, cnk, cmk;
  int **xcombs;
  int **ycombs;
  int flag=1;


  for (r1=0;r1<n;r1++){
    for (r2=0;r2<m;r2++){
      if (imat1[r1][r2]*imat2[r1][r2]<0) 
	flag=0;
    }
  }
  if(!flag){
    fprintf(stderr, "\nThe matrices do not have compatible sign patterns.\n");
    return 0;
  }

  r=min(n,m);

  for(k=1;k<=r;k++){
    xcombs=allcombsgen(n,k);
    ycombs=allcombsgen(m,k);
    cnk=comb(n, k);
    cmk=comb(m, k);
    fprintf(stderr, "\nchecking %.0f minors of size %d...\n", ((double) cnk)*((double) (cmk)), k);
    for(r1=0;r1<cnk;r1++){
      if(r1%100==0)
	fprintf(stderr, ".");
      for(r2=0;r2<cmk;r2++){
	if(!submat_signpat_compat(imat1, imat2, n, m, xcombs[r1], ycombs[r2], k)){
	  fprintf(stderr, "\npair of submatrices which fail to be compatible:\n\n");
	  printsubmat(imat1, xcombs[r1], ycombs[r2], k, k);
	  fprintf(stderr, "*** and ***\n");
	  printsubmat(imat2, xcombs[r1], ycombs[r2], k, k);
	  flag=0; // exit here for speed
	  if(q){
	    free_imatrix(xcombs, 0, cnk-1, 0, k-1);
	    free_imatrix(ycombs, 0, cmk-1, 0, k-1);
	    return flag;
	  }
	}

      }
    }
    free_imatrix(xcombs, 0, cnk-1, 0, k-1);
    free_imatrix(ycombs, 0, cmk-1, 0, k-1);
  }


  return flag;
}

// are two matrices compatible - only check largest minors

int mat_signpat_compat1(int **imat1, int **imat2, int n, int m, int q){
  int k;
  int r;
  long r1, r2, cnk, cmk;
  int **xcombs;
  int **ycombs;
  int flag=1;


  r=min(n,m);

  k=r;
  
  xcombs=allcombsgen(n,k);
  fprintf(stderr, "here\n");
  ycombs=allcombsgen(m,k);
  fprintf(stderr, "here1\n");
  cnk=comb(n, k);
  cmk=comb(m, k);
  fprintf(stderr, "cnk= %.0f, cmk = %.0f\n", (double) cnk, (double)cmk);
  for(r1=0;r1<cnk;r1++){
    fprintf(stderr, "r1/cnk = %.0f/%.0f\n", (double)r1, (double)cnk);
    for(r2=0;r2<cmk;r2++){
      if(!submat_signpat_compat(imat1, imat2, n, m, xcombs[r1], ycombs[r2], k)){
	fprintf(stderr, "\npair of submatrices which fail to be compatible:\n\n");
	printsubmat(imat1, xcombs[r1], ycombs[r2], k, k);
	fprintf(stderr, "*** and ***\n");
	printsubmat(imat2, xcombs[r1], ycombs[r2], k, k);
	flag=0; // exit here for speed
	if(q){
	  free_imatrix(xcombs, 0, cnk-1, 0, k-1);
	  free_imatrix(ycombs, 0, cmk-1, 0, k-1);
	  return flag;
	}
      }

    }
  }
  free_imatrix(xcombs, 0, cnk-1, 0, k-1);
  free_imatrix(ycombs, 0, cmk-1, 0, k-1);


  return flag;
}


/* Check if n X m matrix imat1 and n X m sign patter imat2 */
/* are compatible. Remove empty rows/columns first. This is */
/* done twice because in theory each removal can create a new */
/* row/column of zeros in the other matrix */

int arecompat(int **imat1, int **imat2, int n, int m, int q){
  int **imat1a, **imat2a, **imat1b, **imat2b;
  int n1, m1, n2, m2;
  int flag=0;
  // simplify first
  simppair(imat1, imat2, n, m, &imat1a, &imat2a, &n1, &m1);
  simppair(imat1a, imat2a, n1, m1, &imat1b, &imat2b, &n2, &m2);
  printmat(imat1b, n2, m2);
  printmat(imat2b, n2, m2);
  fprintf(stderr, "Checking if matrices are compatible...\n");
  flag=mat_signpat_compat(imat1b, imat2b, n2, m2, q);
  if(flag)
    fprintf(stderr, "Finished checking compatibility: matrices are compatible.\n---------------------------------------\n");
  else
    fprintf(stderr, "Finished checking compatibility: matrices are not compatible.\n---------------------------------------\n");

  free_imatrix(imat1a, 0, n1-1, 0, m1-1);
  free_imatrix(imat2a, 0, n1-1, 0, m1-1);
  free_imatrix(imat1b, 0, n2-1, 0, m2-1);
  free_imatrix(imat2b, 0, n2-1, 0, m2-1);
  return flag;
}


/* Check if n X m matrix mat is SSD. The routine */
/* is a wrapper for isSSD1 adding a step removing */
/* columns/rows of zeros. */

int isSSD(int **mat, int n, int m, int q){
  int **imat;
  int n1, m1, flag=0;
  imat=simpmat(mat, n, m, &n1, &m1);
  fprintf(stderr, "Checking if the matrix is SSD...\n");
  if(n1+m > 15)
    fprintf(stderr, "...please be patient. This could take time.\n");
  flag=isSSD1(imat, n1, m1, q);
  if(flag)
    fprintf(stderr, "Finished checking SSD: the matrix is SSD.\n---------------------------------------\n");
  else
    fprintf(stderr, "Finished checking SSD: the matrix is not SSD.\n---------------------------------------\n");
  free_imatrix(imat,0, n1-1, 0, m1-1);
  return flag;
}

/* Check if n X m matrix mat is CSD (All square */
/* submatrices are sign nonsingular or sign singular) */


int isCSD1(int **imat, int n, int m, int q){
  int k,r;
  long r1, r2, cnk, cmk;
  int **xcombs;
  int **ycombs;
  int flag=1;
  unsigned long fk;
  int **pms;
  int ret;
  int rtrig=9;
  r=min(n,m);
  for(k=r;k>=2;k--){
    xcombs=allcombsgen(n,k);
    ycombs=allcombsgen(m,k);
    cnk=comb(n, k);
    cmk=comb(m, k);
    fk=factorial(k);
    for(r2=0;r2<cmk;r2++){
      if(k<=rtrig)
	pms=allperms1(ycombs[r2], k);
      for(r1=0;r1<cnk;r1++){
	if(k>rtrig) // large submatrix
	  ret=minorisSNSSS2(imat, n, m, xcombs[r1], ycombs[r2], k);
	else
	  ret=minorisSNSSS(imat, n, m, xcombs[r1], ycombs[r2], k, fk, pms);

	if(!ret){
	  fprintf(stderr, "submatrix which fails to be completely signed:\n");
	  printsubmat(imat, xcombs[r1], ycombs[r2], k, k);
	  flag=0; // exit here for speed
	  if(q){
	    if(k<=rtrig)
	      free_imatrix(pms,0, fk-1, 0, k);
	    free_imatrix(xcombs, 0, cnk-1, 0, k-1);
	    free_imatrix(ycombs, 0, cmk-1, 0, k-1);
	    return flag;
	  }
	}

      }
      if(k<=rtrig)
	free_imatrix(pms,0, fk-1, 0, k);
    }
    free_imatrix(xcombs, 0, cnk-1, 0, k-1);
    free_imatrix(ycombs, 0, cmk-1, 0, k-1);

  }


  return flag;
}

/* Check if n X m matrix imat is SSD */

int isSSD1(int **imat, int n, int m, int q){
  int k, ret=0;
  int r=100;
  long r1, r2, cnk, cmk;
  int **xcombs;
  int **ycombs;
  int flag=1;
  unsigned long fk;
  int **pms;

  r=min(n,m); 
  for(k=r;k>=2;k--){
 
    xcombs=allcombsgen(n,k);
    ycombs=allcombsgen(m,k);
    cnk=comb(n, k);cmk=comb(m, k);
    fprintf(stderr, "\nchecking %.0f minors of size %d...\n", ((double) cnk)*((double) (cmk)), k);
    fk=factorial(k);

    for(r2=0;r2<cmk;r2++){
      if(k<=10)
	pms=allperms1(ycombs[r2], k);
      for(r1=0;r1<cnk;r1++){
	if(k>10) // large submatrix
	  ret=minorisSNSsing1(imat, n, m, xcombs[r1], ycombs[r2], k);
	else
	  ret=minorisSNSsing(imat, n, m, xcombs[r1], ycombs[r2], k, fk, pms);

	if(!ret){
	  fprintf(stderr, "submatrix which fails to be sign-nonsingular or singular:\n");
	  printsubmat(imat, xcombs[r1], ycombs[r2], k, k);
	  flag=0; // exit here for speed
	  if(q){
	    if(k<=10)
	      free_imatrix(pms,0, fk-1, 0, k);
	    free_imatrix(xcombs, 0, cnk-1, 0, k-1);
	    free_imatrix(ycombs, 0, cmk-1, 0, k-1);
	    return flag;
	  }
	}

      }
      if(k<=10)
	free_imatrix(pms,0, fk-1, 0, k);
    }

    free_imatrix(xcombs, 0, cnk-1, 0, k-1);
    free_imatrix(ycombs, 0, cmk-1, 0, k-1);
  }
  return flag;
}

/* Check all minors of size minorsize of an n X m matrix imat */
/* to see if they are SNS or singular */
/* get rid of pms and fk - slow and unnecessary */

int checkminors(char *fname, int **imat, int **vmat, int n, int m, int minorsize){
  int k, ret=0, cutsz;
  long r1, r2, cnk, cmk;
  int **xcombs;
  int **ycombs;
  int flag=1, flag1;
  unsigned long fk;
  int **pms;
  int rtrig=10;
  FILE *fd, *fd1;
  long ct=1, ct1=1, ct2=1;
  int *xc, *yc;
  int failflag=0;
  long lab=1;
  int *cut1,*cut2;
  int totgraphs=0, totgold=0;
  int **allgraphs=NULL;

  fd = fopen(fname, "w");
  if(!fd){
    fprintf(stderr, "ERROR in checkminors: \"%s\" could not be opened for reading.\n", fname);
    return 0;
  }
  fd1=fopen("graphs.dot", "w");
  if(!fd1){
    fprintf(stderr, "ERROR in checkminors: \"graphs.dot\" could not be opened for reading.\n");
    return 0;
  }

  k=minorsize;
  fprintf(fd, "ratmx:true;\n");
  fprintf(stderr, "Dimensions: matrices are %d X %d\n\n", n,m);
  fprintf(stderr, "Expected rank = %d. Checking only minors of size %d\n\n", k, k);

  if(n>20 || m > 20){
    fprintf(stderr, "dimensions %d X %d are large. Routine may not terminate.\n ", n, m);

    xc=(int *)malloc((size_t) (k*sizeof(int)));
    yc=(int *)malloc((size_t) (k*sizeof(int)));

    firstcomb(xc, n, k);
    firstcomb(yc, m, k);
    flag=1;
    flag1=1;


    //    printvec(xc,k);
    //    printvec(yc,k);

    while(flag==1){
      flag1=1;
      while(flag1==1){
	ct1++;
	if(ct1%10000==0){
	  fprintf(stderr, "%ld e04\n", ct2);ct2++;
	  //	  printvec(xc,k);
	  //	  printvec(yc,k);
	}
	ret=minorisSNSsing2(imat, n, m, xc, yc, k);

	if(!ret){
	  //	  fprintf(stderr, "submatrix which fails to be sign-nonsingular or singular:\n");
	  //	  fprintf(stderr, "rows:    ");printvec(xc, k);fprintf(stderr, "columns: ");printvec(yc, k);
	  //	  printsubmat(imat, xc, yc, k, k);
	  //	  printmaximaindsubmat(fd, vmat, xc, yc, k, k);
	  //	  fprintf(stderr, "qualitative determinant = %d\n\n", qualdetsubmat(imat, n, m, xc, yc, k));

	  cutsz=cuttails(vmat, k, k, xc, yc, &cut1, &cut2);
	  totgold=totgraphs;
	  totgraphs=addtwointstr(totgraphs, cut1, cut2, k-cutsz, k-cutsz, &allgraphs);
	  if(totgraphs>totgold){
	    printdotindsubmat(fd1, vmat, cut1, cut2, k-cutsz, k-cutsz, lab);lab++;
	  }
	  free((char *) cut1);free((char *) cut2);

	//	  printdotindsubmat(fd1, vmat, xc, yc, k, k, lab);lab++;
	  failflag=1;
	  /* if(q){ */
	  /*   free((char *) xc);free((char *) yc); */
	  /*   return 0; */
	  /* } */
	}
      
	flag1=nextcomb(yc, m, k);
      }
      flag=nextcomb(xc, n, k);
    }
    free((char *) xc);free((char *) yc);
    
    if(!failflag)
      fprintf(stderr, "The second additive compound of the Jacobian is well behaved.\n");
    else
      fprintf(stderr, "The second additive compound of the Jacobian fails to be well behaved.\n");

    fclose(fd);fclose(fd1);
    return 0;
  }

  xcombs=allcombsgen(n,k);
  ycombs=allcombsgen(m,k);
  cnk=comb(n, k);cmk=comb(m, k);
  fprintf(stderr, "\nchecking %.0f minors of size %d...\n", ((double) cnk)*((double) (cmk)), k);
  fk=factorial(k);

  for(r2=0;r2<cmk;r2++){
    if(k<rtrig)
      pms=allperms1(ycombs[r2], k);
    for(r1=0;r1<cnk;r1++){
      ct1++;
      if(ct1%10000==0){
	fprintf(stderr, "%ld e04\n", ct2);ct2++;
	//	  printvec(xc,k);
	//	  printvec(yc,k);
      }

      if(k>=rtrig) // large submatrix
	ret=minorisSNSsing2(imat, n, m, xcombs[r1], ycombs[r2], k);
      else
	ret=minorisSNSsing(imat, n, m, xcombs[r1], ycombs[r2], k, fk, pms);
      
      if(!ret){
	fprintf(stderr, "submatrix which fails to be sign-nonsingular or singular:\n");
	for(ct=0;ct<k;ct++)
	  fprintf(stderr, "%d,",ycombs[r2][ct]+1);
	fprintf(stderr, "\n");

	printsubmat(imat, xcombs[r1], ycombs[r2], k, k);
	printmaximaindsubmat(fd, vmat, xcombs[r1], ycombs[r2], k, k);

	cutsz=cuttails(vmat, k, k, xcombs[r1], ycombs[r2], &cut1, &cut2);
	totgold=totgraphs;
	totgraphs=addtwointstr(totgraphs, cut1, cut2, k-cutsz, k-cutsz, &allgraphs);
	if(totgraphs>totgold){
	  printdotindsubmat(fd1, vmat, cut1, cut2, k-cutsz, k-cutsz, lab);lab++;
	}
	//	printdotindsubmat(fd1, vmat, xcombs[r1], ycombs[r2], n1, m1, lab);lab++;
       	free((char *) cut1);free((char *) cut2);
	flag=0; // exit here for speed
      }

    }
    if(k<rtrig)
      free_imatrix(pms,0, fk-1, 0, k);
  }

  for(k=0;k<totgraphs;k++)
    free ((char *)(allgraphs[k]));
  free((char *) allgraphs);

  free_imatrix(xcombs, 0, cnk-1, 0, k-1);
  free_imatrix(ycombs, 0, cmk-1, 0, k-1);
  fclose(fd);fclose(fd1);
  fprintf(stderr, "totgraphs = %d\n", totgraphs);
  return flag;
}

// return all permutations of [0...n-1]

int **allperms(int n){
  int **tmp;
  int **tmpmat;
  int i,j,k;
  unsigned long r;
  //  fprintf(stderr, "num = %d\n", factorial(n)-1);
  tmp=imatrix(0, factorial(n)-1, 0, n);

  //  fprintf(stderr, "num = %d\n", factorial(n)-1);
  if(n==1){
    tmp[0][0]=0;
    tmp[0][1]=1; // parity
    //    fprintf(stderr, "%d ", tmp[0][0]);
    //    fprintf(stderr, "\n");
  }
  else{
    tmpmat=allperms(n-1);
    i=0;
 
      
    for(k=0;k<n;k++){

      for(r=0;r<factorial(n-1);r++){
	i=n*r+k;
	for(j=0;j<n-k-1;j++){
	  tmp[i][j]=tmpmat[r][j];
	  //	  fprintf(stderr, "%d ", tmp[i][j]);
	}
	tmp[i][n-k-1]=n-1;
	//	fprintf(stderr, "%d ", tmp[i][n-k-1]);
	for(j=n-k;j<n;j++){
	  tmp[i][j]=tmpmat[r][j-1];
	  //	  fprintf(stderr, "%d ", tmp[i][j]);
	}
	tmp[i][j]=tmpmat[r][j-1]*par(k); // parity
	//	fprintf(stderr, "\n");	


      }
      //      fprintf(stderr, "\n");

    }
    free_imatrix(tmpmat, 0,factorial(n-1)-1,0,n-1);
  }
  return tmp;
}

/* recursive routine to return all permutations of an arbitrary  */
/* n-vector [vec1...vecn]. A final digit represents the parity  */
/* of the permutation: 1 for even, and -1 for odd */


int **allperms1(int *vec, int n){
  int **tmp;
  int **tmpmat;
  int i,j,k;
  unsigned long r;
  tmp=imatrix(0, factorial(n)-1, 0, n);
  if(n==1){
    tmp[0][0]=vec[0];
    tmp[0][1]=1; // parity
    //    fprintf(stderr, "%d ", tmp[0][0]);
    //    fprintf(stderr, "\n");
  }
  else{
    tmpmat=allperms1(vec, n-1);
    i=0;
    for(k=0;k<n;k++){
      for(r=0;r<factorial(n-1);r++){
	i=n*r+k;
	for(j=0;j<n-k-1;j++){
	  tmp[i][j]=tmpmat[r][j];
	  //	  fprintf(stderr, "%d ", tmp[i][j]);
	}
	tmp[i][n-k-1]=vec[n-1];
	//	fprintf(stderr, "%d ", tmp[i][n-k-1]);
	for(j=n-k;j<n;j++){
	  tmp[i][j]=tmpmat[r][j-1];
	  //	  fprintf(stderr, "%d ", tmp[i][j]);
	}
	tmp[i][j]=tmpmat[r][j-1]*par(k);  //parity
	//	fprintf(stderr, "\n");	
      }
      //     fprintf(stderr, "\n");
    }
    free_imatrix(tmpmat, 0,factorial(n-1)-1,0,n-1);
  }
  return tmp;
}

/* int **allcombsgen1(int n, int n1, int m){ */
/*   int **xcombs, **ycombs; */
/*   long tot1; */
/*   tot1=comb(n, n1); */
/*   xcombs=allcombsgen(n, n1); */
/*   ycombs=imatrix(0, tot1*pow(m, n1)-1, 0, n1-1); */
/*   for (i=0;i<tot1;i++){//each original combination */

/*   } */

/*   free_imatrix(xcombs, 0, comb(n, n1)-1, 0, n1-1); */
/*   return ycombs; */
/* } */


// generates all combinations of [0...n-1] of length n1, ordered
// recursive

int **allcombsgen(int n, int n1){
  int **tmp;
  int **tmpmat;
  int i,j,k;
  long cn, r;
  tmp=imatrix(0, comb(n, n1)-1, 0, n1-1);
  if(n1==1){
    for(k=0;k<n;k++){
      tmp[k][0]=k;
      //      fprintf(stderr, "%d ", tmp[k][0]);
      //      fprintf(stderr, "\n");
    }
  }
  
  else{
    tmpmat=allcombsgen(n, n1-1);

    i=0;
    cn=comb(n,n1-1);
    for(r=0;r<cn-1;r++){

      if(tmpmat[r][n1-2]<n-1){
	for(k=tmpmat[r][n1-2]+1;k<n;k++){
	  for(j=0;j<n1-1;j++){
	    tmp[i][j]=tmpmat[r][j];
	    //	    fprintf(stderr, "%d ", tmp[i][j]);
	  }
	  tmp[i][j]=k;
	  //	  fprintf(stderr, "%d\n", tmp[i][j]);
	  //	  fprintf(stderr, "position = %d\n", getpos(tmp[i], n, n1));
	  i++;
	}
      }
      //      fprintf(stderr, "\n");
    }
    free_imatrix(tmpmat, 0,cn-1,0,n1-2);
  }

  return tmp;
  
}

// get the position of a vector [i1, ...,ir] in a list
// of lexicographically ordered r-vectors chosen from {1,...,n}
// starts at position 0


long getpos(int *vec, int n, int r){
  long tot=0;
  int i,j, k1=n-1, k2=r-1, k3=0;
  for(j=0;j<r;j++){
    for(i=0;i<vec[j]-k3;i++){
      tot+=comb(k1, k2);
      k1--;
    }
    k1--;k2--;k3=vec[j]+1;
  }
  return tot;
}

// generates a submatrix of co-dimension 1 from matrix imat
// removing row i1 and column i2		      
 
int **submat1(int **imat, int n, int i1, int i2){
  int **tmp;
  int i,j,k,l;
  tmp=imatrix(0, n-2, 0, n-2);i=0;k=0;
  while(k<n){
    j=0;l=0;
    if(k!=i1){
      while(l<n){
	if(l!=i2){
	  tmp[i][j]=imat[k][l];
	  j++;
	}
	l++;
      }
      i++;
    }
    k++;
  }
  return tmp;
}

int pat(int r){
  if(r==0)
    return 0;
  return r/abs(r);
}

int S2(int **S, int n, int m, int ***imat1, int ***imat2, int *n1, int *m1){
  int **xcombs;
  int **tmp, **tmp1;
  long tot, r1, r2;
  int j;
  xcombs=allcombsgen(n,2);
  tot=comb(n,2);
  tmp=imatrix(0, tot-1, 0, 2*m*tot-1);
  tmp1=imatrix(0, tot-1, 0, 2*m*tot-1);

  for(r1=0;r1<tot;r1++){ //each row
    for(r2=0;r2<tot;r2++){ // each block
      //block one
      if(xcombs[r1][1]==xcombs[r2][1]){
	for(j=0;j<m;j++){
	  tmp[r1][2*m*r2+j]=S[xcombs[r1][0]][j];
	}
      }
      else if(xcombs[r1][0]==xcombs[r2][1]){
	for(j=0;j<m;j++){
	  tmp[r1][2*m*r2+j]=-S[xcombs[r1][1]][j];
	}
      }
      else{
	for(j=0;j<m;j++){
	  tmp[r1][2*m*r2+j]=0;
	}
      }

      // block 2
      if(xcombs[r1][0]==xcombs[r2][0]){
	for(j=0;j<m;j++){
	  tmp[r1][2*m*r2+m+j]=S[xcombs[r1][1]][j];
	}
      }
      else if(xcombs[r1][1]==xcombs[r2][0]){
	for(j=0;j<m;j++){
	  tmp[r1][2*m*r2+m+j]=-S[xcombs[r1][0]][j];
	}
      }
      else{
	for(j=0;j<m;j++){
	  tmp[r1][2*m*r2+m+j]=0;
	}
      }


      //imat2

      if(xcombs[r1][0]==xcombs[r2][0] && xcombs[r1][1]==xcombs[r2][1]){
	for(j=0;j<m;j++){
	  tmp1[r1][2*m*r2+j]=pat(S[xcombs[r1][0]][j]);
	  tmp1[r1][2*m*r2+m+j]=pat(S[xcombs[r1][1]][j]);
	}
      }

    }
  }

  free_imatrix(xcombs, 0, comb(n, 2)-1, 0, 1);


  simppair(tmp, tmp1, comb(n,2), 2*m*comb(n,2), imat1, imat2, n1, m1); 
    free_imatrix(tmp, 0, comb(n,2)-1, 0, 2*m*comb(n,2)-1);
    free_imatrix(tmp1, 0, comb(n,2)-1, 0, 2*m*comb(n,2)-1);


  return 1;
}

// like S2, but with V explicitly

int S2a(int **S, int **V, int n, int m, int ***imat1, int ***imat2, int *n1, int *m1){
  int **xcombs;
  int **tmp, **tmp1;
  long tot, r1, r2;
  int j;
  xcombs=allcombsgen(n,2);
  tot=comb(n,2);
  tmp=imatrix(0, tot-1, 0, 2*m*tot-1);
  tmp1=imatrix(0, tot-1, 0, 2*m*tot-1);

  for(r1=0;r1<tot;r1++){ //each row
    for(r2=0;r2<tot;r2++){ // each block
      //block one
      if(xcombs[r1][1]==xcombs[r2][1]){
	for(j=0;j<m;j++){
	  tmp[r1][2*m*r2+j]=S[xcombs[r1][0]][j];
	}
      }
      else if(xcombs[r1][0]==xcombs[r2][1]){
	for(j=0;j<m;j++){
	  tmp[r1][2*m*r2+j]=-S[xcombs[r1][1]][j];
	}
      }
      else{
	for(j=0;j<m;j++){
	  tmp[r1][2*m*r2+j]=0;
	}
      }

      // block 2
      if(xcombs[r1][0]==xcombs[r2][0]){
	for(j=0;j<m;j++){
	  tmp[r1][2*m*r2+m+j]=S[xcombs[r1][1]][j];
	}
      }
      else if(xcombs[r1][1]==xcombs[r2][0]){
	for(j=0;j<m;j++){
	  tmp[r1][2*m*r2+m+j]=-S[xcombs[r1][0]][j];
	}
      }
      else{
	for(j=0;j<m;j++){
	  tmp[r1][2*m*r2+m+j]=0;
	}
      }


      //imat2

      if(xcombs[r1][0]==xcombs[r2][0] && xcombs[r1][1]==xcombs[r2][1]){
	for(j=0;j<m;j++){
	  tmp1[r1][2*m*r2+j]=pat(V[xcombs[r1][0]][j]);
	  tmp1[r1][2*m*r2+m+j]=pat(V[xcombs[r1][1]][j]);
	}
      }

    }
  }

  free_imatrix(xcombs, 0, comb(n, 2)-1, 0, 1);


  simppair(tmp, tmp1, comb(n,2), 2*m*comb(n,2), imat1, imat2, n1, m1); 
  free_imatrix(tmp, 0, comb(n,2)-1, 0, 2*m*comb(n,2)-1);
  free_imatrix(tmp1, 0, comb(n,2)-1, 0, 2*m*comb(n,2)-1);


  return 1;
}

int S2b(int **S, int **V, int n, int m, int ***imat1, int **valsvec, int **indsvec, int *base, int *basetot, int *n1, int *m1){
  int **xcombs;
  long tot, r1, r2;
  int j, totlen;
  int *tmpvals, *tmpinds;
  int **tmpmat;
  xcombs=allcombsgen(n,2);
  tot=comb(n,2);

  tmpmat=imatrix(0, tot-1, 0, 2*m*tot-1);
  tmpvals=(int *)malloc(sizeof(int*) *(2*m*tot));
  tmpinds=(int *)malloc(sizeof(int*) *(2*m*tot));

  for(r1=0;r1<tot;r1++){ //each row
    for(r2=0;r2<tot;r2++){ // each block
      //block one
      if(xcombs[r1][1]==xcombs[r2][1]){
	for(j=0;j<m;j++){
	  tmpmat[r1][2*m*r2+j]=S[xcombs[r1][0]][j];
	}
      }
      else if(xcombs[r1][0]==xcombs[r2][1]){
	for(j=0;j<m;j++){
	  tmpmat[r1][2*m*r2+j]=-S[xcombs[r1][1]][j];
	}
      }
      else{
	for(j=0;j<m;j++){
	  tmpmat[r1][2*m*r2+j]=0;
	}
      }

      // block 2
      if(xcombs[r1][0]==xcombs[r2][0]){
	for(j=0;j<m;j++){
	  tmpmat[r1][2*m*r2+m+j]=S[xcombs[r1][1]][j];
	}
      }
      else if(xcombs[r1][1]==xcombs[r2][0]){
	for(j=0;j<m;j++){
	  tmpmat[r1][2*m*r2+m+j]=-S[xcombs[r1][0]][j];
	}
      }
      else{
	for(j=0;j<m;j++){
	  tmpmat[r1][2*m*r2+m+j]=0;
	}
      }


      //valsvec and indsvec // change - only do once

      if(xcombs[r1][0]==xcombs[r2][0] && xcombs[r1][1]==xcombs[r2][1]){
	for(j=0;j<m;j++){
	  tmpvals[2*m*r2+j]=pat(V[xcombs[r1][0]][j]);
	  tmpvals[2*m*r2+m+j]=pat(V[xcombs[r1][1]][j]);

	  tmpinds[2*m*r2+j]=xcombs[r1][0]*m+j;
	  tmpinds[2*m*r2+m+j]=xcombs[r1][1]*m+j;
	}
      }

    }
  }

  totlen=0;
  // get rid of zeros
  for(r2=0;r2<tot;r2++){
    for(j=0;j<2*m;j++){
      if(tmpvals[2*m*r2+j]!=0)
	totlen++;
    }
  }



  (*valsvec)=(int *)malloc(sizeof(int*) *(totlen));
  (*indsvec)=(int *)malloc(sizeof(int*) *(totlen));
  (*imat1)=imatrix(0, tot-1, 0, totlen-1);
  (*n1)=tot;(*m1)=totlen;

  totlen=0;
  for(r2=0;r2<tot;r2++){
    base[r2]=0;
    for(j=0;j<2*m;j++){
      if(tmpvals[2*m*r2+j]!=0){
	(*valsvec)[totlen]=tmpvals[2*m*r2+j];
	(*indsvec)[totlen++]=tmpinds[2*m*r2+j];
	(base[r2])++;
      }
    }
  }
  
  basetot[0]=0; // running total
  for(r2=1;r2<tot;r2++){
    basetot[r2]=basetot[r2-1]+base[r2-1];
  }

  for(r1=0;r1<tot;r1++){ //each row
    totlen=0;
    for(r2=0;r2<tot;r2++){ // each block
      for(j=0;j<2*m;j++){
	if(tmpvals[2*m*r2+j]!=0){
	  (*imat1)[r1][totlen++]=tmpmat[r1][2*m*r2+j];
	}

      }
    }
  }

  printvec(tmpvals, 2*m*tot);
  printvec((*valsvec), (*m1));
  printvec(base, tot);
  printvec(basetot, tot);


  free((char *) tmpvals);
  free((char *) tmpinds);
  free_imatrix(tmpmat, 0, tot-1, 0, 2*m*tot-1);
  free_imatrix(xcombs, 0, comb(n, 2)-1, 0, 1);

  return 1;
}




//detsk has dimension (0, comb(n, k-1), 0, comb(m, k-1))

int **detsk(int **imat, int n, int m, int k){
  int **xcombs;
  int **ycombs;
  int **tmp;
  long r1, r2;
  int **submat;

  tmp=imatrix(0, comb(n, k-1), 0, comb(m, k-1));
  xcombs=allcombsgen(n,k);
  ycombs=allcombsgen(m,k);
  for(r1=0;r1<comb(n, k);r1++){
    for(r2=0;r2<comb(m, k);r2++){
      submat = submatgen(imat, n, m, xcombs[r1], ycombs[r2],k);
      tmp[r1][r2]=det(submat,k);
      free_imatrix(submat, 0, k-1, 0, k-1);
    }
  }

  free_imatrix(xcombs, 0, comb(n, k)-1, 0, k-1);
  free_imatrix(ycombs, 0, comb(m, k)-1, 0, k-1);
  return tmp;

}

/* Check if n X m matrix mat is CSD. The routine */
/* is a wrapper for isCSD1 adding a step removing */
/* columns/rows of zeros. */

int isCSD(int **mat, int n, int m, int q){
  int **imat;
  int n1, m1, flag;
  imat=simpmat(mat, n, m, &n1, &m1);
  fprintf(stderr, "Checking if the matrix is CSD...\n");
  printmat(imat, n1, m1);
  flag=isCSD1(imat, n1, m1, q);
  if(flag==1)
    fprintf(stderr, "Finished checking CSD: the matrix is CSD.\n---------------------------------------\n");
  else
    fprintf(stderr, "Finished checking CSD: the matrix is not CSD.\n---------------------------------------\n");
  free_imatrix(imat,0, n1-1, 0, m1-1);
  return flag;
}

/* Convert n X m matrix imat into irreversible form: */
/* namely replace each column C with */
/* the column pair C|-C */


int **doublemat(int **imat, int n, int m){
  int **tmp;
  int i,j;
  tmp=imatrix(0, n-1, 0, 2*m-1);
  for(i=0;i<n;i++){
    for(j=0;j<m;j++){
      tmp[i][2*j]=imat[i][j];
      tmp[i][2*j+1]=-imat[i][j];
    }
  }
  return tmp;
}

int doubleisWSD2(int **imat, int n, int m){
  int **tmp;int flag=0;
  tmp=doublemat(imat, n, m);
  printmat(tmp, n, 2*m);
  if (isWSD2(tmp, n, 2*m))
    flag=1;
  free_imatrix(tmp, 0, n-1, 0, 2*m-1);
  return flag;

}

// WSD1 is faster than WSD2

int doubleisWSD(int **mat, int n, int m, int q){
  int **imat, **tmp;
  int n1,m1,flag=0;
  imat=simpmat(mat, n, m, &n1, &m1); // first simplify and then double
  tmp=doublemat(imat, n1, m1);
  fprintf(stderr, "Checking if this matrix is WSD...\n");
  printmat(tmp, n1, 2*m1);
  flag=isWSD1(tmp, n1, 2*m1, q);
  if(flag==1)
    fprintf(stderr, "Finished checking WSD: the matrix is WSD.\n---------------------------------------\n");
  else
    fprintf(stderr, "Finished checking WSD: the matrix is not WSD.\n---------------------------------------\n");
  free_imatrix(tmp, 0, n-1, 0, 2*m-1);
  free_imatrix(imat,0, n1-1, 0, m1-1);
  return flag;

}

/* A wrapper for isWSD1, first removing rows/columns of zeros */

int isWSD(int **mat, int n, int m, int q){
  int **imat;
  int n1,m1,flag=0;
  imat=simpmat(mat, n, m, &n1, &m1);
  flag=isWSD1(imat, n1, 2*m1, q);
  free_imatrix(imat,0, n1-1, 0, m1-1);
  return flag;

}


int isWSD2(int **imat, int n, int m){
  int k,r;
  long r1,r2,cnk,cmk;
  int **xcombs;
  int **ycombs;
  int **submat;
  int flag=1;
  r=min(n,m);
  for(k=2;k<=r;k++){
    xcombs=allcombsgen(n,k);
    ycombs=allcombsgen(m,k);
    cnk=comb(n,k);cmk=comb(m,k);
    for(r1=0;r1<cnk;r1++){
      for(r2=0;r2<cmk;r2++){
	submat = submatgen(imat, n, m, xcombs[r1], ycombs[r2],k);
	if(!isWS(submat, k)){
	  fprintf(stderr, "submatrix which fails to det(S)det(S-)>=0:\n");
	  printmat(submat, k,k);
	  flag=0; // exit here for speed
	}
	free_imatrix(submat, 0, k-1, 0, k-1);
      }
    }

    free_imatrix(xcombs, 0, cnk-1, 0, k-1);
    free_imatrix(ycombs, 0, cmk-1, 0, k-1);
  }
  return flag;
}

/* Direct computation on n X m matrix imat to see */
/* if it is WSD. The first step (WSpair) is to find the */
/* matrix of exponents, and to remove rows/columns of zeros. */
/* Then compatibility is checked. */

int isWSD1(int **imat, int n, int m, int q){
  int k;
  int r=100;
  long r1,r2,cnk,cmk;
  int **xcombs;
  int **ycombs;
  int flag=1;
  int **imat2, **imat3, n1, m1;

  WSpair(imat, n, m, &imat2, &imat3, &n1, &m1); // imat3 is sparser
  fprintf(stderr, "The simplified matrix pair...\n");
  printmat(imat2, n1, m1);
  printmat(imat3, n1, m1);

  r=min(n1,m1);

  for(k=2;k<=r;k++){
    xcombs=allcombsgen(n1,k);
    ycombs=allcombsgen(m1,k);
    cnk=comb(n1,k);cmk=comb(m1,k);
    for(r1=0;r1<cnk;r1++){
      for(r2=0;r2<cmk;r2++){
	if(fixedminorcompat(imat2, imat3, n1, m1, xcombs[r1], ycombs[r2],k)<0){
	  fprintf(stderr, "submatrix which fails det(S)det(S-) >=0:\n");
	  printsubmat(imat, xcombs[r1], ycombs[r2], k, k);
	  flag=0; 
	  if(q){
	    free_imatrix(xcombs, 0, cnk-1, 0, k-1);
	    free_imatrix(ycombs, 0, cmk-1, 0, k-1);
	    if(imat2)
	      free_imatrix(imat2, 0, n1-1, 0, m1-1);
	    if(imat3)
	      free_imatrix(imat3, 0, n1-1, 0, m1-1);
	    return flag;
	  }
	}
      }
    }

    free_imatrix(xcombs, 0, cnk-1, 0, k-1);
    free_imatrix(ycombs, 0, cmk-1, 0, k-1);
  }

  if(imat2)
    free_imatrix(imat2, 0, n1-1, 0, m1-1);
  if(imat3)
    free_imatrix(imat3, 0, n1-1, 0, m1-1);


  return flag;
}



/* Check if n X m matrices imat1 and imat2 are compatible */

int mats_compat(int **imat1, int **imat2, int n, int m, int q){
  int k;
  int r=100,s1,s2;
  long r1,r2,cnk,cmk;
  int **xcombs;
  int **ycombs;
  int flag=1;
  int **imat1a, **imat2a;
  int n1, m1;

  fprintf(stderr, "Checking if mass-action matrices are compatible.\n\n");

  simppair(imat1, imat2, n, m, &imat1a, &imat2a, &n1, &m1); // imat2a is sparser
  printmat(imat1a, n1, m1);
  printmat(imat2a, n1, m1);

  r=min(n1,m1);

  for(k=2;k<=r;k++){
    xcombs=allcombsgen(n1,k);
    ycombs=allcombsgen(m1,k);
    cnk=comb(n1,k);cmk=comb(m1,k);
    for(r1=0;r1<cnk;r1++){
      for(r2=0;r2<cmk;r2++){
	//	fprintf(stderr, "in here\n");
	if(fixedminorcompat(imat1a, imat2a, n1, m1, xcombs[r1], ycombs[r2],k) < 0){
	  //	  fprintf(stderr, "\npair of submatrices with determinants of opposite sign:\n\n");
	  for(s1=0;s1<k;s1++){
	    for(s2=0;s2<k;s2++){
	      fprintf(stderr, "%2d  ",imat1a[xcombs[r1][s1]][ycombs[r2][s2]]);
	    }
	    fprintf(stderr, "\n");
	  }
	  fprintf(stderr, "*** and ***\n");
	  for(s1=0;s1<k;s1++){
	    for(s2=0;s2<k;s2++){
	      fprintf(stderr, "%2d  ",imat2a[xcombs[r1][s1]][ycombs[r2][s2]]);
	    }
	    fprintf(stderr, "\n");
	  }

	  for(s1=0;s1<k;s1++){
	    fprintf(stderr, "[ ");
	    for(s2=0;s2<k;s2++)
	      fprintf(stderr, "%d,%d  ",xcombs[r1][s1],ycombs[r2][s2]);
	    fprintf(stderr, "]\n");
	  }

	  fprintf(stderr, "_____________________\n\n");
	  flag=0; // could return here for speed
	  if(q){
	    free_imatrix(xcombs, 0, cnk-1, 0, k-1);
	    free_imatrix(ycombs, 0, cmk-1, 0, k-1);
	    if(imat1a)
	      free_imatrix(imat1a, 0, n1-1, 0, m1-1);
	    if(imat2a)
	      free_imatrix(imat2a, 0, n1-1, 0, m1-1);

	    fprintf(stderr, "Finished checking if mass-action matrices are compatible: matrices are not compatible.\n---------------------------------------\n");
	    return flag;
	  }
	}
      }
    }

    free_imatrix(xcombs, 0, cnk-1, 0, k-1);
    free_imatrix(ycombs, 0, cmk-1, 0, k-1);
  }

  if(imat1a)
    free_imatrix(imat1a, 0, n1-1, 0, m1-1);
  if(imat2a)
    free_imatrix(imat2a, 0, n1-1, 0, m1-1);

  if(flag==1)
    fprintf(stderr, "Finished checking if mass-action matrices are compatible: matrices are compatible.\n---------------------------------------\n");
  else
    fprintf(stderr, "Finished checking if mass-action matrices are compatible: matrices are not compatible.\n---------------------------------------\n");

  return flag;
}

// the Mass-action Jacobian is sign nonsingular


int allminorsigns(int **imat1, int **imat2, int n, int m, int q){
  int k, r, flag;
  long r1,r2,cnk,cmk;
  int **xcombs;
  int **ycombs;
  int val;
  int sgn=0;
  long countp, countn;

  fprintf(stderr, "Checking all minors for mass-action Jacobian.\n\n");

  r=min(n,m);

  for(k=1;k<=r;k++){
    sgn=0;
    xcombs=allcombsgen(n,k);
    ycombs=allcombsgen(m,k);
    cnk=comb(n,k);cmk=comb(m,k);
    flag=0;countp=0;countn=0;
    //    for(r1=0;(r1<cnk  && flag==0);r1++){
    //      for(r2=0;(r2<cmk && flag==0);r2++){
    for(r1=0;r1<cnk;r1++){
      for(r2=0;r2<cmk;r2++){
	val=fixedminorcompat(imat1, imat2, n, m, xcombs[r1], ycombs[r2],k);
	if(val>0)
	  countp++;
	else if(val<0)
	  countn++;
	if(sgn==0 && val > 0){
	  sgn=1;
	}
	else if(sgn==0 && val < 0){
	  sgn=-1;
	}
	else if((sgn >0 && val < 0) || (sgn <0 && val > 0)){
	  flag=1;
	}
   
      }
    }
    free_imatrix(xcombs, 0, cnk-1, 0, k-1);
    free_imatrix(ycombs, 0, cmk-1, 0, k-1);

    if(flag==1){
      fprintf(stderr, "Minors of dimension %d: %.2f %% (%.0f:%.0f) of all nonzero minors are positive.\n", k, 100.0*((double) countp)/((double)countp + (double)countn), (double) countp, (double) countn);
    }
    else if(sgn==0){
      fprintf(stderr, "Minors of dimension %d are all zero.\n", k);
    }
    else if(sgn==1){
      fprintf(stderr, "Nonzero minors of dimension %d are all positive (%.0f).\n", k, (double) countp);
    }
    else if(sgn==-1){
      fprintf(stderr, "Nonzero minors of dimension %d are all negative (%.0f).\n", k, (double) countn);
    }

  }


  return 1;
}

//
// convert an integer matrix to a symbolic matrix
//

matrix imattoexmat(int **A, int n, int m){
  matrix C(n,m);
  int i,j;
  for(i=0;i<n;i++){
    for(j=0;j<m;j++)
      C(i,j)=A[i][j];
  }
  return C;
}

int MAisSNS(int **imat1, int **imat2, int n, int m, int q){
  int k;
  long r1,r2,cnk,cmk;
  int **xcombs;
  int **ycombs;
  int val;
  int sgn=0;

  fprintf(stderr, "Checking if mass-action Jacobian is SNS.\n\n");

  if(m<n)
    return 0;

  k=n;
  xcombs=allcombsgen(n,k);
  ycombs=allcombsgen(m,k);
  cnk=comb(n,k);cmk=comb(m,k);
  for(r1=0;r1<cnk;r1++){
    for(r2=0;r2<cmk;r2++){
      //      fprintf(stderr, "%.2f percent\n", ((double)r2)/((double) cmk));
      val=fixedminorcompat(imat1, imat2, n, m, xcombs[r1], ycombs[r2],k);
      //      fprintf(stderr, "here\n");
      if(sgn==0 && val > 0)
	sgn=1;
      else if(sgn==0 && val < 0)
	sgn=-1;
      else if((sgn >0 && val < 0) || (sgn <0 && val > 0)){
	free_imatrix(xcombs, 0, cnk-1, 0, k-1);
	free_imatrix(ycombs, 0, cmk-1, 0, k-1);
	fprintf(stderr, "Finished checking if mass-action Jacobian is SNS: it has terms of opposite sign.\n---------------------------------------\n");
	return 0;
      }
   
    }
  }

  free_imatrix(xcombs, 0, cnk-1, 0, k-1);
  free_imatrix(ycombs, 0, cmk-1, 0, k-1);


  if(sgn==0){
    fprintf(stderr, "Finished checking if mass-action Jacobian is SNS: it is sign singular.\n---------------------------------------\n");
    return 0;
  }
  else
    fprintf(stderr, "Finished checking if mass-action Jacobian is SNS: it is.\n---------------------------------------\n");

  return 1;
}






int isWS(int **mat, int n){
  int **tmp;
  tmp=submatminus(mat, n, n);


/*   printmat(mat, n, n); */
/*   fprintf(stderr, "\n"); */
/*   printmat(tmp, n, n); */
/*   fprintf(stderr, "det1 = %d,det2 = %d\n", det(mat, n),det(tmp, n)); */
/*   fprintf(stderr, "\n  ****  \n"); */


  if(det(mat, n)*det(tmp, n) >=0){
    free_imatrix(tmp, 0, n-1, 0, n-1);
    return 1;
  }
  free_imatrix(tmp, 0, n-1, 0, n-1);
  return 0;
}

// get rid of empty rows and columns

void simppair(int **imat, int **imat1, int n, int m, int ***imat2, int ***imat3, int *n1, int *m1){
  int i, j, s1, s2, flag;
  int *emptycols;
  int *emptyrows;
  (*n1)=0;(*m1)=0;



  emptycols=(int *)malloc((size_t) ((m)*sizeof(int)));
  emptyrows=(int *)malloc((size_t) ((n)*sizeof(int)));

  for(i=0;i<n;i++)
    emptyrows[i]=1;

  for(i=0;i<m;i++)
    emptycols[i]=1;


  i=0;j=0;

  while(j<m){
    i=0;flag=0;
    while(i<n && flag==0){
      if(imat[i][j] !=0){
	flag=1; // not a column of zeros
      }
      i++;
    }
    if(flag==1){
      i=0;
      while(i<n && emptycols[j] == 1){
	if(imat1[i][j] !=0){
	  emptycols[j]=0; // not a column of zeros
	  (*m1)++;
	}
	i++;
      }
    }

    j++;
  }
  fprintf(stderr, "numcols = %d --> %d\n", m, (*m1));

  i=0;j=0;

  while(i<n){
    j=0;flag=0;
    while(j<m && flag==0){
      if(imat[i][j] !=0)
	flag=1; // not a row of zeros
      j++;
    }
    if(flag==1){
      j=0;
      while(j<m && emptyrows[i] == 1){
	if(imat1[i][j] !=0){
	  emptyrows[i]=0; // not a row of zeros
	  (*n1)++;
	}
	j++;
      }
    }
    i++;
  }
  fprintf(stderr, "numrows = %d --> %d\n", n, (*n1));

  if((*n1)==0 || (*m1)==0){ // one matrix is a matrix of zeros
    (*imat2)=NULL;(*imat3)=NULL;
    free((FREE_ARG) (emptyrows));free((FREE_ARG) (emptycols));
    return;
  }


  (*imat2)=imatrix(0, (*n1)-1, 0, (*m1)-1);
  (*imat3)=imatrix(0, (*n1)-1, 0, (*m1)-1);

  s1=0;
  for(i=0;i<n;i++){
    s2=0;
    if(emptyrows[i]==0){ // not empty row
      for(j=0;j<m;j++){
	if(emptycols[j]==0){ // not empty column
	  (*imat2)[s1][s2]=imat[i][j];
	  (*imat3)[s1][s2]=imat1[i][j];
	  s2++;
	}
      }
      s1++;
    }
  }


  free((FREE_ARG) (emptyrows));
  free((FREE_ARG) (emptycols));

  return;
}

// get rid of empty rows and columns
// from an integer matrix and a matrix of expressions

void simppair1(int **imat, matrix imat1, int n, int m, int ***imat2, matrix *imat3, int *n1, int *m1){
  int i, j, s1, s2, flag;
  int *emptycols;
  int *emptyrows;
  (*n1)=0;(*m1)=0;



  emptycols=(int *)malloc((size_t) ((m)*sizeof(int)));
  emptyrows=(int *)malloc((size_t) ((n)*sizeof(int)));

  for(i=0;i<n;i++)
    emptyrows[i]=1;

  for(i=0;i<m;i++)
    emptycols[i]=1;


  i=0;j=0;

  while(j<m){
    i=0;flag=0;
    while(i<n && flag==0){
      if(imat[i][j] !=0){
	flag=1; // not a column of zeros
      }
      i++;
    }
    if(flag==1){
      i=0;
      while(i<n && emptycols[j] == 1){
	if(imat1(i,j)!=0){
	  emptycols[j]=0; // not a column of zeros
	  (*m1)++;
	}
	i++;
      }
    }
    j++;
  }
  fprintf(stderr, "numcols = %d --> %d\n", m, (*m1));

  i=0;j=0;

  while(i<n){
    j=0;flag=0;
    while(j<m && flag==0){
      if(imat[i][j] !=0)
	flag=1; // not a row of zeros
      j++;
    }

    if(flag==1){
      j=0;
      while(j<m && emptyrows[i] == 1){
	if(imat1(i,j)!=0){
	  emptyrows[i]=0; // not a row of zeros
	  (*n1)++;
	}
	j++;
      }
    }
    i++;
  }
  fprintf(stderr, "numrows = %d --> %d\n", n, (*n1));

  if((*n1)==0 || (*m1)==0){ // one matrix is a matrix of zeros
    (*imat2)=NULL;(*imat3)=NULL;
    free((FREE_ARG) (emptyrows));free((FREE_ARG) (emptycols));
    return;
  }


  (*imat2)=imatrix(0, (*n1)-1, 0, (*m1)-1);
  (*imat3)=matrix(*n1, *m1);

  s1=0;
  for(i=0;i<n;i++){
    s2=0;
    if(emptyrows[i]==0){ // not empty row
      for(j=0;j<m;j++){
	if(emptycols[j]==0){ // not empty column
	  (*imat2)[s1][s2]=imat[i][j];
	  (*imat3)(s1,s2)=imat1(i,j);
	  s2++;
	}
      }
      s1++;
    }
  }
  free((FREE_ARG) (emptyrows));
  free((FREE_ARG) (emptycols));

  return;
}

void WSpair(int **imat, int n, int m, int ***imat2, int ***imat3, int *n1, int *m1){
  int i, j;
  int **tmpmat;
  (*n1)=0;(*m1)=0;

  tmpmat=imatrix(0, n-1, 0, m-1);

  for(i=0;i<n;i++)
    for(j=0;j<m;j++)
      tmpmat[i][j]=min(imat[i][j], 0);

  simppair(tmpmat, imat, n, m, imat3, imat2, n1, m1); // tmpmat --> imat3 (sparser matrix)

  free_imatrix(tmpmat, 0, n-1, 0, m-1);
  return;

}



int minorisWS(int **imat1, int **imat2, int n, int m, int *vec1, int *vec2, int k, unsigned long fk){
  int **pms;
  int dt=0, dt2=0;
  unsigned long i;

  // imat2 is sparser
  if(minorhas0rc(imat2, n, m, vec1, vec2, k))
    return 1;

  pms=allperms1(vec2, k);

  for(i=0;i<fk;i++)
    dt2+=unsterm1(imat2, n, m, vec1, pms[i],k)*pms[i][k];

  if(dt2!=0){
    for(i=0;i<fk;i++)
      dt+=unsterm1(imat1, n, m, vec1, pms[i],k)*pms[i][k];

    if(dt*dt2<0){
      //      fprintf(stderr, "dt = %d, dt2=%d\n", dt, dt2);
      free_imatrix(pms,0, fk-1, 0, k);
      return 0;

    }
  }
  free_imatrix(pms,0, fk-1, 0, k);
  return 1;
}

/* This function returns the sign of the product of the  */
/* two k X k submatrices imat1(vec1|vec2) and imat2(vec1|vec2)  */
/* of n X m matrices imat1 and imat2 */

int fixedminorcompat(int **imat1, int **imat2, int n, int m, int *vec1, int *vec2, int k){
  int dt=0, dt2=0, prod;

  if(minorhas0rc(imat2, n, m, vec1, vec2, k) || minorhas0rc(imat1, n, m, vec1, vec2, k))
    return 0;

  dt2=detsubmat(imat2, n, m, vec1, vec2, k);
  if(dt2!=0){
    dt=detsubmat(imat1, n, m, vec1, vec2, k);
    if((prod=dt*dt2)<0)
      return -1;
    else if(prod>0)
      return 1;
    else
      return 0;
  }
  return 0;
}

int minorisnonzero(int **imat1, int n, int m, int *vec1, int *vec2, int k){
  int dt=0;
  if(minorhas0rc(imat1, n, m, vec1, vec2, k) || ((dt=detsubmat(imat1, n, m, vec1, vec2, k))==0))
    return 0;
  if(dt>0)
    return 1;
  return -1;

}
/* copy an integer matrix */

int **cpmat(int **mat, int n, int m){
  int i,j;
  int **tmp;
  tmp=imatrix(0, n-1, 0, m-1);
  for(i=0;i<n;i++){
    for(j=0;j<m;j++){tmp[i][j]=mat[i][j];}
  }
  return tmp;
}

/* change an integer matrix into a subscript matrix */
/* columns first and then rows */

int **mattoind(int **mat, int n, int m){
  int i,j;
  int **tmp;
  int ind=0;
  tmp=imatrix(0, n-1, 0, m-1);
  for(j=0;j<m;j++){
    for(i=0;i<n;i++){
      if(mat[i][j] > 0){
	ind++;tmp[i][j] = ind;
      }
      else if(mat[i][j] < 0){
	ind++;tmp[i][j] = -ind;
      }
      else
	tmp[i][j]=0;
    }
  }
  return tmp;
}

int pairsareequal(int *p1, int *p2, int s1, int *q1, int *q2, int s2){
  int i;
  for (i=0;i<s1;i++){
    if (p1[i]!=p2[i])
      return 0;
  }
  for (i=0;i<s2;i++){
    if (q1[i]!=q2[i])
      return 0;
  }
  return 1;
}



int addtwointstr(int k, int *s, int *s1, int len1, int len2, int ***t)
     /* routine to check if an integer string already belongs to a list of strings, and if not to add it to the end. Returns new free position in the list. */
{
  int i=0,j=0, flag;

  if(k==0){
    (*t) = (int**) malloc(sizeof(int*) * 1);
    (*t)[k] = (int*) malloc(sizeof(int) * (len1+len2+2));
    (*t)[k][0] = len1;
    (*t)[k][1] = len2;
    for(i=2;i<2+len1;i++)
      (*t)[k][i] = s[i-2];
    for(i=2+len1;i<2+len1+len2;i++)
      (*t)[k][i] = s1[i-len1-2];
    return k+1;
  }

  flag=0;
  while(j<k){
    //    printvec(s, len1); printvec(s1, len2);printvec((*t)[j], (*t)[j][0]+(*t)[j][1]+2);
    flag=1;
    if(len1!=(*t)[j][0] || len2!=(*t)[j][1]){
      j++;flag=0;continue;
    }
    for(i=2;i<2+len1;i++){
      if((*t)[j][i]!=s[i-2]){
	flag=0;
	break;
      }
    }
    if(flag){
      for(i=2+len1;i<2+len1+len2;i++){
	if((*t)[j][i]!=s1[i-len1-2]){
	  flag=0;
	  break;
	}
      }
    }
    if(flag) // found
      break;
    j++;
  }



  if(!flag){//not found
    //    fprintf(stderr, "not found: %d\n", k);
    //    printvec(s, len1);printvec(s1,len2);
    (*t)=(int**) realloc((*t), sizeof(int*) *(k+1));
    (*t)[k] = (int*) malloc(sizeof(int) * (len1+len2+2));
    (*t)[k][0] = len1;
    (*t)[k][1] = len2;
    for(i=2;i<2+len1;i++)
      (*t)[k][i] = s[i-2];
    for(i=2+len1;i<2+len1+len2;i++)
      (*t)[k][i] = s1[i-len1-2];
    return k+1;
  }
  return k;
}



int addtwosignedintstr(int k, int *s, int *s1, int len1, int len2, int sgn, int ***t)
     /* routine to check if an integer string already belongs to a list of strings, and if not to add it to the end. Returns new free position in the list. */
{
  int i=0,j=0, flag;

  if(k==0){
    (*t) = (int**) malloc(sizeof(int*) * 1);
    (*t)[k] = (int*) malloc(sizeof(int) * (len1+len2+3));
    (*t)[k][0] = sgn;
    (*t)[k][1] = len1;
    (*t)[k][2] = len2;
    for(i=3;i<3+len1;i++)
      (*t)[k][i] = s[i-3];
    for(i=3+len1;i<3+len1+len2;i++)
      (*t)[k][i] = s1[i-len1-3];
    return k+1;
  }

  flag=0;
  while(j<k){
    //    printvec(s, len1); printvec(s1, len2);printvec((*t)[j], (*t)[j][0]+(*t)[j][1]+2);
    flag=1;
    if(len1!=(*t)[j][0] || len2!=(*t)[j][1]){
      j++;flag=0;continue;
    }
    for(i=3;i<3+len1;i++){
      if((*t)[j][i]!=s[i-3]){
	flag=0;
	break;
      }
    }
    if(flag){
      for(i=3+len1;i<3+len1+len2;i++){
	if((*t)[j][i]!=s1[i-len1-3]){
	  flag=0;
	  break;
	}
      }
    }
    if(flag) // found
      break;
    j++;
  }



  if(!flag){//not found
    //    fprintf(stderr, "not found: %d\n", k);
    //    printvec(s, len1);printvec(s1,len2);
    (*t)=(int**) realloc((*t), sizeof(int*) *(k+1));
    (*t)[k] = (int*) malloc(sizeof(int) * (len1+len2+3));
    (*t)[k][0] = sgn;
    (*t)[k][1] = len1;
    (*t)[k][2] = len2;
    for(i=3;i<3+len1;i++)
      (*t)[k][i] = s[i-3];
    for(i=3+len1;i<3+len1+len2;i++)
      (*t)[k][i] = s1[i-len1-3];
    return k+1;
  }
  return k;
}

// remove vertices of degree 1 and their unique neighbour.
// do this recursively until there is no change
// The matrix is assumed to be a submatrix of the original, described
// by row-set rw and column set col, and the output is a reduced
// row set rwout and a reduced column set colout

int cuttails(int **mat, int n, int m, int *rw, int *col, int **rwout, int **colout){
  int i, j, i1, j1;
  int rwcount, colcount, val;
  int ch=1;
  int rwtmp[n];
  int coltmp[m];
  int del=0;

  for(i=0;i<n;i++)
    rwtmp[i]=rw[i];
  for(j=0;j<m;j++)
    coltmp[j]=col[j];


  while(ch){
    ch=0;
    for(i=0;i<n;i++){
      if(rwtmp[i]>=0){
	rwcount=0;
	for(j=0;j<m;j++){ // count entries in row
	  if(coltmp[j]>=0){
	    if(mat[rw[i]][col[j]]!=0){rwcount++;val=j;}
	  }
	}
	if(rwcount==1){rwtmp[i]=-1;coltmp[val]=-1;ch=1;del++;/* fprintf(stderr, "cut: %d, %d\n", i, val); */} //singleton
      }
    }
    for(j=0;j<m;j++){
      if(coltmp[j]>=0){
	colcount=0;
	for(i=0;i<n;i++){ // count entries in column
	  if(rwtmp[i]>=0){
	    if(mat[rw[i]][col[j]]!=0){colcount++;val=i;}
	  }
	}
	if(colcount==1){rwtmp[val]=-1;coltmp[j]=-1;ch=1;del++;/* fprintf(stderr, "cut1: %d, %d\n", val, j); */}//singleton
      }
    }
  }
  //  fprintf(stderr, "del = %d\n", del);

  (*rwout)=(int *)malloc((size_t) ((n-del)*sizeof(int)));
  (*colout)=(int *)malloc((size_t) ((m-del)*sizeof(int)));
  i1=0;j1=0;
  for(i=0;i<n;i++){
    if(rwtmp[i]>=0){
      (*rwout)[i1]=rw[i];i1++;
    }
  }
  //  fprintf(stderr, "i1 = %d\n", i1);
  for(j=0;j<m;j++){
    if(coltmp[j]>=0){
      (*colout)[j1]=col[j];j1++;
    }
  }
  //  fprintf(stderr, "j1 = %d\n", j1);
  return del;

}

// remove rows/columns with one or fewer nonzero entries

int **simpmat(int **mat, int n, int m, int *n1, int *m1){
  int i,j,k,l,flag;
  int **tmp;
  int x[n], y[m];
  *n1=0;*m1=0;
  for(i=0;i<n;i++){
    flag=0;
    for(j=0;j<m;j++){ // count entries in row
      if(mat[i][j]!=0){flag++;}
    }
    if(flag>=2){x[i]=1;(*n1)++;}
    else{x[i]=0;}
  }

  for(j=0;j<m;j++){
    flag=0;
    for(i=0;i<n;i++){ // count entries in column
      if(mat[i][j]!=0){flag++;}
    }
    if(flag>=2){y[j]=1;(*m1)++;}
    else{y[j]=0;}
  }
  tmp=imatrix(0, (*n1)-1, 0, (*m1)-1);

  k=0;
  for(i=0;i<n;i++){
    if(x[i]>0){
      l=0;
      for(j=0;j<m;j++){
	if(y[j] > 0){
	  tmp[k][l] = mat[i][j];
	  l++;
	}
      }
      k++;
    }
  }

  return tmp;
}


// remove empty rows/columns with one or fewer nonzero entries

int **vsimpmat(int **mat, int n, int m, int *n1, int *m1){
  int i,j,k,l,flag;
  int **tmp;
  int x[n], y[m];
  *n1=0;*m1=0;
  for(i=0;i<n;i++){
    flag=0;
    for(j=0;j<m;j++){ // count entries in row
      if(mat[i][j]!=0){flag++;}
    }
    if(flag>=1){x[i]=1;(*n1)++;}
    else{x[i]=0;}
  }

  for(j=0;j<m;j++){
    flag=0;
    for(i=0;i<n;i++){ // count entries in column
      if(mat[i][j]!=0){flag++;}
    }
    if(flag>=1){y[j]=1;(*m1)++;}
    else{y[j]=0;}
  }
  tmp=imatrix(0, (*n1)-1, 0, (*m1)-1);

  k=0;
  for(i=0;i<n;i++){
    if(x[i]>0){
      l=0;
      for(j=0;j<m;j++){
	if(y[j] > 0){
	  tmp[k][l] = mat[i][j];
	  l++;
	}
      }
      k++;
    }
  }

  return tmp;
}

int **submatminus(int **mat, int n, int m){
  int i,j;
  int **tmp;
  tmp=imatrix(0, n-1, 0, m-1);
  for(i=0;i<n;i++){
    for(j=0;j<m;j++){
      if (mat[i][j]>0)
	tmp[i][j]=0;
      else
	tmp[i][j]=mat[i][j];
    }
  }
  return tmp;
}

int **submatgen(int **imat, int n, int m, int *i1, int *i2, int dim){
  int **tmp;
  int i,j;
  tmp=imatrix(0, dim-1, 0, dim-1);
  for(i=0;i<dim;i++){
    for(j=0;j<dim;j++){
      tmp[i][j]=imat[i1[i]][i2[j]];
    }
  }
  return tmp;
}


int det(int **imat, int n){
  int i;
  int tot=0;
  int **tmpmat;
  if(n==1)
    return imat[0][0];
  for (i=0;i<n;i++){
    tmpmat=submat1(imat, n, 0, i);
    tot+=par(i)*imat[0][i]*det(tmpmat, n-1);
    free_imatrix(tmpmat, 0, n-1, 0, n-1);
  }
  return tot;

}

int detsubmat(int **imat, int n, int m, int *i1, int *i2, int dim){
  int i, j, r;
  int tot=0;
  int i2a[dim-1];
  if(dim==1)
    return imat[i1[0]][i2[0]];

  //  i2a=(int*) malloc(sizeof(int) * (dim-1));
  for (i=0;i<dim;i++){
    if(imat[i1[0]][i2[i]]!=0){
      r=0;
      for(j=0;j<dim;j++){
	if(j!=i)
	  i2a[r++]=i2[j];
      }
      tot+=par(i)*imat[i1[0]][i2[i]]*detsubmat(imat, n, m, i1+1, i2a, dim-1);
    }

  }
  //  free((FREE_ARG) (i2a));
  return tot;

}



// the symbolic determinant of a symbolic matrix

ex symbdetsubmat(matrix imat, int n, int m, int *i1, int *i2, int dim){
  int i, j, r;
  ex tot=0;
  //ex extmp;
  int i2a[dim-1];
  if(dim==1)
    return imat(i1[0], i2[0]);
  /* if(dim>4) */
  /*   cout << "dim= " << dim << endl; */
  for (i=0;i<dim;i++){

    if(imat(i1[0], i2[i])!=0){
      r=0;
      for(j=0;j<dim;j++){
	if(j!=i)
	  i2a[r++]=i2[j];
      }
      tot+=par(i)*imat(i1[0],i2[i])*symbdetsubmat(imat, n, m, i1+1, i2a, dim-1);
    /*   extmp=symbdetsubmat(imat, n, m, i1+1, i2a, dim-1); */
    /*   if(!extmp.info(info_flags::positive)) */
    /* 	tot+=par(i)*imat(i1[0],i2[i])*expand(extmp); */
    /*   else{ */
    /* 	tot+=par(i)*imat(i1[0],i2[i])*extmp; */
    /* 	//cout << "signed expr:\n" << extmp << endl; */
    /*   } */
    /* 	//      tot+=par(i)*imat(i1[0],i2[i])*symbdetsubmat(imat, n, m, i1+1, i2a, dim-1); */
    /*   //      cout << tot << endl; */
    /* } */

    }
  }

  return expand(tot);

}

ex mydet(matrix imat, int n){
int *xc=(int *)malloc((size_t) (n*sizeof(int)));
 ex detex;
 firstcomb(xc, n, n);
 detex=symbdetsubmat(imat,n,n,xc,xc,n);
 free((char*)xc);
 return detex;
}

// the symbolic determinant of a symbolic matrix
// without expansion of the answer

ex symbdetsubmat_noexp(matrix imat, int n, int m, int *i1, int *i2, int dim){
  int i, j, r;
  ex tot=0;
  //ex extmp;
  int i2a[dim-1];
  if(dim==1)
    return imat(i1[0], i2[0]);
  /* if(dim>4) */
  /*   cout << "dim= " << dim << endl; */
  for (i=0;i<dim;i++){

    if(imat(i1[0], i2[i])!=0){
      r=0;
      for(j=0;j<dim;j++){
	if(j!=i)
	  i2a[r++]=i2[j];
      }
      tot+=par(i)*imat(i1[0],i2[i])*symbdetsubmat(imat, n, m, i1+1, i2a, dim-1);
    /*   extmp=symbdetsubmat(imat, n, m, i1+1, i2a, dim-1); */
    /*   if(!extmp.info(info_flags::positive)) */
    /* 	tot+=par(i)*imat(i1[0],i2[i])*expand(extmp); */
    /*   else{ */
    /* 	tot+=par(i)*imat(i1[0],i2[i])*extmp; */
    /* 	//cout << "signed expr:\n" << extmp << endl; */
    /*   } */
    /* 	//      tot+=par(i)*imat(i1[0],i2[i])*symbdetsubmat(imat, n, m, i1+1, i2a, dim-1); */
    /*   //      cout << tot << endl; */
    /* } */

    }
  }

  return tot;

}

//
// monomial to a list
//

lst monotolist(ex e){
  lst tmp;
  int j;
  int coe=-1;
  //int r;

 if(is_a<add>(e)){
    cout << "failure in monotolist: not expecting a sum:\n";
    cout << e << endl;
    exit(0);
  }

  if(e!=0){
    for (size_t k=0; k!=e.nops(); ++k){
      if(is_a<power>(e.op(k))){
	for(j=0;j<e.op(k).op(1);j++){
	  if(is_a<symbol>(e.op(k).op(0)))
	    tmp.append(e.op(k).op(0));
	  else{
	    cout << "failure 1 in monotolist: expecting a symbol but got:\n";
	    cout << e.op(k).op(0) << endl;
	    exit(0);
	  }

	  //	   cout << e.op(k).op(0) << endl;
	}
      }
      else if(is_a<numeric>(e.op(k))){
	coe=k;
      }
      else{
	if(is_a<symbol>(e.op(k)))
	  tmp.append(e.op(k));
	else{
	  cout << "failure 2 in monotolist: expecting a symbol but got:\n";
	  cout << e.op(k) << endl;
	  exit(0);
	}
	//	 cout << e.op(k) << endl;
      }
    }
  }
  /* if(coe >=0) */
  /*   tmp.append(coe); */
  /* else */
  /*   tmp.append(1); */

  return tmp;
}

void iswap(int *v, int i, int j){
  int temp = v[i];
  v[i]=v[j];
  v[j]=temp;
}





void qsortt(int *v, int left, int right)
{
  int i, last;

  if (left >= right)
    return;
  iswap(v,left,(left + right)/2);
  last = left;
  for (i=left+1;i<=right;i++)
    if (v[i]<v[left])
      iswap(v,++last,i);
  iswap(v,left,last);
  qsortt(v, left, last-1);
  qsortt(v, last+1, right);
}

//
// swap ith and jth columns of matrix *v
// n is the column length
//

void colswap(int ***v, int n, int i, int j){
  int k;
  int temp;
  for(k=0;k<n;k++){
    temp=(*v)[k][i];
    (*v)[k][i]=(*v)[k][j];
    (*v)[k][j]=temp;
  }
  return;
}

void colswapb(bool ***v, int n, int i, int j){
  int k;
  bool temp;
  for(k=0;k<n;k++){
    temp=(*v)[k][i];
    (*v)[k][i]=(*v)[k][j];
    (*v)[k][j]=temp;
  }
  return;
}

void colswapmat(matrix *v, int n, int i, int j){
  int k;
  ex temp;
  for(k=0;k<n;k++){
    temp=(*v)(k,i);
    (*v)(k,i)=(*v)(k,j);
    (*v)(k,j)=temp;
  }
  return;
}


int cmpare(int *a, int *b, int n){
  int i;
  for(i=0;i<n;i++){
    if(a[i]<b[i])
      return -1;
    else if(a[i]>b[i])
      return 1;
  }
  return 0;
}



//
// reverse lexicographic ordering
//

int cmparerev(int *a, int *b, int n){
  int i;
  for(i=n-1;i>=0;i--){
    if(a[i]<b[i])
      return -1;
    else if(a[i]>b[i])
      return 1;
  }
  return 0;
}

int cmparebrev(bool *a, bool *b, int n){
  int i;
  for(i=n-1;i>=0;i--){
    if(a[i]<b[i])
      return -1;
    else if(a[i]>b[i])
      return 1;
  }
  return 0;
}

// n is the column length

void lexsortcol(int ***v, int n, int left, int right)
{
  int i, last;

  if (left >= right)
    return;
  colswap(v,n,left,(left + right)/2);
  last = left;
  for (i=left+1; i<=right;i++)
    if (cmparerev((*v)[i],(*v)[left],n)<0)
      colswap(v,n,++last,i);
  colswap(v,n,left,last);
  lexsortcol(v,n,left,last-1);
  lexsortcol(v,n,last+1,right);
  return;
}

void lexsortcol4(bool ***v, bool ***v1, int ***imat, matrix *mmat, int n, int left, int right)
{
  int i,k,last;
  bool vleft[n];
  bool vi[n];

  if (left >= right)
    return;
  colswapb(v,n,left,(left + right)/2);
  colswapb(v1,n,left,(left + right)/2);
  colswap(imat,n,left,(left + right)/2);
  colswapmat(mmat,n,left,(left + right)/2);
  last = left;


  for (i=left+1; i<=right;i++){
    for(k=0;k<n;k++)
      vi[k]=(*v)[k][i];
  for(k=0;k<n;k++)
    vleft[k]=(*v)[k][left];
    if (cmparebrev(vi,vleft,n)<0){
      last++;
      colswapb(v,n,last,i);
      colswapb(v1,n,last,i);
      colswap(imat,n,last,i);
      colswapmat(mmat,n,last,i);
    }
  }

  colswapb(v,n,left,last);
  colswapb(v1,n,left,last);
  colswap(imat,n,left,last);
  colswapmat(mmat,n,left,last);

  lexsortcol4(v,v1,imat,mmat,n,left,last-1);
  lexsortcol4(v,v1,imat,mmat,n,last+1,right);
  return;
}

// Reorder the columns of a bool matrix, an int matrix and 
// an ex matrix putting them in upwards lexicographic order

void multireordercollex(bool **bmatin, bool ***bmatout, bool **bmat2in, bool ***bmat2out, int **imatin, int ***imatout, matrix mmatin, matrix *mmatout, int n, int m){
  int i,j;
  (*bmatout)=bmatrix(0, n-1, 0, m-1);
  (*bmat2out)=bmatrix(0, n-1, 0, m-1);
  (*imatout)=imatrix(0, n-1, 0, m-1);
  (*mmatout)=matrix(n,m);
  //copy
  for(i=0;i<n;i++){
    for(j=0;j<m;j++){
      (*bmatout)[i][j]=bmatin[i][j];
      (*bmat2out)[i][j]=bmat2in[i][j];
      (*imatout)[i][j]=imatin[i][j];
      (*mmatout)(i,j)=mmatin(i,j);
    }
  }
  lexsortcol4(bmatout,bmat2out,imatout,mmatout,n,0,m-1);
  return;
}



//
// convert an expression, assumed to be a monomial
// into a list of integers. Each symbol in the monomial is 
// assumed to be of the form vk, where "v" is a letter and "k"
// an integer. For example 3*x2^2*x4 becomes [2,2,4]
// The coefficient is discarded.
//



int monotointlist(ex e, int **lst, int *r){
  
  int j;
  int ret=1;
  string str;
  *r=0;
  (*lst)=NULL; // NULL returned if expression is zero

  if(is_a<add>(e)){
    cout << "failure in monotointlist: not expecting a sum:\n";
    cout << e << endl;
    exit(0);
  }

  if(is_a<power>(e)){
    for(j=0;j<e.op(1);j++){
      (*r)++;
      if((*r)==1)
	(*lst)=(int *)malloc((size_t) (sizeof(int)));
      else
	(*lst)=(int*) realloc((*lst), sizeof(int)*((*r)));
      //	   cout << e.op(k).op(0) << endl;
      if(is_a<symbol>(e.op(0))){
	str=(ex_to<symbol>(e.op(0))).get_name();
	(*lst)[(*r)-1] = atoi(str.c_str()+1);
      }
      else{
	cout << "failure 1 in monotointlist: expecting a symbol but got:\n";
	cout << e.op(0) << endl;
	exit(0);
      }
	 
    }
  }
  else if(is_a<symbol>(e)){ // single symbol - not a product
    str=(ex_to<symbol>(e)).get_name();
    (*r)++;
    if((*r)==1)
      (*lst)=(int *)malloc((size_t) (sizeof(int)));
    else
      (*lst)=(int*) realloc((*lst), sizeof(int)*((*r)));
    (*lst)[(*r)-1] = atoi(str.c_str()+1);
  }
  else if(e!=0){
    for (size_t k=0; k!=e.nops(); ++k){
      if(is_a<power>(e.op(k))){
	for(j=0;j<e.op(k).op(1);j++){
	  (*r)++;
	  if((*r)==1)
	    (*lst)=(int *)malloc((size_t) (sizeof(int)));
	  else
	    (*lst)=(int*) realloc((*lst), sizeof(int)*((*r)));
	  //	   cout << e.op(k).op(0) << endl;
	  if(is_a<symbol>(e.op(k).op(0))){
	    str=(ex_to<symbol>(e.op(k).op(0))).get_name();
	    (*lst)[(*r)-1] = atoi(str.c_str()+1);
	  }
	  else{
	    cout << "failure 1 in monotointlist: expecting a symbol but got:\n";
	    cout << e.op(k).op(0) << endl;
	    exit(0);
	  }
	 
	}
      }
      else if(is_a<numeric>(e.op(k))){
	ret=ex_to<numeric>(e.op(k)).to_int();
      }// discard numbers
      else{
	(*r)++;
	if((*r)==1)
	  (*lst)=(int *)malloc((size_t) (sizeof(int)));
	else
	  (*lst)=(int*) realloc((*lst), sizeof(int) *((*r)));
	//	   cout << e.op(k).op(0) << endl;
	if(is_a<symbol>(e.op(k))){
	  str=(ex_to<symbol>(e.op(k))).get_name();
	  (*lst)[(*r)-1] = atoi(str.c_str()+1);
	}
	else{
	  cout << "failure 2 in monotointlist: expecting a symbol but got:\n";
	  cout << e.op(k) << endl;
	  exit(0);
	}

      }
    }
  }
  qsortt((*lst),0,(*r)-1);
  return ret;
}

// Does not assume prior knowledge of the order
// the first entry in lst is the order of the monomial

int monotointlist0(ex e, int **lst){
  int i,p,r,*q;
  p=monotointlist(e,&q,&r);
  (*lst)=(int *)malloc((size_t) (sizeof(int)*(r+1)));
  (*lst)[0]=r;
  for(i=1;i<=r;i++)
    (*lst)[i]=q[i-1];
  free((char*)q);

  return p;
}

// r is the number of monomials
// the return is the degree. (assumed homogeneous)
// the switch hom tells us whether to expect a homogeneous polynomial

int polytointlist(ex e, int ***lst, int **cflst, long *r){
  int monord,ord1;
  if(e==0){// a single term
    (*r)=1;
    (*lst)=(int **)malloc((size_t) (sizeof(int*)));
    (*lst)[0]=(int *)malloc((size_t) (sizeof(int)));
    (*cflst)=(int *)malloc((size_t) (sizeof(int)));
    (*cflst)[0]=0;
    (*lst)[0][0]=0;
    monord=0;
    return monord;
  }
  if(!(is_a<add>(e))){// a single term
    (*r)=1;
    (*lst)=(int **)malloc((size_t) (sizeof(int*)));
    (*cflst)=(int *)malloc((size_t) (sizeof(int)));
    (*cflst)[0]=monotointlist(e,lst[0],&monord); //allocates
    return monord;
  }
  (*r)=e.nops(); // terms in sum
  //cout << "r= "<< *r << endl;
  (*lst)=(int **)malloc((size_t) (sizeof(int*)*(*r)));
  (*cflst)=(int *)malloc((size_t) (sizeof(int)*(*r)));
  for (size_t k=0; k!=e.nops(); ++k){
    //cout << "k= "<< k << endl;
    //cout << "first term= "<< e.op(k) << endl;
    (*cflst)[k]=monotointlist(e.op(k),&((*lst)[k]),&ord1);
    //cout << "firstcoeff= "<< (*cflst)[k] << endl;
    if(k==0)
      monord=ord1;
    else if(ord1!=monord){
      cout << "polynomial is not homogeneous" << endl;
      cout << "expected ord = " << monord << " term = " << e.op(k) <<" actual ord = " << ord1 << endl;
      exit(0);
    }
  }
  return monord;

}

// r is the number of monomials
// the return is the degree. 
// not assumed homogeneous. The first entry in each 
// output vector is the degree of the the monomial

void polytointlist0(ex e, int ***lst, int **cflst, long *r){
  if(e==0){// empty
    (*r)=1;
    (*lst)=(int **)malloc((size_t) (sizeof(int*)));
    (*lst)[0]=(int *)malloc((size_t) (sizeof(int)));
    (*cflst)=(int *)malloc((size_t) (sizeof(int)));
    (*cflst)[0]=0;
    (*lst)[0][0]=0;
    return;
  }
  if(!(is_a<add>(e))){// not a sum
    (*r)=1;
    (*lst)=(int **)malloc((size_t) (sizeof(int*)));
    (*cflst)=(int *)malloc((size_t) (sizeof(int)));
    (*cflst)[0]=monotointlist0(e,&((*lst)[0])); //allocates
    return;
  }
  (*r)=e.nops(); // terms in sum
  //cout << "r= "<< *r << endl;
  (*lst)=(int **)malloc((size_t) (sizeof(int*)*(*r)));
  (*cflst)=(int *)malloc((size_t) (sizeof(int)*(*r)));
  for (size_t k=0; k!=e.nops(); ++k){
    (*cflst)[k]=monotointlist0(e.op(k),&((*lst)[k]));
  }
  return;
}


//
// to convert a symbolic matrix to an integer matrix
//

int sumtointlist(ex e, int **lst){
  
  string str;
  (*lst)=NULL; // NULL returned if expression is zero

  if(e==0){
    (*lst)=(int *)malloc((size_t) (sizeof(int)));
    (*lst)[0]=0;
  }
  else if(is_a<add>(e)){ // a sum
    (*lst)=(int *)malloc((size_t) (sizeof(int)*(e.nops()+1)));
    (*lst)[0] = e.nops();
    for (size_t k=0; k!=e.nops(); ++k){
      if(is_a<symbol>(e.op(k))){
	str=(ex_to<symbol>(e.op(k))).get_name();
	(*lst)[k+1] = atoi(str.c_str()+1);
      }
      else if(is_a<symbol>(-(e.op(k)))){
	str=(ex_to<symbol>(-(e.op(k)))).get_name();
	(*lst)[k+1] = -atoi(str.c_str()+1);
      }
      else{ // something unexpected in the sum
	cout << "failure in sumtointlist: not expecting this:\n";
	cout << e.op(k) << endl;
	exit(0);
      }
    }
    return e.nops()+1;
  }
  else if(is_a<symbol>(e)){ // a positive symbol
    (*lst)=(int *)malloc((size_t) (sizeof(int)*(2)));
    (*lst)[0] = 1;
    str=(ex_to<symbol>(e)).get_name();
    (*lst)[1] = atoi(str.c_str()+1);
    return 2;
  }
  else if(is_a<symbol>(-e)){ // a negative symbol
    (*lst)=(int *)malloc((size_t) (sizeof(int)*(2)));
    (*lst)[0] = 1;
    str=(ex_to<symbol>(-e)).get_name();
    (*lst)[1] = -atoi(str.c_str()+1);
    return 2;
  }
  else{ // something else
    cout << "failure 2 in sumtointlist: not expecting this:\n";
    cout << e << endl;
    exit(0);
  }
 
  return 0;
}

void checkmonotolist(){
  lst l;
  symbol x("x"), y("y");
  ex e=2*x*x*pow(y, 2)-pow(y,3);
  l=monotolist(e);
  cout << l << endl;
  return;
}

void checkmonotointlist(){
  int *l;
  int r, k;
  symbol x1("x1"), x2("x2");
  ex e=-32*x1*x1*pow(x2, 5);
  r=monotointlist(e, &l, &k);
  printvec(l, k);
  free ((char *)l);
  return;
}


//
// extract the coefficient of a term from the determinant of a 
// matrix. The routine is very inefficient: 
// for a dense matrix, can be almost as bad as computing a determinant
//

int coeff(matrix J2, int *i1, int *i2, int dim, lst trm){
  int tot=0, j, l, m,r;
  int i2a[dim-1];
  lst trm1;
  int sgn;
  if(dim>12)
  cout<<"calling: dim = " << dim << endl;

  if(dim==1){
    //fprintf(stderr, "dim=1\n");
    if(is_a<add>(J2(i1[0],i2[0]))){ // a sum of terms
      //fprintf(stderr, "multiterm\n");
      for (size_t k=0; k!=J2(i1[0],i2[0]).nops(); ++k){
	if(trm[0]==J2(i1[0],i2[0]).op(k)){
	  //cout << "found" << endl;
	  return 1;
	}
	else if(trm[0]==-J2(i1[0],i2[0]).op(k))
	  return -1;
      }
    }
    else{ // a single term
      //fprintf(stderr, "single termm\n");
      if(trm[0]==J2(i1[0],i2[0]))
	return 1;
      else if(trm[0]==-J2(i1[0],i2[0]))
	return -1;
    }
    return 0;
  }

  //fprintf(stderr, "dim=%d\n", dim);
  for(j=0;j<dim;j++){
    if(J2(i1[0],i2[j])==0){}
    else if(is_a<add>(J2(i1[0],i2[j]))){
      //fprintf(stderr, "multiterm\n");
      for (size_t k = 0; k != J2(i1[0],i2[j]).nops(); ++k){
      	//cout << "k=" <<k<< endl;
      	l=0;sgn=0;
      	while(l<dim && sgn==0){
      	  if(trm[l]==J2(i1[0],i2[j]).op(k))
      	    sgn=1;
      	  else if(trm[l]==-J2(i1[0],i2[j]).op(k))
      	    sgn=-1;
	  if(sgn){
	    r=0;
	    for(m=0;m<dim;m++){
	      if(m!=j)
		i2a[r++]=i2[m];
	    }
	    trm1=lst();
	    for(m=0;m<dim;m++){
	      if(m!=l)
		trm1.append(trm[m]);
	    }
	    //cout << "trm1 = " << trm1 << endl;
	    tot+=par(j)*sgn*coeff(J2, i1+1, i2a, dim-1, trm1);
	    break; // found in list
	  }
	  l++;
	}
      }
    }
    else{
      //fprintf(stderr, "single term\n");
      l=0;sgn=0;
      while(l<dim && sgn==0){
	//cout << "trm[l] = " << trm[l] << endl;
	//cout << "J2(i1[0],i2[j]) = " << J2(i1[0],i2[j]) << endl;
	if(trm[l]==J2(i1[0],i2[j])){
	  //cout << "found_single\n";
	  sgn=1;
	}
	else if(trm[l]==-J2(i1[0],i2[j]))
	  sgn=-1;
	if(sgn){// lth member found
	  r=0;
	  for(m=0;m<dim;m++){
	    if(m!=j)
	      i2a[r++]=i2[m];
	  }
	  trm1=lst();
	  for(m=0;m<dim;m++){
	    if(m!=l)
	      trm1.append(trm[m]);
	  }
	  //printvec(i2a,dim-1);
	  //cout << "trm1 = " << trm1 << endl;
	  tot+=par(j)*sgn*coeff(J2, i1+1, i2a, dim-1, trm1);
	  break; // found in list
	}
	l++;
      }

    }

  }
  return tot;
}

// recursive function to:
// extract the coefficient of a term from the determinant of a matrix
// each entry in the matrix is a list of integers
// the first entry in each list, is the length
// subsequent entries can be seen as symbols

int coeff_int(int ****J2, int *i1, int *i2, int dim, int *trm){
  int tot=0, k, j, l, m,r;
  int i2a[dim-1];
  int trm1[dim-1];
  int sgn;
  //cout<<"calling\n";

  if(dim==1){
    //cout << "in here\n";
    for (k=1; k<=(*J2)[i1[0]][i2[0]][0]; k++){
      if(trm[0]==(*J2)[i1[0]][i2[0]][k]){
	//cout << "found" << endl;
	return 1;
      }
      else if(trm[0]==-(*J2)[i1[0]][i2[0]][k])
	return -1;
    }
    return 0;
  }


  //fprintf(stderr, "dim=%d\n", dim);
  for(j=0;j<dim;j++){
    //cout << "in here1 dim=" << dim << endl;
    if((*J2)[i1[0]][i2[j]][0]==0){}
    else{
      //cout << "i,j " << i1[0] << "," << i2[j] << endl;
      for (k=1; k<=(*J2)[i1[0]][i2[j]][0]; k++){
	//cout << trm[0] << endl;
	//cout << "i,j,k= " << i1[0] << "," << i2[j] << "," << k << endl;
	//cout << "in here2 (*J2)[i1[0]][i2[j]][k]=" << (*J2)[i1[0]][i2[j]][k] << endl;
      	//cout << "k=" <<k<< endl;
      	l=0;sgn=0;
      	while(l<dim && sgn==0){
      	  if(trm[l]==(*J2)[i1[0]][i2[j]][k])
      	    sgn=1;
      	  else if(trm[l]==-(*J2)[i1[0]][i2[j]][k])
      	    sgn=-1;
	  if(sgn){
	    r=0;
	    for(m=0;m<dim;m++){
	      if(m!=j)
		i2a[r++]=i2[m];
	    }
	    r=0;
	    for(m=0;m<dim;m++){
	      if(m!=l)
		trm1[r++]=trm[m];
	    }

	    //cout << "(*J2)[i1[0]][i2[j]][k] =" << (*J2)[i1[0]][i2[j]][k] << endl;
	    //cout << "trm[l] =" << trm[l] << endl;
	    //cout << "trm1 = " << trm1 << endl;
	    tot+=par(j)*sgn*coeff_int(J2, i1+1, i2a, dim-1, trm1);
	    break; // found in list
	  }
	  l++;
	}
      }
    }
  }
  return tot;
}

// routine to check the coeff function on a large matrix 
// whose determinant has been evaluated with GINAC

void check_coeff_fn(){
  int k=15;
  matrix trial(k,k);
  int *xc, *yc;
  lst trm;
  possymbol v1("v1"), v2("v2"), v3("v3"), v4("v4"), v5("v5"), v6("v6"), v7("v7"), v8("v8"), v9("v9"), v10("v10"), v11("v11");
  trial(0,0)=v2+v1+v3; trial(0,1)=-v4; trial(0,2)=v6; trial(0,5)=v6;
  trial(1,0)=-v3;trial(1,1)= v4+v5+v1;trial(1,2)=-v7;trial(1,5)=-v2;trial(1,8)=v6;
  trial(2,0)=v2;trial(2,1)=-v5;trial(2,2)=v7+v8+v1+v6;trial(2,3)=-v9;trial(2,6)=-v2;
  trial(3,2)=-v8;trial(3,3)=v10+v1+v9;trial(3,4)=-v11;trial(3,7)=-v2;trial(3,12)=-v6;
  trial(4,3)=-v10;trial(4,4)=v11+v1;trial(4,8)=-v2;trial(4,13)= -v6;
  trial(5,1)=-v1;trial(5,5)=v4+v2+v5+v3;trial(5,6)=-v7;trial(5,9)=-v6;
  trial(6,0)=v1;trial(6,2)=-v1;trial(6,5)=-v5;trial(6,6)=v7+v2+v8+v6+v3;trial(6,7)=-v9;trial(6,9)=-v4;
  trial(7,3)=-v1;trial(7,6)=-v8;trial(7,7)=v10+v2+v9+v3;trial(7,8)=-v11;trial(7,10)=-v4;trial(7,12)=v6;
  trial(8,4)=-v1;trial(8,7)=-v10;trial(8,8)=v2+v11+v3;trial(8,11)=-v4;trial(8,13)=v6;
  trial(9,1)=v1;trial(9,5)=-v2;trial(9,6)=-v3;trial(9,9)=v7+v4+v8+v5+v6;trial(9,10)=-v9;
  trial(10,7)=-v3;trial(10,9)=-v8;trial(10,10)=v10+v4+v5+v9;trial(10,11)=-v11;trial(10,12)=-v7;
  trial(11,8)=-v3;trial(11,10)=-v10;trial(11,11)=v4+v11+v5;trial(11,13)=-v7;
  trial(12,3)=-v1;trial(12,7)=v2;trial(12,10)=-v5;trial(12,12)=v10+v7+v8+v9+v6;trial(12,13)=-v11;
  trial(13,4)=-v1;trial(13,8)=v2;trial(13,11)=-v5;trial(13,12)=-v10;trial(13,13)=v7+v11+v8+v6;trial(13,14)=-v9; 
  trial(14,13)=-v8;trial(14,14)=v10+v11+v9;
  printexmat(trial, k, k);
  xc=(int *)malloc((size_t) (k*sizeof(int)));
  yc=(int *)malloc((size_t) (k*sizeof(int)));
  firstcomb(xc, k, k);
  firstcomb(yc, k, k);
  printvec(xc,k);
  printvec(yc,k);
  trm = monotolist(pow(v10,2)*pow(v7,3)*pow(v4,2)*pow(v2,3)*v11*pow(v8,2)*v1*v9);
  cout << trm << endl;
  cout << "coeff = " << coeff(trial, xc, yc, k, trm) << endl;
  trm.remove_all();
  //120*v10^2*v7*v11^2*v8^2*v5^3*v1^2*v3^3
  trm = monotolist(-v10*v10*v7*v11*v11*pow(v8,2)*pow(v5,3)*pow(v1,2)*pow(v3,3));
  //trm = v10,v10,v7,v11,v11,v8,v8,v5,v5,v5,v1,v1,v3,v3,v3; // expect 120
  cout << trm << endl;
  cout << "coeff = " << coeff(trial, xc, yc, k, trm) << endl;
  //48*v7*v4*v2*v11^2*v5^4*v9^2*v6^2*v3^2
  trm.remove_all();
  trm = v7,v4,v2,v11,v11,v5,v5,v5,v5,v9,v9,v6,v6,v3,v3; // expect 48
  cout << "coeff = " << coeff(trial, xc, yc, k, trm) << endl;
  //v8^2*v5^2*v1^3*v9*v3*v10*v7^2*v4*v2*v11
  trm.remove_all();
  trm = v8,v8,v5,v5,v1,v1,v1,v9,v3,v10,v7,v7,v4,v2,v11; // expect 267
  cout << "coeff = " << coeff(trial, xc, yc, k, trm) << endl;
  free((char *) xc);free((char *) yc);

  return;
}

void check_coeff_int(){
  int k=15;
  int i, j, m;
  int ***trial;
  trial =(int ***) malloc((size_t)(k*sizeof(int**)));
  for(i=0;i<k;i++){
    trial[i]=(int **) malloc((size_t)(k*sizeof(int*)));
    for(j=0;j<k;j++){
      trial[i][j]=(int *) malloc((size_t)(6*sizeof(int)));
    }
  }
 
  int *xc, *yc;
  int trm[15]={10,10,7,11,11,8,8,5,5,5,1,1,3,3,3};
  for(i=0;i<k;i++){
    for(j=0;j<k;j++){
      for(m=0;m<6;m++){
	trial[i][j][m]=0;
      }
    }
  }
 
  trial[0][0][0]=3;trial[0][0][1]=1; trial[0][0][2]=2; trial[0][0][3]=3; 
  trial[0][1][0]=1;trial[0][1][1]=-4; 
  trial[0][2][0]=1;trial[0][2][1]=6; 
  trial[0][5][0]=1;trial[0][5][1]=6;

  trial[1][0][0]=1;trial[1][0][1]=-3;
  trial[1][1][0]= 3;trial[1][1][1]= 1;trial[1][1][2]= 4;trial[1][1][3]= 5;
  trial[1][2][0]=1;trial[1][2][1]=-7;
  trial[1][5][0]=1;trial[1][5][1]=-2;
  trial[1][8][0]=1;trial[1][8][1]=6;

  trial[2][0][0]=1;trial[2][0][1]=2;
  trial[2][1][0]=1;trial[2][1][1]=-5;
  trial[2][2][0]=4;trial[2][2][1]=7;trial[2][2][2]=8;trial[2][2][3]=1;trial[2][2][4]=6;
  trial[2][3][0]=1;trial[2][3][1]=-9;
  trial[2][6][0]=1;trial[2][6][1]=-2;

  trial[3][2][0]=1;trial[3][2][1]=-8;
  trial[3][3][0]=3;trial[3][3][1]=10;trial[3][3][2]=1;trial[3][3][3]=9;
  trial[3][4][0]=1;trial[3][4][1]=-11;
  trial[3][7][0]=1;trial[3][7][1]=-2;
  trial[3][12][0]=1;trial[3][12][1]=-6;

  trial[4][3][0]=1;trial[4][3][1]=-10;
  trial[4][4][0]=2;trial[4][4][1]=11;trial[4][4][2]=1;
  trial[4][8][0]=1;trial[4][8][1]=-2;
  trial[4][13][0]=1;trial[4][13][1]= -6;

  trial[5][1][0]=1;trial[5][1][1]=-1;
  trial[5][5][0]=4;trial[5][5][1]=4;trial[5][5][2]=2;trial[5][5][3]=5;trial[5][5][4]=3;
  trial[5][6][0]=1;trial[5][6][1]=-7;
  trial[5][9][0]=1;trial[5][9][1]=-6;

  trial[6][0][0]=1;trial[6][0][1]=1;
  trial[6][2][0]=1;trial[6][2][1]=-1;
  trial[6][5][0]=1;trial[6][5][1]=-5;
  trial[6][6][0]=5;trial[6][6][1]=7;trial[6][6][2]=2;trial[6][6][3]=8;trial[6][6][4]=6;trial[6][6][5]=3;
  trial[6][7][0]=1;trial[6][7][1]=-9;
  trial[6][9][0]=1;trial[6][9][1]=-4;

  trial[7][3][0]=1;trial[7][3][1]=-1;
  trial[7][6][0]=1;trial[7][6][1]=-8;
  trial[7][7][0]=4;trial[7][7][1]=10;trial[7][7][2]=2;trial[7][7][3]=9;trial[7][7][4]=3;
  trial[7][8][0]=1;trial[7][8][1]=-11;
  trial[7][10][0]=1;trial[7][10][1]=-4;
  trial[7][12][0]=1;trial[7][12][1]=6;

  trial[8][4][0]=1;trial[8][4][1]=-1;
  trial[8][7][0]=1;trial[8][7][1]=-10;
  trial[8][8][0]=3;trial[8][8][1]=2;trial[8][8][2]=11;trial[8][8][3]=3;
  trial[8][11][0]=1;trial[8][11][1]=-4;
  trial[8][13][0]=1;trial[8][13][1]=6;

  trial[9][1][0]=1;trial[9][1][1]=1;
  trial[9][5][0]=1;trial[9][5][1]=-2;
  trial[9][6][0]=1;trial[9][6][1]=-3;
  trial[9][9][0]=5;trial[9][9][1]=7;trial[9][9][2]=4;trial[9][9][3]=8;trial[9][9][4]=5;trial[9][9][5]=6;
  trial[9][10][0]=1;trial[9][10][1]=-9;

  trial[10][7][0]=1;trial[10][7][1]=-3;
  trial[10][9][0]=1;trial[10][9][1]=-8;
  trial[10][10][0]=4;trial[10][10][1]=10;trial[10][10][2]=4;trial[10][10][3]=5;trial[10][10][4]=9;
  trial[10][11][0]=1;trial[10][11][1]=-11;
  trial[10][12][0]=1;trial[10][12][1]=-7;

  trial[11][8][0]=1;trial[11][8][1]=-3;
  trial[11][10][0]=1;trial[11][10][1]=-10;
  trial[11][11][0]=3;trial[11][11][1]=4;trial[11][11][2]=11;trial[11][11][3]=5;
  trial[11][13][0]=1;trial[11][13][1]=-7;

  trial[12][3][0]=1;trial[12][3][1]=-1;
  trial[12][7][0]=1;trial[12][7][1]=2;
  trial[12][10][0]=1;trial[12][10][1]=-5;
  trial[12][12][0]=5;trial[12][12][1]=10;trial[12][12][2]=7;trial[12][12][3]=8;trial[12][12][4]=9;trial[12][12][5]=6;
  trial[12][13][0]=1;trial[12][13][1]=-11;

  trial[13][4][0]=1;trial[13][4][1]=-1;
  trial[13][8][0]=1;trial[13][8][1]=2;
  trial[13][11][0]=1;trial[13][11][1]=-5;
  trial[13][12][0]=1;trial[13][12][1]=-10;
  trial[13][13][0]=4;trial[13][13][1]=7;trial[13][13][2]=11;trial[13][13][3]=8;trial[13][13][4]=6;
  trial[13][14][0]=1;trial[13][14][1]=-9;

  trial[14][13][0]=1;trial[14][13][1]=-8;
  trial[14][14][0]=3;trial[14][14][1]=10;trial[14][14][2]=11;trial[14][14][3]=9;

  xc=(int *)malloc((size_t) (k*sizeof(int)));
  yc=(int *)malloc((size_t) (k*sizeof(int)));
  firstcomb(xc, k, k);
  firstcomb(yc, k, k);
  printvec(xc,k);
  printvec(yc,k);
  //120*v10^2*v7*v11^2*v8^2*v5^3*v1^2*v3^3
  //trm={10,10,7,11,11,8,8,5,5,5,1,1,3,3,3};
  //trm = v10,v10,v7,v11,v11,v8,v8,v5,v5,v5,v1,v1,v3,v3,v3; // expect 120
  //cout << *trm << endl;
  cout << "coeff = " << coeff_int(&trial, xc, yc, k, trm) << endl;
  /* //48*v7*v4*v2*v11^2*v5^4*v9^2*v6^2*v3^2 */
  /* trm.remove_all(); */
  /* trm = v7,v4,v2,v11,v11,v5,v5,v5,v5,v9,v9,v6,v6,v3,v3; // expect 48 */
  /* cout << "coeff = " << coeff(trial, xc, yc, k, trm) << endl; */
  /* //v8^2*v5^2*v1^3*v9*v3*v10*v7^2*v4*v2*v11 */
  /* trm.remove_all(); */
  /* trm = v8,v8,v5,v5,v1,v1,v1,v9,v3,v10,v7,v7,v4,v2,v11; // expect 267 */
  /* cout << "coeff = " << coeff(trial, xc, yc, k, trm) << endl; */
  free((char *) xc);free((char *) yc);

  return;
}



// qualitative addition (inputs/outputs are 
// 0,1,-1 and 2 where 2 means unsigned)

int qualp(int a, int b){
  if(a==0 || b==0)
    return a+b;
  if(a==1 && b==1)
    return 1;
  if(a==-1 && b==-1)
    return -1;
  if(a==2 || b==2 || (a*b == -1))
    return 2;
  return a+b;

}


// qualitative multiplication (inputs/outputs are 
// 0,1,-1 and 2 where 2 means unsigned)

int qualt(int a, int b){
  if(a==0 || b==0) // takes priority
    return 0;
  if(a==2 || b==2)
    return 2;
  return a*b;
}

int sgnf(int j){
  if(j==0)
    return 0;
  if(j>0)
    return 1;
  return -1;
}

// the qualitative determinant (outputs are 
// 0,1,-1 and 2 where 2 means unsigned)

int qualdetsubmat(int **imat, int n, int m, int *i1, int *i2, int dim){
  int i, j, r;
  int tot=0;
  int i2a[dim-1];
  /* assume this check has already been done. */
  /* if(minorhas0rc(imat, n, m, i1, i2, dim)) */
  /*   return 0; */
  if(dim==1)
    return sgnf(imat[i1[0]][i2[0]]);
  for (i=0;i<dim;i++){
    if(imat[i1[0]][i2[i]]!=0){
      r=0;
      for(j=0;j<dim;j++){
	if(j!=i)
	  i2a[r++]=i2[j];
      }
      tot=qualp(tot,qualt(par(i),qualt(imat[i1[0]][i2[i]],qualdetsubmat(imat, n, m, i1+1, i2a, dim-1))));
    }

  }
  return tot;

}






// here vec must be one larger than r

long getpos1(int *vec, int omit, int n, int r){
  long tot=0;
  int i,j, k1=n-1, k2=r-2, k3=0;
  for(j=0;j<r;j++){
    if(j!=omit){
      for(i=0;i<vec[j]-k3;i++){
	tot+=comb(k1, k2);
	k1--;
      }
      k1--;k2--;k3=vec[j]+1;
    }
  }
  return tot;
}

int ***allminors(int **imat, int n, int m){
  int i, s;//largest minors
  int ***tmp;
  s=min(n,m);


  tmp=(int ***) malloc((size_t)((s)*sizeof(int**)));
  if (!tmp) {fprintf(stderr, "allocation failure.\n");return NULL;}

  tmp[0]=cpmat(imat, n, m);
  for (i=1;i<s;i++){
    fprintf(stderr, "s = %d, i=%d***\n", s, i);
    tmp[i]=detsk1(imat, n, m, i+1, tmp[i-1]);
  }
  return tmp;

}

void free_allminors(int ***t, int n, int m)
{
  int i, s;
  s=min(n,m);
  for (i=0;i<s;i++){
    free_imatrix(t[i], 0, comb(n, i), 0, comb(m, i));
  }
  free((FREE_ARG) (t));
}


// calculates the minor of n times m matrix imat 
// indexed by xcombs and ycombs (which are of sign k)
// assumes that there is a matrix of minors (dets) of
// submatrices of imat of size k-1

int minor1(int **imat, int *xcombs, int *ycombs, int n, int m, int k, int **dets){
  int i;
  long i1,i2;
  int tot=0;
  fprintf(stderr, "xcombs = ");
  for(i=0;i<k;i++){
    fprintf(stderr, "%d ", xcombs[i]);
  }
  fprintf(stderr, "ycombs = ");
  for(i=0;i<k;i++){
    fprintf(stderr, "%d ", ycombs[i]);
  }
  fprintf(stderr, "\n");

  i1=getpos1(xcombs, 0, n, k);
  fprintf(stderr, "i1 = %li\n", i1);
  for (i=0;i<k;i++){
    i2=getpos1(ycombs, i, m, k);
    fprintf(stderr, "i1, i2 = %li, %li, entry=%d, dets = %d\n", i1, i2, imat[xcombs[0]][ycombs[i]], dets[i1][i2]);
    tot+=par(i)*imat[xcombs[0]][ycombs[i]]*dets[i1][i2];
  }
  fprintf(stderr, "tot = %d\n", tot);
  return tot;

}

const symbol & get_symbol(const string & s)
{
  static map<string, symbol> directory;
  map<string, symbol>::iterator i = directory.find(s);
  if (i != directory.end())
    return i->second;
  else
    return directory.insert(make_pair(s, symbol(s))).first->second;
}

const possymbol & get_possymbol(const string & s)
{
  static map<string, possymbol> directory;
  map<string, possymbol>::iterator i = directory.find(s);
  if (i != directory.end())
    return i->second;
  else
    return directory.insert(make_pair(s, possymbol(s))).first->second;
}

/* generate both factors in the clever factorisation of J[2] */
/* The second factor is stored as a set of signed integers */
/* which can be converted to strings when needed */

int genS2(int **imat, int n, int m, int ***mat1, int ***mat2){
  int **xcombs;
  long r1, r2,cnk,cmk;
  int **tmp;

  tmp = mattoind(imat, n, m);

  cnk=n*(n-1)/2;cmk=n*m;
  (*mat1)=imatrix(0, cnk-1, 0, cmk-1);
  (*mat2)=imatrix(0, cnk-1, 0, cmk-1);
  xcombs=allcombsgen(n,2);
  for(r1=0;r1<cnk;r1++){
    for(r2=0;r2<cmk;r2++){
      /* fprintf(stderr, "%d, %d, %d\n", r2%n,xcombs[r1][0],xcombs[r1][1]); */
      if(r2%n==xcombs[r1][0]){
	(*mat1)[r1][r2] = -imat[xcombs[r1][1]][r2/n];
	(*mat2)[r1][r2] = -tmp[xcombs[r1][1]][r2/n];
      }
      else if (r2%n==xcombs[r1][1]){
	(*mat1)[r1][r2] = imat[xcombs[r1][0]][r2/n];
	(*mat2)[r1][r2] = tmp[xcombs[r1][0]][r2/n];
      }
      else{
	(*mat1)[r1][r2] = 0;
	(*mat2)[r1][r2] = 0;
      }
    }
  }
  free_imatrix(xcombs, 0, cnk-1, 0, 1);
  free_imatrix(tmp, 0, n-1, 0, m-1);
  return 1;
}

/* generate both factors in the clever factorisation of J[2] */
/* The second factor is stored as a set of signed integers */
/* which can be converted to strings when needed */
/* Here it is assumed that both factors of J are given */
/* and may be different */

int genS2from2(int **imatS, int **imatV, int n, int m, int ***mat1, int ***mat2){
  int **xcombs;
  long r1, r2,cnk,cmk;

  cnk=n*(n-1)/2;cmk=n*m;
  (*mat1)=imatrix(0, cnk-1, 0, cmk-1);
  (*mat2)=imatrix(0, cnk-1, 0, cmk-1);
  xcombs=allcombsgen(n,2);
  for(r1=0;r1<cnk;r1++){
    for(r2=0;r2<cmk;r2++){
      /* fprintf(stderr, "%d, %d, %d\n", r2%n,xcombs[r1][0],xcombs[r1][1]); */
      if(r2%n==xcombs[r1][0]){
	(*mat1)[r1][r2] = -imatS[xcombs[r1][1]][r2/n];
	(*mat2)[r1][r2] = -imatV[xcombs[r1][1]][r2/n];
      }
      else if (r2%n==xcombs[r1][1]){
	(*mat1)[r1][r2] = imatS[xcombs[r1][0]][r2/n];
	(*mat2)[r1][r2] = imatV[xcombs[r1][0]][r2/n];
      }
      else{
	(*mat1)[r1][r2] = 0;
	(*mat2)[r1][r2] = 0;
      }
    }
  }
  free_imatrix(xcombs, 0, cnk-1, 0, 1);
  return 1;
}

int genS2from2symb(int **imatS, matrix imatV, int n, int m, int ***mat1, matrix *mat2){
  int **xcombs;
  long r1, r2,cnk,cmk;

  cnk=n*(n-1)/2;cmk=n*m;
  (*mat1)=imatrix(0, cnk-1, 0, cmk-1);
  xcombs=allcombsgen(n,2);
  for(r1=0;r1<cnk;r1++){
    for(r2=0;r2<cmk;r2++){
      /* fprintf(stderr, "%d, %d, %d\n", r2%n,xcombs[r1][0],xcombs[r1][1]); */
      if(r2%n==xcombs[r1][0]){
	(*mat1)[r1][r2] = -imatS[xcombs[r1][1]][r2/n];
	(*mat2)(r1,r2) = -imatV(xcombs[r1][1],r2/n);
      }
      else if (r2%n==xcombs[r1][1]){
	(*mat1)[r1][r2] = imatS[xcombs[r1][0]][r2/n];
	(*mat2)(r1,r2) = imatV(xcombs[r1][0],r2/n);
      }
      else{
	(*mat1)[r1][r2] = 0;
	(*mat2)(r1,r2) = 0;
      }
    }
  }
  free_imatrix(xcombs, 0, cnk-1, 0, 1);
  return 1;
}

int genS2symb(matrix mat_in, int n, int m, matrix *mat2){
  int **xcombs;
  long r1, r2,cnk,cmk;

  cnk=n*(n-1)/2;cmk=n*m;
  xcombs=allcombsgen(n,2);
  for(r1=0;r1<cnk;r1++){
    for(r2=0;r2<cmk;r2++){
      if(r2%n==xcombs[r1][0])
	(*mat2)(r1,r2) = -mat_in(xcombs[r1][1],r2/n);
      else if (r2%n==xcombs[r1][1])
	(*mat2)(r1,r2) = mat_in(xcombs[r1][0],r2/n);
      else
	(*mat2)(r1,r2) = 0;
    }
  }
  free_imatrix(xcombs, 0, cnk-1, 0, 1);
  return 1;
}



int **detsk1(int **imat, int n, int m, int k, int **dets){
  int **xcombs;
  int **ycombs;
  int **tmp;
  long r1, r2,cnk,cmk;


  tmp=imatrix(0, comb(n, k)-1, 0, comb(m, k)-1);
  xcombs=allcombsgen(n,k);
  ycombs=allcombsgen(m,k);
  cnk=comb(n,k);cmk=comb(m,k);
  for(r1=0;r1<cnk;r1++){
    for(r2=0;r2<cmk;r2++){
      tmp[r1][r2]=minor1(imat, xcombs[r1], ycombs[r2], n, m, k, dets);
    }
  }
  free_imatrix(xcombs, 0, cnk-1, 0, k-1);
  free_imatrix(ycombs, 0, cmk-1, 0, k-1);
  return tmp;

}

bool **bmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a bool matrix with subscript range m[nrl..nrh][ncl..nch] */
{
  long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
  bool **m;

  /* allocate pointers to rows */
  m=(bool **) malloc((size_t)((nrow+NR_END)*sizeof(bool*)));
  if (!m) {fprintf(stderr, "allocation failure 1 in matrix()"); exit(0);}
  m += NR_END;
  m -= nrl;


  /* allocate rows and set pointers to them */
  m[nrl]=(bool *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(bool)));
  if (!m[nrl]) {fprintf(stderr, "allocation failure 1 in matrix()"); exit(0);}
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

  /* return pointer to array of pointers to rows */
  return m;
}

int **imatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
  long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
  int **m;

  /* allocate pointers to rows */
  m=(int **) malloc((size_t)((nrow+NR_END)*sizeof(int*)));
  if (!m) {fprintf(stderr, "allocation failure 1 in matrix()"); exit(0);}
  m += NR_END;
  m -= nrl;


  /* allocate rows and set pointers to them */
  m[nrl]=(int *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(int)));
  if (!m[nrl]) {fprintf(stderr, "allocation failure 1 in matrix()"); exit(0);}
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

  /* return pointer to array of pointers to rows */
  return m;
}

void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch)
/* free an int matrix allocated by imatrix() */
{
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
}

void free_bmatrix(bool **m, long nrl, long nrh, long ncl, long nch)
/* free an int matrix allocated by imatrix() */
{
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
}

ex **exmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a symb matrix with subscript range m[nrl..nrh][ncl..nch] */
{
  long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
  ex **m;

  /* allocate pointers to rows */
  m=(ex **) malloc((size_t)((nrow+NR_END)*sizeof(ex*)));
  if (!m) {fprintf(stderr, "allocation failure 1 in matrix()"); exit(0);}
  m += NR_END;
  m -= nrl;


  /* allocate rows and set pointers to them */
  m[nrl]=(ex *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(ex)));
  if (!m[nrl]) {fprintf(stderr, "allocation failure 1 in matrix()"); exit(0);}
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

  /* return pointer to array of pointers to rows */
  return m;
}

void free_exmatrix(ex **m, long nrl, long nrh, long ncl, long nch)
/* free an ex matrix allocated by exmatrix() */
{
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
}

char *readfileintostr(const char fname[]){
  FILE *fd;
  long i=0;
  int c=0;
  char *str;

  // first parse - find length

  fd = fopen(fname, "r");
  if(!fd){
    fprintf(stderr, "ERROR in readfileintostr: \"%s\" could not be opened for reading.\n", fname);
    str=strdup("");
    return str;
  }

  while((c=getc(fd)) != EOF)
    i++;
 
  i++;
  str = (char*) malloc(sizeof(char) * (i));
  // second parse - get number of non-comment lines

  rewind(fd);
  i=0;
  while((c=getc(fd)) != EOF)
    str[i++]=c;
 
  str[i++]=0;
 
  fclose(fd);
  return str;
}

int getline0(FILE *fp, char s[], int lim)
{
  int c=0, i;

  for(i=0; i<lim-1 && (c=getc(fp))!=EOF && c!=13 && c!='\n';++i)
    s[i] = c;
  s[i] = '\0';
  if(i==lim-1 && c!='\n' && c!=EOF)
    fprintf(stderr, "WARNING in getline0: line length exceeded. %s\n", s);
  if (c == '\n' || c == 13) { // don't do anything
    return i+1;
  }
  else
    return i;
  // we return the length it would have if it had a newline character at the end
}


char *getlinefromstr(long *pos, char *block)
{
  int c=0, i=0;
  long len;
  char *str, *p;
  if(*pos > (int)strlen(block) || (*pos < 0)){
    str=strdup("");
    return str;
  }
  p=strchr(block+(*pos), '\n');
  if(p){
    len=strlen(block+(*pos)) - strlen(p);
    str = (char*) malloc(sizeof(char) * (len + 2));
    i=0;
    while((c=block[(*pos)++])!='\n')
      str[i++]=c;
    str[i++] = c;
    str[i] = 0;
  }
  else{
    str=strdup(block+(*pos));
    (*pos)=strlen(block);
  }
  return str;
}

char *strchop2(char *t, int n, int n1){
  // memory allocating version of strchop1
  char *s;
  int i=0;
  if((n<= (int)strlen(t)) && n1>0){
    s = (char*) malloc(sizeof(char) * (n1+1));
    while (t[n+i] && (i< n1)){s[i] = t[n+i];i++;}
    s[i] = '\0';
  }
  else
    s = strdup("");
  
  return s;
}

int getnumints(char *s){
  int i, k;
  i=0, k=0;
  while(s[k]){
    while(s[k] && !isdigit(s[k])&& (s[k]!='-')&& (s[k]!='+')){k++;}
    if(isdigit(s[k])|| (s[k]=='-')|| (s[k]=='+')){
      i++;
      while(isdigit(s[k])|| (s[k]=='-')|| (s[k]=='+')){k++;}
    }
  }
  return i;
}

int getnumints1(char *s){
  int i, k;
  i=0, k=0;
  while(s[k]){
    while(s[k] && !isdigit(s[k])&& (s[k]!='-')&& (s[k]!='+')){k++;}
    if(isdigit(s[k])|| (s[k]=='-')|| (s[k]=='+')){
      i++;
      while((s[k]=='-')||(s[k]=='+')){k++;}
      while(isdigit(s[k])){k++;}
    }
  }
  return i;
}

// number of segments in a line, separated by space
int getnumsegs(char *s){
  int i=0,k=0;
  while(s[k]){
    while(s[k] && isspace(s[k])){k++;} // skip initial space
    if(s[k] && !isspace(s[k])){
      i++;
      while(s[k] && !isspace(s[k])){k++;}//skip nonspace
    }
  }
  return i;
}

char *getnthwd(char *s, int n){
  int i, j, k=0;
  char *v=NULL;
  while(s[k]){
    while(s[k] && isspace(s[k])){k++;} // skip space
    for(i=0;i<n-1;i++){
      while(s[k] && !isspace(s[k])){k++;} // skip word
      while(s[k] && isspace(s[k])){k++;} // skip nonword
    }
    j=0;
    while(!isspace(s[k])){j++;k++;} // get the word
    v = strchop2(s, k-j, j);
    return v;
  }
  if(!v)
    v=strdup("");
  return v;
}

char *getnthint(char *s, int n){
  int i, j, k;
  i=0, k=0;
  char *v=NULL;
  while(s[k] != 0){
    while(s[k] && !isdigit(s[k]) && (s[k]!='-')&& (s[k]!='+')){k++;} // skip nonwords
    for(i=0;i<n-1;i++){
      while(isdigit(s[k]) || (s[k]=='-')|| (s[k]=='+')){k++;} // skip word
      while(s[k] && !isdigit(s[k]) && (s[k]!='-')&& (s[k]!='+')){k++;} // skip nonwords
    }
    j=0;
    while((s[k]=='-')||(s[k]=='+')||isdigit(s[k])){j++;k++;} // get the word
    v = strchop2(s, k-j, j);
    return v;
  }
  if(!v)
    v=strdup("");
  return v;
}

char *getnthint1(char *s, int n){
  int i=0, j, k=0;
  char *v=NULL;
  while(s[k]){
    while(s[k] && !isdigit(s[k]) && (s[k]!='-')&& (s[k]!='+')){k++;} // skip nonwords
    for(i=0;i<n-1;i++){
      while((s[k]=='-')||(s[k]=='+')){j++;k++;} //skip sign
      while(isdigit(s[k])){j++;k++;} // skip int
      while(s[k] && !isdigit(s[k]) && (s[k]!='-')&& (s[k]!='+')){k++;} // skip nonwords
    }
    j=0;
    while((s[k]=='-')||(s[k]=='+')){j++;k++;} // get the word
    while(isdigit(s[k])){j++;k++;} // get the word
    v = strchop2(s, k-j, j);
    return v;
  }
  if(!v)
    v=strdup("");
  return v;
}

// get the first integer in a string. Return pointer to just after

int getint(char **s){
  char *v;
  int j, k=0, r;

  while((*s)[k] && !isdigit((*s)[k]) && ((*s)[k]!='-')&& ((*s)[k]!='+')){k++;} // skip nonwords
  j=0;
  while(((*s)[k]=='-')||((*s)[k]=='+')){j++;k++;} // get the word
  while(isdigit((*s)[k])){j++;k++;} // get the word
  v = strchop2((*s), k-j, j);
  r=atoi(v);
  //cout << "line = \"" << *s << "\"" << endl;
  //cout << "v = " << v << ", k = " << k << endl;
  free(v);
  (*s)=(*s)+k;
  //cout << "newline = \"" << *s << "\"" << endl;
  return r;
}

unsigned long factorial(int x){
  int i;
  unsigned long factx = 1;
  for(i=1; i<=x ; i++ )
    factx *= i;
  return factx;
}

//not pretty, but works to reasonably large values

long comb(int n, int k){
  int i,j,ind1=2;
  long combx=1;
  int inds[n-k];
  for(i=0;i<n-k;i++){
    inds[i]=1;
  }
  j=ind1;
  for(i=k+1; i<=n ; i++ ){//fprintf(stderr, "%d*%d\t", combx, i);
    combx *= i;
    j=ind1;
    while(j<=n-k){ // try to remove factors
      if(combx%j==0 && inds[j-1]==1 ){
	combx/=j; inds[j-1]=0; if(j==ind1){ind1++;}
      }
      j++;
    }
  }
  return combx;
}



int par(int k){
  if (k%2==0)
    return 1;
  return -1;
}

//
// is the integer "i" in the list lst?
//

int isinlist(int i, int ilst[], int tot){
  int j;
  for(j=0;j<tot;j++){
    if(i==ilst[j])
      return 1;
  }
  return 0;
}

bool isinllist(long i, long ilst[], int tot){
  int j;
  for(j=0;j<tot;j++){
    if(i==ilst[j])
      return 1;
  }
  return 0;
}

int isonlyspace(char *s){
  int c, k=0;
  while((c=s[k++]) != '\0'){
    if(!(isspace(c)))
      return 0;
  }
  return 1;
}

int iscomline(char s[]){ // a comment line or separator line or empty
  int i=0;
  while((isspace((int) s[i]))){i++;}
  if (!s[i] || ((s[i] == '/') && (s[i+1] == '/')) || ((s[i] == '/') && (s[i+1] == '*')) || (s[i] == '#')){
    return 1;
  }
  return 0;
}

void printmat(int **imat, int n, int m){
  int i,j;
  for(i=0;i<n;i++){
    for(j=0;j<m;j++)
      fprintf(stderr, "%2d ", imat[i][j]);
    fprintf(stderr, "\n");
  }
  fprintf(stderr, "\n");
  return;
}

void printtens(int ***tens, int n, int m){
  int i,j,k;
  for(i=0;i<n;i++){
    for(j=0;j<m;j++){
      if(!tens[i][j][0])
	fprintf(stderr, "%2d ", 0);
      else{
	fprintf(stderr, " [");
	for(k=1;k<=tens[i][j][0];k++)
	  fprintf(stderr, "%2d ", tens[i][j][k]);
	fprintf(stderr, "] ");
      }
    }
    fprintf(stderr, "\n");
  }
  fprintf(stderr, "\n");
  return;
}

void printbmat(bool **imat, int n, int m){
  int i,j;
  for(i=0;i<n;i++){
    for(j=0;j<m;j++)
      fprintf(stderr, "%2d ", imat[i][j]);
    fprintf(stderr, "\n");
  }
  fprintf(stderr, "\n");
  return;
}

void printvec(int *ivec, int n){
  int i;
  for(i=0;i<n;i++)
    fprintf(stderr, "%d ", ivec[i]);
  fprintf(stderr, "\n");
  return;
}

void printbvec(bool *ivec, int n){
  int i;
  for(i=0;i<n;i++)
    fprintf(stderr, "%d ", ivec[i]);
  fprintf(stderr, "\n");
  return;
}

void printlvec(long *ivec, int n){
  int i;
  for(i=0;i<n;i++)
    fprintf(stderr, "%ld ", ivec[i]);
  fprintf(stderr, "\n");
  return;
}

void printvec1(int *ivec, int n){
  int i;
  for(i=0;i<n;i++)
    printf("%d ", ivec[i]);
  printf("\n");
  return;
}

void printexmat(matrix imat, int n, int m){
  int i,j;
  cout << "symbolic matrix:\n";
  for(i=0;i<n;i++){
    for(j=0;j<m;j++){
      cout << imat(i, j) << " ";
    }
    cout << endl;
  }
  cout << endl;
  return;
}

void printsubmat(int **imat, int *vec1, int *vec2, int k1, int k2){
  int i,j;
  for(i=0;i<k1;i++){
    for(j=0;j<k2;j++){
      fprintf(stderr, "%2d ",imat[vec1[i]][vec2[j]]);
    }
    fprintf(stderr, "\n");
  }
  fprintf(stderr, "\n");
  return;
}

void printexsubmat(matrix imat, int *vec1, int *vec2, int n, int m){
  int i,j;
  cout << "symbolic matrix:\n";
  for(i=0;i<n;i++){
    for(j=0;j<m;j++)
      cout << imat(vec1[i], vec2[j]) << "\t";
    cout << endl;
  }
  cout << endl;
  return;
}

void printindsubmat(int **imat, int *vec1, int *vec2, int k1, int k2){
  int i,j;
  for(i=0;i<k1;i++){
    for(j=0;j<k2;j++){
      if(imat[vec1[i]][vec2[j]]==0)
	fprintf(stderr, "0  ");
      else if(imat[vec1[i]][vec2[j]]<0)
	fprintf(stderr, "-v%d ",-imat[vec1[i]][vec2[j]]);
      else
	fprintf(stderr, "v%d ",imat[vec1[i]][vec2[j]]);
    }
    fprintf(stderr, "\n");
  }
  fprintf(stderr, "\n");
  return;
}

void printmaximaindsubmat(FILE *fd, int **imat, int *vec1, int *vec2, int k1, int k2){
  int i,j;
  fprintf(fd, "MM:matrix(");
  for(i=0;i<k1;i++){
    fprintf(fd, "[");
    for(j=0;j<k2;j++){
      if(imat[vec1[i]][vec2[j]]==0)
	fprintf(fd, "0");
      else if(imat[vec1[i]][vec2[j]]<0)
	fprintf(fd, "-v%d",-imat[vec1[i]][vec2[j]]);
      else
	fprintf(fd, "v%d",imat[vec1[i]][vec2[j]]);
      if(j<(k2-1))
	fprintf(fd, ",");
    }
    if(i<k1-1)
      fprintf(fd, "],");
    else
      fprintf(fd, "]");
  }
  fprintf(fd, ");\n");
  fprintf(fd, "expand(determinant(MM));\n");
  return;
}

// Print a submatrix as an (undirected) graph in the "dot" language

void printdotindsubmat(FILE *fd, int **imat, int *vec1, int *vec2, int k1, int k2, long lab){
  int i,j;
  fprintf(fd, "graph mat%ld {\n", lab);
  for(i=0;i<k1;i++){
    for(j=0;j<k2;j++){
      if(imat[vec1[i]][vec2[j]]==0){}
      else if(imat[vec1[i]][vec2[j]]<0)
	fprintf(fd, "\tR%d -- C%d [label=v%d] [color=red]\n",vec1[i],vec2[j],-imat[vec1[i]][vec2[j]]);
      else
	fprintf(fd, "\tR%d -- C%d [label=v%d]\n",vec1[i],vec2[j],imat[vec1[i]][vec2[j]]);
    }
  }
  fprintf(fd, "}\n\n");
  return;
}

// print a pair of submatrices as a directed graph in the dot language

void printdotpairsubmat(FILE *fd, int **imat, int **imat1, int *vec1, int *vec2, int k1, int k2, long lab){
  int i,j;
  fprintf(fd, "digraph mat%ld {\n", lab);
  for(i=0;i<k1;i++){
    for(j=0;j<k2;j++){
      if(imat[vec1[i]][vec2[j]]==0 && imat1[vec1[i]][vec2[j]]==0){}
      else if(imat[vec1[i]][vec2[j]]<0 && imat1[vec1[i]][vec2[j]]<0)
	fprintf(fd, "\tR%d -> C%d [label=v%d] [color=red] [dir=none]\n",vec1[i],vec2[j],-imat1[vec1[i]][vec2[j]]);
      else if(imat[vec1[i]][vec2[j]]>0 && imat1[vec1[i]][vec2[j]]>0)
	fprintf(fd, "\tR%d -> C%d [label=v%d] [dir=none]\n",vec1[i],vec2[j],imat1[vec1[i]][vec2[j]]);
      else{
	if(imat[vec1[i]][vec2[j]]<0)
	  fprintf(fd, "\tC%d -> R%d [label=v%d] [color=red]\n",vec1[i],vec2[j],-imat1[vec1[i]][vec2[j]]);
	if(imat1[vec1[i]][vec2[j]]<0)
	  fprintf(fd, "\tR%d -> C%d [label=v%d] [color=red]\n",vec1[i],vec2[j],-imat1[vec1[i]][vec2[j]]);
	if(imat[vec1[i]][vec2[j]]>0)
	  fprintf(fd, "\tC%d -> R%d [label=v%d]\n",vec1[i],vec2[j],imat1[vec1[i]][vec2[j]]);
	if(imat1[vec1[i]][vec2[j]]>0)
	  fprintf(fd, "\tR%d -> C%d [label=v%d]\n",vec1[i],vec2[j],imat1[vec1[i]][vec2[j]]);
      }
    }
  }
  fprintf(fd, "}\n\n");
  return;
}

void printmaximamat(int **imat, int n, int m){
  int i,j;
  fprintf(stderr, "SS:matrix(");
  for(i=0;i<n;i++){
    fprintf(stderr, "[");
    for(j=0;j<m;j++){
      fprintf(stderr, "%d",imat[i][j]);
      if(j<(m-1))
	fprintf(stderr, ",");
    }
   if(i<n-1)
      fprintf(stderr, "],");
    else
      fprintf(stderr, "]");
  }
  fprintf(stderr, ");\n");
  return;
}

void printmaximaindmat(int **imat, int k1, int k2){
  int i,j;
  fprintf(stderr, "MM:matrix(");
  for(i=0;i<k1;i++){
    fprintf(stderr, "[");
    for(j=0;j<k2;j++){
      if(imat[i][j]==0)
	fprintf(stderr, "0");
      else if(imat[i][j]<0)
	fprintf(stderr, "-v%d",-imat[i][j]);
      else
	fprintf(stderr, "v%d",imat[i][j]);
      if(j<(k2-1))
	fprintf(stderr, ",");
    }
    if(i<k1-1)
      fprintf(stderr, "],");
    else
      fprintf(stderr, "]");
  }
  fprintf(stderr, ");\n");
  return;
}




void printmat1(char **cmat, int **imat, int n, int m){
  int i,j;
  for(i=0;i<n;i++){
    fprintf(stderr, "%s\t", cmat[i]);
    for(j=0;j<m;j++){
      fprintf(stderr, "%2d  ", imat[i][j]);
    }
    fprintf(stderr, "\n");
  }
  fprintf(stderr, "\n");
  return;
}

void printmatpat(char **cmat, int **imat, int n, int m){
  int i,j;
  for(i=0;i<n;i++){
    fprintf(stderr, "%s\t", cmat[i]);
    for(j=0;j<m;j++){
      if(imat[i][j]==2)
	fprintf(stderr, " x  ");
      else
	fprintf(stderr, "%2d  ", imat[i][j]);
    }
    fprintf(stderr, "\n");
  }
  fprintf(stderr, "\n");
  return;
}

void chop(char *str, char **str1, char **str2, const char str3[]){
  char *line;
  long pos=0;
  int pass=0;
  (*str1)=strdup(str);
  (*str2)=strdup(str);
  line = getlinefromstr(&pos, str);
  while(strcmp(line, "")!=0 && !strstr(line, str3)){
    if(pass==0){strcpy((*str1), line);pass=1;}else{strcat((*str1), line);}
    free(line);line = getlinefromstr(&pos, str);
  }
  free(line);
  line = getlinefromstr(&pos, str);
  pass=0;
  while(strcmp(line, "")!=0){
    if(pass==0){strcpy((*str2), line);pass=1;}else{strcat((*str2), line);}
    free(line);line = getlinefromstr(&pos, str);
  }
  free(line);

}

/* read a pair of matrices from a string */
/* matrices are separated by ***** */

int readmatpairfromstr(char *str, int *nlen, int *mlen, int ***mat1, int ***mat2){
  char *str1, *str2;
  int nl1, ml1, nl2, ml2;

  chop(str, &str1, &str2, "*****");

  (*mat1)=readmatrixfromstr(str1, &nl1, &ml1);

  if(!(*mat1)){
    free(str1);free(str2);
    return 0;
  }
  (*mat2)=readmatrixfromstr(str2, &nl2, &ml2);
  if(!(*mat2) || (nl1!=nl2) || (ml1!=ml2)){
    fprintf(stderr, "ERROR: matrices could not be read.\n");
    free_imatrix((*mat1), 0, nl1, 0, ml1);
    (*mat1)=NULL;
    if((*mat2)){
      free_imatrix((*mat2), 0, nl2, 0, ml2);
      (*mat2)=NULL;
    }

    free(str1);free(str2);
    return 0;
  }
  (*nlen)=nl1;(*mlen)=ml1;
  free(str1);free(str2);
  return 1;

}

/* read a matrix and its rank from a string */
/* matrices are separated by ***** */

int readmatrankfromstr(char *str, int *nlen, int *mlen, int ***mat1, int *rk){
  char *str1, *str2;
  int nl1, ml1, nl2=0, ml2=0;
  int **mat2;

  chop(str, &str1, &str2, "*****");

  (*mat1)=readmatrixfromstr(str1, &nl1, &ml1);

  if(!(*mat1)){
    free(str1);free(str2);
    fprintf(stderr, "ERROR: matrix could not be read.\n");
    return 0;
  }
  mat2=readmatrixfromstr(str2, &nl2, &ml2);
  if(!mat2 || nl2!=1 || ml2!=1){
    fprintf(stderr, "ERROR: rank could not be read. (It must be an integer)\n");
    free_imatrix((*mat1), 0, nl1-1, 0, ml1-1);
    (*mat1)=NULL;
    if(mat2){
      free_imatrix(mat2, 0, nl2-1, 0, ml2-1);
      mat2=NULL;
    }
    free(str1);free(str2);
    return 0;
  }
  (*nlen)=nl1;(*mlen)=ml1;
  (*rk) = mat2[0][0];
  /*  fprintf(stderr, "rk = %d\n", *rk); */
  free_imatrix(mat2, 0, nl2-1, 0, ml2-1);
  free(str1);free(str2);
  return 1;

}

int readmatpairrankfromstr(char *str, int *nlen, int *mlen, int ***mat1, int ***mat2, int *rk){
  char *str1, *strtmp, *str2, *str3;
  int nl1, ml1, nl2, ml2, nl3, ml3;
  int **mat3;

  chop(str, &str1, &strtmp, "*****");
  fprintf(stderr, "str1 = %s\n", str1);
  fprintf(stderr, "strtmp = %s\n", strtmp);
  chop(strtmp, &str2, &str3, "*****");
  fprintf(stderr, "str2 = %s\n", str2);
  fprintf(stderr, "str3 = %s\n", str3);

  (*mat1)=readmatrixfromstr(str1, &nl1, &ml1);

  if(!(*mat1)){
    free(str1);free(str2);free(strtmp);free(str3);
    return 0;
  }
  (*mat2)=readmatrixfromstr(str2, &nl2, &ml2);
  if(!(*mat2) || (nl1!=nl2) || (ml1!=ml2)){
    fprintf(stderr, "ERROR: matrices could not be read.\n");
    free_imatrix((*mat1), 0, nl1, 0, ml1);
    (*mat1)=NULL;
    if((*mat2)){
      free_imatrix((*mat2), 0, nl2, 0, ml2);
      (*mat2)=NULL;
    }

    free(str1);free(str2);free(strtmp);free(str3);
    return 0;
  }
  mat3=readmatrixfromstr(str3, &nl3, &ml3);
  if(!mat3 || nl3!=1 || ml3!=1){
    //    fprintf(stderr, "ttnl3 = %d, ml3 = %d\n", nl3, ml3);
    fprintf(stderr, "ERROR: rank could not be read. (It must be an integer)\n");
    free_imatrix((*mat1), 0, nl1-1, 0, ml1-1);
    free_imatrix((*mat2), 0, nl2-1, 0, ml2-1);
    (*mat1)=NULL;(*mat2)=NULL;
    if(mat3){
      free_imatrix(mat3, 0, nl3-1, 0, ml3-1);
      mat3=NULL;
    }
    free(str1);free(str2);free(strtmp);free(str3);
    return 0;
  }
  (*nlen)=nl1;(*mlen)=ml1;
  (*rk) = mat3[0][0];
  /*  fprintf(stderr, "rk = %d\n", *rk); */
  free_imatrix(mat3, 0, nl3-1, 0, ml3-1);


  (*nlen)=nl1;(*mlen)=ml1;
  free(str1);free(str2);free(strtmp);free(str3);
  return 1;

}

int ***readiimatrixfromstr(char *str, int *nlen, int *mlen){
  int numlines=0, numtmp, i,j,nt;
  char *line;
  long pos=0;
  int ***imat;
  char *wd,*wd0;
  line = getlinefromstr(&pos, str);
  (*mlen)=0;(*nlen)=0;
  while(strcmp(line, "")!=0){
    if(!iscomline(line)){
      numlines++;
      numtmp=getnumsegs(line);
      if(numlines!=1 && numtmp!=(*mlen)){fprintf(stderr, "WARNING: matrix does not seem well defined.\n");}
      if(numtmp>(*mlen)){(*mlen)=numtmp;}
    }
    free(line);
    line = getlinefromstr(&pos, str);
  }

  free(line);
  (*nlen)=numlines;

  // no matrix found

  if((*nlen)==0 || (*mlen)==0){
    fprintf(stderr, "No matrix found in string %s.\n", str);
    return NULL;
  }

  imat=(int***) malloc(sizeof(int**)*(*nlen));
  for(i=0;i<(*nlen);i++)
    imat[i]=(int**) malloc(sizeof(int*)*(*mlen));

  pos=0;numlines=0;
  line = getlinefromstr(&pos, str);
  while(strcmp(line, "")!=0){
    if(!iscomline(line)){
      numlines++;
      for(i=1;i<(*mlen)+1;i++){
	wd0=getnthwd(line, i);
	nt=getnumints1(wd0);
	imat[numlines-1][i-1]=(int*) malloc(sizeof(int)*(nt+1));
	imat[numlines-1][i-1][0]=nt;
	for(j=1;j<=nt;j++){
	  wd=getnthint1(wd0, j);
	  imat[numlines-1][i-1][j]=atoi(wd);
	  free(wd);
	}
	free(wd0);
      }
    }
    free(line);
    line = getlinefromstr(&pos, str);
  }
  free(line);
  return imat;

}

void freeiim(int ***mat, int n, int m){
  int i,j;
  for(i=0;i<n;i++){
    for(j=0;j<m;j++){
      free((char *)mat[i][j]);
    }
    free((char *)mat[i]);
  }
  free((char *)mat);
  return;
}


int readiimatpairrankfromstr(char *str, int *nlen, int *mlen, int ***mat1, int ****mat2, int *rk){
  char *str1, *strtmp, *str2, *str3;
  int nl1, ml1, nl2, ml2, nl3, ml3;
  int **mat3;

  chop(str, &str1, &strtmp, "*****");
  fprintf(stderr, "str1 = %s\n", str1);
  fprintf(stderr, "strtmp = %s\n", strtmp);
  chop(strtmp, &str2, &str3, "*****");
  fprintf(stderr, "str2 = %s\n", str2);
  fprintf(stderr, "str3 = %s\n", str3);

  (*mat1)=readmatrixfromstr(str1, &nl1, &ml1);

  if(!(*mat1)){
    free(str1);free(str2);free(strtmp);free(str3);
    return 0;
  }
  (*mat2)=readiimatrixfromstr(str2, &nl2, &ml2);
  if(!(*mat2) || (nl1!=nl2) || (ml1!=ml2)){
    fprintf(stderr, "ERROR: matrices could not be read.\n");
    free_imatrix((*mat1), 0, nl1, 0, ml1);
    (*mat1)=NULL;
    if((*mat2)){
      freeiim((*mat2), nl2, ml2);
      (*mat2)=NULL;
    }

    free(str1);free(str2);free(strtmp);free(str3);
    return 0;
  }
  mat3=readmatrixfromstr(str3, &nl3, &ml3);
  if(!mat3 || nl3!=1 || ml3!=1){
    //    fprintf(stderr, "ttnl3 = %d, ml3 = %d\n", nl3, ml3);
    fprintf(stderr, "ERROR: rank could not be read. (It must be an integer)\n");
    free_imatrix((*mat1), 0, nl1-1, 0, ml1-1);
    freeiim((*mat2), nl2, ml2);
    (*mat1)=NULL;(*mat2)=NULL;
    if(mat3){
      free_imatrix(mat3, 0, nl3-1, 0, ml3-1);
      mat3=NULL;
    }
    free(str1);free(str2);free(strtmp);free(str3);
    return 0;
  }
  (*nlen)=nl1;(*mlen)=ml1;
  (*rk) = mat3[0][0];
  /*  fprintf(stderr, "rk = %d\n", *rk); */
  free_imatrix(mat3, 0, nl3-1, 0, ml3-1);


  (*nlen)=nl1;(*mlen)=ml1;
  free(str1);free(str2);free(strtmp);free(str3);
  return 1;

}






int **readmatrixfromstr(char *str, int *nlen, int *mlen){
  int numlines=0, numtmp, i;
  char *line;
  long pos=0;
  int **imat;
  char *wd;
  line = getlinefromstr(&pos, str);
  (*mlen)=0;(*nlen)=0;
  while(strcmp(line, "")!=0){
    if(!iscomline(line)){
      numlines++;

      numtmp=getnumints(line);
      if(numlines!=1 && numtmp!=(*mlen)){fprintf(stderr, "WARNING: matrix does not seem well defined.\n");}
      if(numtmp>(*mlen)){(*mlen)=numtmp;}
    }
    free(line);
    line = getlinefromstr(&pos, str);
  }

  free(line);
  (*nlen)=numlines;

  // no matrix found

  if((*nlen)==0 || (*mlen)==0){
    fprintf(stderr, "No matrix found in string %s.\n", str);
    return NULL;
  }

  // otherwise create the matrix

  imat=imatrix(0, numlines-1, 0, (*mlen)-1);
  
  pos=0;numlines=0;
  line = getlinefromstr(&pos, str);
  while(strcmp(line, "")!=0){
    if(!iscomline(line)){
      numlines++;
      for(i=1;i<(*mlen)+1;i++){
	wd=getnthint(line, i);
	imat[numlines-1][i-1]=atoi(wd);
	free(wd);
      }
    }
    free(line);
    line = getlinefromstr(&pos, str);
  }
  free(line);
  return imat;

}

int isreac(char *str){
  if(strstr(str, "<-->") || strstr(str, "<==>") || strstr(str, "<->")|| strstr(str, "<>") || strstr(str, "<=>")|| strstr(str, "-->")|| strstr(str, "==>")|| strstr(str, "->")|| strstr(str, ">")|| strstr(str, "=>")|| strstr(str, "<--")|| strstr(str, "<==")|| strstr(str, "<-")|| strstr(str, "<")|| strstr(str, "<="))
    return 1;
  return 0;
}

int getallreacs(char *str, int ***imat1, int ***imat2, int ***imat3, int ***imat4, char ***chems, int *n, int *m, int *cols3, int *allrev, int *allgood){
  // imat3 is the irreversible stoichiometric matrix
  // imat4 is the corresponding matrix of powers (a mass action like system)
  long pos=0;
  int pos1;
  char *line;
  int totchems=0;
  int i, ind;
  int irrcols=0, numcols=0;
  int lstoic, rstoic;

  char **leftchems, **rightchems;
  int *leftstoics,*rightstoics;
  int numlines=0, numleft, numright, rev=0;
  (*allrev)=1;
  (*allgood)=1;
  (*chems)=NULL;

  fprintf(stderr, "%s\n", str);
  line = getlinefromstr(&pos, str);
  while(strcmp(line, "")!=0){
    if(!iscomline(line) && isreac(line)){
      if(!getreac(line, &leftchems, &leftstoics, &numleft, &rightchems, &rightstoics, &numright, &rev)){
	freearraydat(leftchems, numleft);
	freearraydat(rightchems, numright);
	freearraydat((*chems), totchems);
	free((FREE_ARG) (leftstoics));
	free((FREE_ARG) (rightstoics));
	return 0;
      }
      if(rev==0){
	irrcols++;(*allrev)=0;
      }
      else
	irrcols+=2;

      for(i=0;i<numleft;i++){
	totchems= addv1(totchems, leftchems[i], chems);
      }
      for(i=0;i<numright;i++){
	totchems= addv1(totchems, rightchems[i], chems);
      }
      freearraydat(leftchems, numleft);
      freearraydat(rightchems, numright);
      free((FREE_ARG) (leftstoics));
      free((FREE_ARG) (rightstoics));
      numlines++;
    }
    free(line);
    line = getlinefromstr(&pos, str);
  }
  free(line);
  (*m)=numlines;
  (*n)=totchems;
  (*cols3)=irrcols;

  (*imat1)=imatrix(0, (*n)-1, 0, (*m)-1);
  (*imat2)=imatrix(0, (*n)-1, 0, (*m)-1);
  (*imat3)=imatrix(0, (*n)-1, 0, irrcols-1);
  (*imat4)=imatrix(0, (*n)-1, 0, irrcols-1);


  //pass two

  pos=0;numlines=0;
  line = getlinefromstr(&pos, str);
  while(strcmp(line, "")!=0){
    if(!iscomline(line) && isreac(line)){
      getreac(line, &leftchems, &leftstoics, &numleft, &rightchems, &rightstoics, &numright, &rev);

      for(i=0;i<(*n);i++){

	lstoic=0;pos1=0;
	while(pos1<numleft){
	  if((ind=isinarray(leftchems+pos1, numleft-pos1, (*chems)[i]))>=0){
	    pos1+=ind;
	    lstoic+=leftstoics[pos1];
	    pos1++;
	  }
	  else
	    pos1=numleft;
	}
	rstoic=0;pos1=0;
	while(pos1<numright){
	  if((ind=isinarray(rightchems+pos1, numright-pos1, (*chems)[i]))>=0){
	    pos1+=ind;
	    rstoic+=rightstoics[pos1];
	    pos1++;
	    //	    fprintf(stderr, "pos1 = %d, chem = %s, stoic = %d\n", ind, rightchems[ind], rstoic);
	  }
	  else
	    pos1=numright;
	}

	if(lstoic>0 && rstoic > 0){ // on both sides of reac
	  (*allgood)=0;
	  (*imat1)[i][numlines]=rstoic-lstoic;
	  if(rev==1){
	    (*imat2)[i][numlines]=2;
	    (*imat3)[i][numcols]=rstoic-lstoic;
	    (*imat4)[i][numcols]=-lstoic;
	    (*imat3)[i][numcols+1]=-rstoic+lstoic;
	    (*imat4)[i][numcols+1]=-rstoic;
	    //	    (*imat4)[i][numcols+1]=-rightstoics[ind];
	  }
	  else{
	    (*imat2)[i][numlines]=-1;
	    (*imat3)[i][numcols]=rstoic-lstoic;
	    (*imat4)[i][numcols]=-lstoic;
	  }

	}
	else if(lstoic>0){ // on left
	  (*imat1)[i][numlines]=-lstoic;
	  (*imat2)[i][numlines]=-1;
	  if(rev==1){
	    (*imat3)[i][numcols]=-lstoic;
	    (*imat4)[i][numcols]=-lstoic;
	    (*imat3)[i][numcols+1]=lstoic;
	    (*imat4)[i][numcols+1]=0;
	  }
	  else{
	    (*imat3)[i][numcols]=-lstoic;
	    (*imat4)[i][numcols]=-lstoic;
	  }
	}
	else if(rstoic>0){ // on right
	  (*imat1)[i][numlines]=rstoic;
	  if(rev==1){
	    (*imat2)[i][numlines]=1;
	    (*imat3)[i][numcols]=rstoic;
	    (*imat4)[i][numcols]=0;
	    (*imat3)[i][numcols+1]=-rstoic;
	    (*imat4)[i][numcols+1]=-rstoic;
	  }
	  else{
	    (*imat2)[i][numlines]=0;
	    (*imat3)[i][numcols]=rstoic;
	    (*imat4)[i][numcols]=0;
	  }
	}

	else{ // doesn't figure in the reaction
	  (*imat1)[i][numlines]=0;
	  (*imat2)[i][numlines]=0;

	  (*imat3)[i][numcols]=0;
	  (*imat4)[i][numcols]=0;
	  if(rev==1){
	    (*imat3)[i][numcols+1]=0;
	    (*imat4)[i][numcols+1]=0;
	  }

	}

      }

      freearraydat(leftchems, numleft);
      freearraydat(rightchems, numright);
      free((FREE_ARG) (leftstoics));
      free((FREE_ARG) (rightstoics));
      numlines++;if(rev==1){numcols+=2;}else{numcols++;}
    }
    free(line);
    line = getlinefromstr(&pos, str);
  }
  free(line);
  return 1;

}


int genMAXMAreacs(char *fname, int **imat1, int **imat2, int n, int m){
  int i, j;
  FILE *fd;

 fd = fopen(fname, "w");
  if(!fd){
    fprintf(stderr, "ERROR in getMAXMAreacs: \"%s\" could not be opened for writing.\n", fname);
    return 0;
  }

  for (j=0;j<m;j++){
    fprintf(fd, "f%d:k%d", j, j);
    for(i=0;i<n;i++){
      if(imat2[i][j]!=0){
	if(imat2[i][j]==-1){
	  fprintf(fd, "*S%d", i);
	}
	else{
	  fprintf(fd, "*(S%d**%d)", i, -imat2[i][j]);
	}
      }
    }
    fprintf(fd, ";\n");
  }
    fprintf(fd, "\n\n");

  for(i=0;i<n;i++){
    fprintf(fd, "F%d: ", i);
    for(j=0;j<m;j++){
      if(imat1[i][j]!=0){
	if(imat1[i][j]==1){
	  fprintf(fd, "+f%d", j);
     
	}
	else if(imat1[i][j]==-1){
	  fprintf(fd, "-f%d", j);
	}
	else if (imat1[i][j]<0){
	  fprintf(fd, "-%d*f%d", -imat1[i][j], j);

	}
	else{
	  fprintf(fd, "+%d*f%d", imat1[i][j], j);

	}

      }
    }
    fprintf(fd, ";\n");

  }
    fprintf(fd, "\n\n");

  fprintf(fd, "jj:matrix(");
  for(i=0;i<n;i++){
    fprintf(fd, "[");
    for(j=0;j<n;j++){
      fprintf(fd, "diff(F%d, S%d)", i, j);
      if(j<n-1){
	fprintf(fd, ", ");
      }
    }
    fprintf(fd, "]");
    if(i<n-1){
      fprintf(fd, ", ");
    }

  }
  fprintf(fd, ");\n\n");

  fclose(fd);
  return 1;


}

int split1(char *s){
  char *p;
  int j=0, ret=1;
  // get the numeric part
  p=strdup(s);
  while(s[j] && isdigit(s[j])){j++;}
  p[j]=0;
  if(!isonlyspace(p))
    ret=atoi(p);
  free(p);
  return ret;
}

char *split2(char *s){
  char *p;
  int j=0;
  // get the character part
  while(s[j] && isdigit(s[j])){j++;}
  p=strdup(s+j);
  return p;
}

int ispureint(char *s){
  // decide if string s is a pure integer (no sign allowed)
  int k=0;
  while(s[k] && isspace(s[k])) //skip space
    k++;
  if(!isdigit(s[k])) // starts with a noninteger or empty (apart from spaces)
    return 0;

  while(s[k] && isdigit(s[k]))
    k++;
  if(s[k] && !isspace(s[k]))
    return 0;
  while(s[k] && isspace(s[k])) //skip space
    k++;
  if(s[k]) // something more
    return 0;

  return 1;
}

int getreac(char *str, char ***leftchems, int **leftstoics, int *numleft, char ***rightchems, int **rightstoics, int *numright, int *rev){
  char *p, *left, *right;
  int j;
  char **v1;
  char **tmp;
  int num1;
  int flag=1;
  if((p=strstr(str, "<-->")) || (p=strstr(str, "<==>"))){
    j=strlen(str)-strlen(p);
    right=strdup(str+j+4);
    left=strdup(str);
    left[j]=0;
    (*rev)=1;
  }
  else if((p=strstr(str, "<->"))|| (p=strstr(str, "<=>"))){
    j=strlen(str)-strlen(p);
    right=strdup(str+j+3);
    left=strdup(str);
    left[j]=0;
    (*rev)=1;
  }
  else if((p=strstr(str, "<>"))){
    j=strlen(str)-strlen(p);
    right=strdup(str+j+2);
    left=strdup(str);
    left[j]=0;
    (*rev)=1;
  }
  else if((p=strstr(str, "-->")) || (p=strstr(str, "==>"))){
    j=strlen(str)-strlen(p);
    right=strdup(str+j+3);
    left=strdup(str);
    left[j]=0;
    (*rev)=0;
  }
  else if((p=strstr(str, "->")) || (p=strstr(str, "=>"))){
    j=strlen(str)-strlen(p);
    right=strdup(str+j+2);
    left=strdup(str);
    left[j]=0;
    (*rev)=0;
  }
  else if((p=strstr(str, ">"))){
    j=strlen(str)-strlen(p);
    right=strdup(str+j+1);
    left=strdup(str);
    left[j]=0;
    (*rev)=0;
  }
  else if((p=strstr(str, "<--")) || (p=strstr(str, "<=="))){
    j=strlen(str)-strlen(p);
    left=strdup(str+j+3);
    right=strdup(str);
    right[j]=0;
    (*rev)=0;
  }
  else if((p=strstr(str, "<-")) || (p=strstr(str, "<="))){
    j=strlen(str)-strlen(p);
    left=strdup(str+j+2);
    right=strdup(str);
    right[j]=0;
    (*rev)=0;
  }
  else if((p=strstr(str, "<"))){
    j=strlen(str)-strlen(p);
    left=strdup(str+j+1);
    right=strdup(str);
    right[j]=0;
    (*rev)=0;
  }
  else{
    fprintf(stderr, "ERROR in reaction %s\n", str);
    flag=0;
    left=strdup("");right=strdup("");
  }
  //left of reaction


  (*numleft)=chemgts2(left, &v1, '+');
  (*leftstoics)=(int*) malloc(sizeof(int) * (*numleft));
  (*leftchems)=(char**) malloc(sizeof(char*) * (*numleft));
  for(j=0;j<(*numleft);j++){
    num1=chemgts2(v1[j], &tmp, ' ');

    if(num1==1){
      (*leftstoics)[j]=split1(v1[j]);(*leftchems)[j]=split2(v1[j]);
    }
    else if(num1==2){
      if(!ispureint(tmp[0])){
	fprintf(stderr, "ERROR in reaction %s\n", str);flag=0;
      }
      (*leftstoics)[j]=atoi(tmp[0]);(*leftchems)[j]=strdup(tmp[1]);
    }
    else{
      (*leftstoics)[j]=0;(*leftchems)[j]=strdup("");
      fprintf(stderr, "ERROR in reaction %s\n", str);flag=0;
    }
    freearraydat(tmp, num1);
  }

  freearraydat(v1, (*numleft));

  //right of reaction

  (*numright)=chemgts2(right, &v1, '+');
  (*rightstoics)=(int*) malloc(sizeof(int) * (*numright));
  (*rightchems)=(char**) malloc(sizeof(char*) * (*numright));
  for(j=0;j<(*numright);j++){
    num1=chemgts2(v1[j], &tmp, ' ');  
  
    if(num1==1){
      (*rightstoics)[j]=split1(v1[j]);(*rightchems)[j]=split2(v1[j]);
      //fprintf(stderr, "%s: rstoic = %d, chem = %s\n", v1[j], (*rightstoics)[j], (*rightchems)[j]);
    }
    else if(num1==2){
     if(!ispureint(tmp[0])){
	fprintf(stderr, "ERROR in reaction %s\n", str);flag=0;
      }
      (*rightstoics)[j]=atoi(tmp[0]);(*rightchems)[j]=strdup(tmp[1]);
    }
    else{
      (*rightstoics)[j]=0;(*rightchems)[j]=strdup("");
      fprintf(stderr, "ERROR in reaction %s\n", str);flag=0;
    }
    freearraydat(tmp, num1);
  }

  freearraydat(v1, (*numright));


  free(left);
  free(right);
  return flag;

}





int isend(char c){
  if((c=='\n') || (c==13) || (c=='\0'))
    return 1;
  return 0;
}

char *lrtrim(char s[]){
  char *p;
  int n;
  while((*s) && isspace(*s)){s++;}// left trim
  p=strdup(s);
  for(n=strlen(p)-1;n>=0;n--)
    if(!(isspace(p[n])))
      break;
  p[n+1] = '\0';
  return p;
}

int chemgts2(char *s, char ***v, char sep){
  // sep is the separator
  int i, j, k;
  int numgets=0;
  i=0, k=0;

  (*v)=NULL;
  while(s[k]){
    j=0;
    while((s[k] == sep) || isspace((int) s[k])){k++;} // skip white space
    while((s[k] != sep) && !isend(s[k])){j++;k++;}
    if(j>0)
      numgets++;
  }
  if(numgets>0){
    (*v)=(char**) malloc(sizeof(char*) * numgets);
    i=0;k=0;
    while(s[k]){
      j=0;
      while((s[k] == sep) || isspace((int) s[k])){k++;} // skip white space
      while((s[k] != sep) && !isend(s[k])){j++;k++;}
      if(j>0)
	(*v)[i++] = lrtrim(strchop2(s, k-j, j));
    }
  }
  return numgets;
}


int freearraydat(char **array, int lim){
  int i=lim;
  if(!array){return 0;}
  while(i-- > 0)
    free(array[i]);
  free((FREE_ARG) (array));
  //free((char*)array);
  return 0;
}

int addv1(int k, char *s, char ***t)
     /* routine to check if a string s already belongs to a list of strings, and if not to add it to the end. Returns new free position in the list. */
{
  int i=0;

  if(k==0){
    (*t) = (char**) malloc(sizeof(char*) * 1);
    (*t)[k] = strdup(s);
    return k+1;
  }

  while(i<k)
    if (strcmp(s, (*t)[i++]) == 0){return k;}

  (*t)=(char**) realloc((*t), sizeof(char*) *(k+1));
  (*t)[k] = strdup(s);
  return k+1;
}

int isinarray(char *v[], int numv, char *s){
  // is the string s a member of the array v? If so return its index, if not return -1
  int i;
  for(i=0;i<numv;i++)
    if(strcmp(s, v[i]) == 0)
      return i;
  return -1;
}

/* check if the matrix stored in file fname */
/* is sign nonsingular or sign singular */

int strmatisSNSSS(char *fname){

  int **imat1;
  char *str;
  int mlen=0, nlen=0;
  int retval;

  str=readfileintostr(fname);
  if(isonlyspace(str)){
    fprintf(stderr, "\n        **ERROR**\nNothing found in file \"%s\". EXITING.\n          *****\n", fname);
    free(str);
    return -1;
  }

  imat1=readmatrixfromstr(str, &nlen, &mlen);
  if(!(*imat1)){
    fprintf(stderr, "ERROR: Couldn't read a matrix from the data in file \"%s\". EXITING. \n", fname);
    free(str);
    return -1;
  }
  fprintf(stderr, "matrix found: \n");
  printmat(imat1, nlen, mlen);

  if(nlen!=mlen){
    fprintf(stderr, "The matrix is not square. Exiting.\n\n");
    exit(-1);
  }

  if((retval=matrixisSNSSS(imat1, nlen))==1){
    fprintf(stderr, "The matrix is sign nonsingular with positive determinant\n");
  }
  else if(retval==-1)
    fprintf(stderr, "The matrix is sign nonsingular with negative determinant\n");
  else if(retval==2)
    fprintf(stderr, "The matrix is sign singular\n");
  else{
    fprintf(stderr, "The matrix is neither sign nonsingular nor sign singular\n");
  }
  free(str);
  free_imatrix(imat1, 0, nlen-1, 0, mlen-1);
  return 0;

}

/* check if the matrix stored in file fname */
/* is SSD */


int strmatisSSD(char *fname, int q){

  int **imat1;
  char *str;
  int mlen=0, nlen=0;
  int retval;

  str=readfileintostr(fname);
  if(isonlyspace(str)){
    fprintf(stderr, "\n        **ERROR**\nNothing found in file \"%s\". EXITING.\n          *****\n", fname);
    free(str);
    return -1;
  }

  imat1=readmatrixfromstr(str, &nlen, &mlen);
  if(!(*imat1)){
    fprintf(stderr, "ERROR: Couldn't read a matrix from the data in file \"%s\". EXITING. \n", fname);
    free(str);
    return -1;
  }
  fprintf(stderr, "matrix: \n");
  printmat(imat1, nlen, mlen);


  if((retval=isSSD(imat1, nlen, mlen, q))){
    fprintf(stderr, "The matrix is SSD.\n");
  }
  else{
    fprintf(stderr, "The matrix is not SSD\n");
  }
  free(str);
  free_imatrix(imat1, 0, nlen-1, 0, mlen-1);
  return 0;

}




/* check if the matrix stored in file fname */
/* is CSD */


int strmatisCSD(char *fname, int q){

  int **imat1;
  char *str;
  int mlen=0, nlen=0;
  int retval;

  str=readfileintostr(fname);
  if(isonlyspace(str)){
    fprintf(stderr, "\n        **ERROR**\nNothing found in file \"%s\". EXITING.\n          *****\n", fname);
    free(str);
    return -1;
  }

  imat1=readmatrixfromstr(str, &nlen, &mlen);
  if(!(*imat1)){
    fprintf(stderr, "ERROR: Couldn't read a matrix from the data in file \"%s\". EXITING. \n", fname);
    free(str);
    return -1;
  }
  fprintf(stderr, "matrix: \n");
  printmat(imat1, nlen, mlen);


  if((retval=isCSD(imat1, nlen, mlen, q))){
    fprintf(stderr, "The matrix is CSD.\n");
  }
  else{
    fprintf(stderr, "The matrix is not CSD\n");
  }
  free(str);
  free_imatrix(imat1, 0, nlen-1, 0, mlen-1);
  return 0;

}


int analysereacs(const char fname[], int q){

  int **imat1, **imat2, **imat3, **imat4;
  char *str;
  int mlen=0, nlen=0;
  char **chems;
  long pos=0;
  char *line;
  int type=-1;
  int allrev, allgood;
  int cols3;
  int csdflag=0;
  int ssdflag=0;

  str=readfileintostr(fname);
  if(isonlyspace(str)){
   fprintf(stderr, "\n        **ERROR**\nNothing found in file \"%s\". EXITING.\n          *****\n", fname);
    free(str);
    return -1;
  }


  line = getlinefromstr(&pos, str);
  while(type==-1 && strcmp(line, "")!=0){
    if(!iscomline(line) && isreac(line))
      type=0; //reaction file
    else if (iscomline(line) && strstr(line, "*****")) // separator character
      type=1;
    free(line);
    line = getlinefromstr(&pos, str);
  }
  free(line);
  pos=0;

  // 
  // The file contains reactions
  //

  if(type==0){
    if(!getallreacs(str, &imat1, &imat2, &imat3, &imat4, &chems, &nlen, &mlen, &cols3, &allrev, &allgood)){
      free(str);
      fprintf(stderr, "ERROR: Couldn't read the reactions in file \"%s\". EXITING. \n", fname);
      return -1;
    }

    if(nlen==0 || mlen ==0){
      fprintf(stderr, "ERROR: Couldn't read the reactions in file \"%s\". EXITING. \n", fname);
      free(str);
      free_imatrix(imat1, 0, nlen-1, 0, mlen-1);
      free_imatrix(imat2, 0, nlen-1, 0, mlen-1);
      free_imatrix(imat3, 0, nlen-1, 0, cols3-1);
      free_imatrix(imat4, 0, nlen-1, 0, cols3-1);
      return -1;
    }
    free(str);


    fprintf(stderr, "The stoichiometric matrix:\n\n");
    printmat1(chems, imat1, nlen, mlen);
    fprintf(stderr, "The pattern matrix for -V^T:\n\n");
    printmatpat(chems, imat2, nlen, mlen);
    fprintf(stderr, "The irreversible stoichiometric matrix:\n\n");
    printmat1(chems, imat3, nlen, cols3);
    fprintf(stderr, "The matrix of powers (for an MA system):\n\n");
    printmat1(chems, imat4, nlen, cols3);

    fprintf(stderr, "\n_________________________________\n\n");

    if(allgood==1){ // no reactants on both sides
      if((csdflag=isCSD(imat1, nlen, mlen, q))){
	fprintf(stderr, "The S-matrix is CSD. The system does not admit MPNE for any kinetics. In fact, this reaction structure does not admit MPNE with *any stoichiometries*, and regardless of whether all reactions are reversible or irreversible.\n");
      }
      else if((ssdflag=isSSD(imat1, nlen, mlen, q))){
	fprintf(stderr, "The S-matrix is SSD. The system does not admit MPNE for any kinetics, regardless of whether the reactions are reversible or irreversible.\n");
      }
      else if(doubleisWSD(imat1, nlen, mlen, q)){
	if(allrev==0 && arecompat(imat1, imat2, nlen, mlen, q)){
	  fprintf(stderr, "The system does not admit MPNE\nfor any kinetics. The system with mass action kinetics does not admit MPNE regardless of whether the reactions are reversible or irreversible.\n");
	}
	else if (allrev==1)
	  fprintf(stderr, "The system with general kinetics may admit MPNE. However, the system with mass action kinetics does not admit MPNE.\n");
      }
      else{
	if(allrev==0){ // some irreversible
	  if(arecompat(imat1, imat2, nlen, mlen, q)){
	    fprintf(stderr, "The S-matrix and pattern -V^T are compatible. The system does not admit MPNE for any kinetics.\n");
	  }
	  else if(mats_compat(imat3, imat4, nlen, cols3, q)){
	    fprintf(stderr, "The system may admit MPNE for general kinetics, but\nwith mass-action kinetics it does not admit MPNE.\n");
	  }
	  else{
	    allminorsigns(imat3, imat4, nlen, cols3, q);
	    fprintf(stderr, "The system may admit MPNE both for mass-action kinetics and for general kinetics.\n");
	  }

	}
	else{
	  fprintf(stderr, "The system may admit MPNE both for mass-action kinetics and for general kinetics.\n");
	}

      }
    }
    else{ // some reactants on both sides (doesn't matter if reversible or not)
      if(arecompat(imat1, imat2, nlen, mlen, q)){
	fprintf(stderr, "The S-matrix and pattern -V^T are compatible. The system does not admit MPNE for any kinetics.\n");
      }
      else if(mats_compat(imat3, imat4, nlen, cols3, q)){
	fprintf(stderr, "The system may admit MPNE for general kinetics, but\nwith mass-action kinetics it does not admit MPNE.\n");
      }
      else{
	allminorsigns(imat3, imat4, nlen, cols3, q);
	fprintf(stderr, "The system may admit MPNE both for mass-action kinetics and for general kinetics.\n");
      }


    }




    free_imatrix(imat1, 0, nlen-1, 0, mlen-1);
    free_imatrix(imat2, 0, nlen-1, 0, mlen-1);
    free_imatrix(imat3, 0, nlen-1, 0, cols3-1);
    free_imatrix(imat4, 0, nlen-1, 0, cols3-1);
  }
  else if (type==1){
    if(!readmatpairfromstr(str, &nlen, &mlen, &imat1, &imat2)){ 
      fprintf(stderr, "ERROR: Expecting two matrices in file \"%s\". Couldn't find these. EXITING. \n", fname);
      free(str);
      return -1;
    }
    fprintf(stderr, "Assuming that these are the stoichiometric matrix and the pattern matrix for -V^T:\n\n");
    printmat(imat1, nlen, mlen);
    printmat(imat2, nlen, mlen);

    if(arecompat(imat1, imat2, nlen, mlen, q)){
      fprintf(stderr, "The system does not admit MPNE\nfor any kinetics.\n");
    }
    else{
      allminorsigns(imat1, imat2, nlen, mlen, q);
      fprintf(stderr, "The system may admit MPNE for some kinetics.\n");
    }
    free_imatrix(imat1, 0, nlen-1, 0, mlen-1);
    free_imatrix(imat2, 0, nlen-1, 0, mlen-1);

  }
  else{
    imat1=readmatrixfromstr(str, &nlen, &mlen);
    if(!(*imat1)){
      fprintf(stderr, "ERROR: Couldn't find reactions, a pair of matrices or a matrix in file \"%s\". EXITING. \n", fname);
      free(str);
      free_imatrix(imat1, 0, nlen-1, 0, mlen-1);
      return -1;
    }
    fprintf(stderr, "Assuming that this is a stoichiometric matrix, and that all reactions may or may not be irreversible.\n");
    printmat(imat1, nlen, mlen);
    if((csdflag=isCSD(imat1, nlen, mlen, q))){
      fprintf(stderr, "The system does not admit MPNE for any kinetics. In fact, this reaction structure does not admit MPNE with *any stoichiometries*, and regardless of whether all reactions are reversible or irreversible.\n");
    }
    else if((ssdflag=isSSD(imat1, nlen, mlen, q))){
      fprintf(stderr, "The system does not admit MPNE for any kinetics, regardless of whether the reactions are reversible or irreversible.\n");
    }
    else if(doubleisWSD(imat1, nlen, mlen, q)){
	  
      fprintf(stderr, "The system with general kinetics may admit MPNE. However, the system with mass action kinetics does not admit MPNE.\n");
    }
    else{
      fprintf(stderr, "Without further information on reversibility, etc. the system may admit MPNE.\n");
    }


    free_imatrix(imat1, 0, nlen-1, 0, mlen-1);
  }



  return 0;


}


bool unordllistareeq(long *a, int na, long *b, int nb){
  int i;
  if(na!=nb)
    return 0;
  for(i=0;i<na;i++){
    if(!isinllist(a[i],b,nb))
       return 0;
  }
  return 1;
}

bool unordlistareeq(int *a, int na, int *b, int nb){
  int i;
  if(na!=nb)
    return 0;
  for(i=0;i<na;i++){
    if(!isinlist(a[i],b,nb))
       return 0;
  }
  return 1;
}


int areequal(int *vec1, int *vec2, int n){
  int i;
  for(i=0;i<n;i++){
    if(vec2[i]!=vec1[i])
      return 0;
  }
  return 1;
}


// is vec2 > vec1 read backwards

int isgrtr(int *vec1, int *vec2, int n){
  int i;
  for(i=0;i<n;i++){
    if(vec2[n-1-i]>vec1[n-1-i])
      return 1;
    else if(vec2[n-1-i]<vec1[n-1-i])
      return -1;
  }
  return 0;
}

// is vec2 > vec1 read forwards

int isgrtrf(int *vec1, int *vec2, int n){
  int i;
  for(i=0;i<n;i++){
    if(vec2[i]>vec1[i])
      return 1;
    else if(vec2[i]<vec1[i])
      return -1;
  }
  return 0;
}

void iswapi(int **v, int i, int j){
  int *temp = v[i];
  v[i]=v[j];
  v[j]=temp;
}

// sort a list of integer vectors (monomials)

void qsorti(int **v, int dim, long left, long right)
{
  long i, last;

  if (left >= right)
    return;
  iswapi(v,left,(left + right)/2);
  last = left;
  for (i=left+1;i<=right;i++)
    if (isgrtrf(v[i],v[left],dim)>0)
      iswapi(v,++last,i);
  iswapi(v,left,last);
  qsorti(v, dim, left, last-1);
  qsorti(v, dim, last+1, right);
}

// sort a list of integer vectors (monomials) and a list of
// coefficients

void qsorti2(int **v, int *cfs, int dim, long left, long right)
{
  long i, last;

  if (left >= right)
    return;
  /* cout << "swapping v[" << left << "] and v[" << (left + right)/2 << endl; */
  /* printvec1(v[left], dim); */
  /* printvec1(v[(left + right)/2], dim); */
  iswapi(v,left,(left + right)/2);
  /* cout << "now: " << endl; */
  /* printvec1(v[left], dim); */
  /* printvec1(v[(left + right)/2], dim); */
 

  iswap(cfs,left,(left + right)/2);
  last = left;
  for (i=left+1;i<=right;i++){
    if (isgrtrf(v[i],v[left],dim)>0){
      ++last;
      iswapi(v,last,i);
      iswap(cfs,last,i);
    }
  }
  iswapi(v,left,last);
  iswap(cfs,left,last);
  qsorti2(v, cfs, dim, left, last-1);
  qsorti2(v, cfs, dim, last+1, right);
}

// v is an ordered list of integer vectors each of length vlen
// there are n vectors in v. The algorithm perfoms a 
// binary search on v to find vector x

long binisearch(int **v, long n, int *x, int vlen){
  long low, high, mid;
  int cmp;
  low=0;
  high=n-1;
  while(low<=high){
    mid=(low+high)/2;
    if((cmp=isgrtrf(x,v[mid],vlen))>0){ // v[mid] > x
      /* printvec1(x,vlen); */
      /* printvec1(v[mid],vlen); */
      /* cout << "greater" << endl; */
      high=mid-1;
    }
    else if(cmp<0) // x > v[mid]
      low=mid+1;
    else // equal
      return mid;
  }
  return -1;
}


int isinarray1(int **v, int numv, int *s, int n){
  // is the integer vector s a member of the array v? If so return its index, if not return -1
  int i;
  // assume first two entries are other data
  for(i=0;i<numv;i++){
    if(areequal(v[i]+2, s, n))
      return i;
  }
  return -1;
}

long isinarray2(int **v, long numv, int *s, int n){
  // is the integer vector s a member of the array v? 
  long i;
  for(i=0;i<numv;i++){
    if(areequal(v[i], s, n))
      return i;
  }
  return -1;
}

int isinarray3(int **v, int numv, int *s, int n){
  // is the integer vector s a member of the array v? 
  int i;
  for(i=0;i<numv;i++){
    if(areequal(v[i], s, n))
      return i;
  }
  return -1;
}

void firstcomb(int *vec, int n, int n1){
  int i;
  if(n1>n){
    for(i=0;i<n1;i++)
      vec[i]=0;
  }
  else{
    for(i=0;i<n1;i++)
      vec[i]=i;
  }
  return;
}

int firstcombfrom(int *vec, int n, int n1, int i_init){
  int i;
  //fprintf(stderr, "n=%d, n1=%d, i_init = %d\n", n, n1, i_init);
  if(n<i_init+n1){
    for(i=0;i<n1;i++)
      vec[i]=0;
    //fprintf(stderr, "exiting here\n");
    return 0;
  }
  else{
    for(i=i_init;i<i_init+n1;i++)
      vec[i-i_init]=i;
  }
  //fprintf(stderr, "or here\n");
  return 1;
}


int nextcomb(int *vec, int n, int n1){
  int i, j;
  for(i=0;i<n1;i++){
    if(vec[n1-1-i]< n-1-i){
      vec[n1-1-i]++;
      for(j=n1-i;j<n1;j++){
	vec[j]=vec[j-1]+1;
      }
      return 1;
    }
  }
  // only get here if we fail
  firstcomb(vec, n, n1);
  return 0;
}

void nextcombk(int *vec, int n, int n1, long k){
  long i;
  for(i=0;i<k;i++){
    nextcomb(vec, n, n1);
  }
  return;
}

void nextnum(int *vec, int n, int base){
  int i;
  for(i=0;i<n;i++){
    if(vec[n-1-i]< base-1){
      vec[n-1-i]++;
      return;
    }
    vec[n-1-i]=0;
  }
}

// each digit has a different base

void nextnumb(int *vec, int n, int *base){
  int i;
  for(i=0;i<n;i++){
    if(vec[n-1-i]< base[n-1-i]-1){
      vec[n-1-i]++;
      return;
    }
    vec[n-1-i]=0;
  }
}



// a number of length n in some base to base 10
long base10(int *vec, int n, int base){
  long j=0;
  int i=0;
  for(i=0;i<n;i++){
    j+=vec[i]*base^(n-i-1);
  }
  return j;
}

// My implementation of Johnson-Trotter algorithm

int nextperm(int **vec, int **veclr, int *par, int n){
  int tmp=-1, tmpi=-1, tmplr;
  int i, flag=0;
  for(i=0;i<n;i++){
    if((*vec)[i]>tmp){
      if((*veclr)[i]==1 && i<n-1 && ((*vec)[i] > (*vec)[i+1])){ // right
	tmp = (*vec)[i];
	tmpi=i;
	tmplr=(*veclr)[i];
	flag=1;
      }
      else if ((*veclr)[i]==-1 && i>0 && ((*vec)[i] > (*vec)[i-1])){ // left
	tmp = (*vec)[i];
	tmpi=i;
	tmplr=(*veclr)[i];
	flag=2;
      }
    }
  }
  if(flag==1){
    (*vec)[tmpi]=(*vec)[tmpi+1];
    (*vec)[tmpi+1]=tmp;
    (*veclr)[tmpi]=(*veclr)[tmpi+1];
    (*veclr)[tmpi+1]=tmplr;
    *par = -(*par);
    for(i=0;i<n;i++){
      if((*vec)[i]>tmp)
	(*veclr)[i]=-(*veclr)[i];
    }
    return 1;
  }
  if(flag==2){
    (*vec)[tmpi]=(*vec)[tmpi-1];
    (*vec)[tmpi-1]=tmp;
    (*veclr)[tmpi]=(*veclr)[tmpi-1];
    (*veclr)[tmpi-1]=tmplr;
    *par = -(*par);
    for(i=0;i<n;i++){
      if((*vec)[i]>tmp)
	(*veclr)[i]=-(*veclr)[i];
    }
    return 1;
  }
  return 0;

}

/* Just for testing purposes: prints out the */
/* permutations on [1, ..., n] in transposition order */

int tr(){
  int *vec;
  int *veclr;
  int i, par=1, flag=1,t=0;
  int n=4;
  vec=(int *)malloc((size_t) ((n)*sizeof(int)));
  veclr=(int *)malloc((size_t) ((n)*sizeof(int)));

  for(i=0;i<n;i++){
    vec[i]=i;
    veclr[i]=-1;
  }

  while(flag){
    for(i=0;i<n;i++){
      if(veclr[i]==-1)
	fprintf(stderr, "<%d  ", vec[i]+1);
      else
	fprintf(stderr, "%d>  ", vec[i]+1);
    }
    fprintf(stderr, "(%d)\n", par);
    flag=nextperm(&vec, &veclr, &par, n);
    t++;
  }
  return 0;

}

// a number in base 10 to a number of length n in some other base
void basek(int *vec, int n, int base, long num){
  int i;
  for(i=0;i<n;i++){
    vec[i]=num/base^(n-i-1);
    num-=vec[i]*base^(n-i-1);
    if(i==0)
      vec[i]=vec[i]%base;
  }
  return;
}

// quick version of adding in base "base", k is in base 10.

void nextnumk(int *vec, int n, int base, long k){
  long tmp;
  tmp=base10(vec, n, base);
  tmp+=k;
  basek(vec, n, base, tmp);
}

void order(int *vec, int n){
  int i, tmp, flag=1;
  while(flag==1){
    flag=0;
    for(i=0;i<n-1;i++){
      if(vec[i]>vec[i+1]){
	tmp=vec[i];vec[i]=vec[i+1];vec[i+1]=tmp;
	flag=1;
      }
    }
  }

}

// n is the length of the vector
// m is twice the number of reactions
// the vector output has the following structure:
// the first entry is the number of entries
// the second entry is the sign of the relevant minor
// the next m entries are a numerical representation of the 
// symbolic term: e.g. if we write V as one long column vector,
// (first column, followed by second column, etc.)
// from x_0 to x_nm, then 2 2 4 means (x_2)^2 x_4
// each subsequent block of m entries is a minor

int **getallinds(int *longv, int *longv1, int *base, int *basetot, long numtot, int *comb, int n, long *k){
  long i;
  int j, ind, l;
  int vec[n];
  int term[n];
  int val;
  int **outmat=NULL;
  int bs[n];
  //initialise
  for(i=0;i<n;i++)
    vec[i]=0;

  for(i=0;i<n;i++) //the base for increment...
    bs[i]=base[comb[i]];


  for(i=0;i<numtot;i++){
    //    if(i%1000==0){fprintf(stderr, "i/1000 = %.0f\n", (double) (i/1000));}


    for(j=0;j<n;j++){
      term[j]=longv[basetot[comb[j]] + vec[j]];
    }
    order(term, n);
    //      printvec(term, n);

    if((ind=isinarray1(outmat, (*k), term, n))>=0){ // in array
      outmat[ind]=(int *)realloc(outmat[ind], sizeof(int*) *(outmat[ind][0]+n));
  
      for(l=0;l<n;l++)
	outmat[ind][outmat[ind][0]+l]=basetot[comb[l]]+vec[l];
      outmat[ind][0]+=n; // length of vector

    }
    else{ // new entry
      if((*k)==0){ // first entry
	outmat = (int**) malloc(sizeof(int*) * 1);
	(*k)++;
      }
      else{
	outmat=(int**) realloc(outmat, sizeof(int*) *((*k)+1));
	(*k)++;
      }
      outmat[(*k)-1]=(int *)malloc(sizeof(int*) *(2*n+2));
      outmat[(*k)-1][0]=2*n+2; // first entry is number of entries

      val=1;
      for(l=0;l<n;l++){
	val*=longv1[basetot[comb[l]] + vec[l]];
      }
      outmat[(*k)-1][1]=val; // second entry is value
      // next n entries are the "symbolic" quantities
      for(l=0;l<n;l++)
	outmat[(*k)-1][l+2]=term[l];
      // next n entries are the indices
      for(l=0;l<n;l++)
	outmat[(*k)-1][n+l+2]=basetot[comb[l]]+vec[l];

      if((*k)%1000==0){fprintf(stderr, "k/1000 = %.0f\n", (double) ((*k)/1000));}
    }

    nextnumb(vec, n, bs);
  }
  return outmat;
}

void printindvec(char *s, int *v1, int *v2, int k){
  int j1, j2;

  for(j1=0;j1<k;j1++){
    for(j2=0;j2<k;j2++){
      fprintf(stderr, "%s[%d,%d] ", s, v1[j1], v2[j2]);
    }
    fprintf(stderr, "\n");
  }
  fprintf(stderr, "\n");
}


// check whether two matrices have sign patterns which are compatible

int checkpat(int **imat1, int **imat2, int n, int m){
// sub: imat2 is a subpattern of imat1
// sup: imat1 is a subpattern of imat2
  int i,j,flageq=1,flagsub=1,flagsup=1; 
  for(i=0;i<n;i++){
    for(j=0;j<m;j++){
      if(imat1[i][j]*imat2[i][j]<0)
	return 0; // short cycles: a special case
    }
  }
  for(i=0;i<n;i++){
    for(j=0;j<m;j++){
      if(imat2[i][j]==0 && imat1[i][j]!=0){
	flageq=0;flagsup=0;
      }
      else if(imat1[i][j]==0 && imat2[i][j]!=0){
	flageq=0;flagsub=0;
      }
      if(!flageq&&!flagsub&&!flagsup)
	return -2; // neither superset not subset
    }
  }
  if(flageq)
    return 1;
  if(flagsub)
    return 2;
  if(flagsup)
    return -1;
  return -2;
}

int **joinmats(int **imat1, int **imat2, int n, int m){
  int **newmat=NULL;
  int i, j;
  newmat=imatrix(0, n-1, 0, m-1);
  for(i=0;i<n;i++){
    for(j=0;j<m;j++){
      if(imat2[i][j]>0 || imat1[i][j]>0)
	newmat[i][j]=1;
      else if(imat1[i][j]<0 || imat2[i][j]<0)
	newmat[i][j]=-1;
      else
	newmat[i][j]=0;
    }
  }
  return newmat;
}

bool insetint(int i, int *j, int sz, int k){
  int m;
  if(k==i)
    return 1;
  for(m=0;m<sz;m++){
    if(k==j[m])
      return 1;
  }
  return 0;
}
void addon(int *v1, int *v2, int n, int k){
  int i;
  for(i=0;i<n;i++)
    v2[i]=v1[i];
  v2[n]=k;
  return;
}

bool rowadjto(bool **p, int j, int *c, int dim){
  int r;
  for(r=0;r<dim;r++)
    if(p[j][c[r]])
      return 1;
  return 0;
}

bool coladjto(bool **p, int j, int *c, int dim){
  int r;
  for(r=0;r<dim;r++)
    if(p[c[r]][j])
      return 1;
  return 0;
}

bool **inttobool(int **imat, int n, int m){
  int i, j;
  bool **bmat=bmatrix(0,n-1,0,m-1);
  for(i=0;i<n;i++){
    for(j=0;j<m;j++){
      if(imat[i][j])
	bmat[i][j]=1;
      else
	bmat[i][j]=0;
    }
  }
  return bmat;
}

bool **mattobool(matrix mmat, int n, int m){
  int i, j;
  bool **bmat=bmatrix(0,n-1,0,m-1);
  for(i=0;i<n;i++){
    for(j=0;j<m;j++){
      if(mmat(i,j)!=0)
	bmat[i][j]=1;
      else
	bmat[i][j]=0;
    }
  }
  return bmat;
}

int **mattoint(matrix mmat, int n, int m){
  int i, j;
  int **imat=imatrix(0,n-1,0,m-1);
  for(i=0;i<n;i++){
    for(j=0;j<m;j++){
      if(mmat(i,j)!=0)
	imat[i][j]=1;
      else
	imat[i][j]=0;
    }
  }
  return imat;
}



// take a list of integers and a 
// list of integer lists
// for each existing list check if each integer belongs to it
// and if not create a new list which is the concatenation
// of the existing list and the new integer
// The aim is to count matchings

int **addintstostrs(long *k, int *newint, int numnew, int **t, int numstrs, int len)  
{
  int i=0,j=0, m=0,flag;
  long r=0;
  int **out=NULL;

  if(numstrs==0){ // no strings
    out = (int**) malloc(sizeof(int*) * numnew); 
    for(i=0;i<numnew;i++){ // add each new int
      out[i] = (int*) malloc(sizeof(int) * 1);
      out[i][0] = newint[i];
    }
    (*k)=numnew;
    return out;
  }

  while(j<numstrs){ //existing strings
    for(m=0;m<numnew;m++){ // new integers
      flag=1;
      for(i=0;i<len;i++){ // is newint[m] in existing string j?
	if(t[j][i]==newint[m]){
	  flag=0; // yes it is
	  break;
	}
      }
      if(flag){// no it isn't so create concatenated str
	if(r==0){
	  out = (int**) malloc(sizeof(int*) * 1);
	  out[r] = (int*) malloc(sizeof(int) * (len+1));
	}
	else{
	  out=(int**) realloc(out, sizeof(int*) *(r+1));
	  out[r] = (int*) malloc(sizeof(int) * (len+1));
	}
	for(i=0;i<len;i++)
	  out[r][i]=t[j][i];
	out[r][len]=newint[m];
	r++;
      }
    }
    j++;
  }
  (*k)=r;
  return out;
}

// level is the size of the lists
// numlists is the number of lists
// k is the column number 
// r is the number of new lists

int **growlists(int **lists, int level, int numlists, int *rowset, int setlen, int k, bool **bmat, long *r){
  int j,p=0;
  int **newlists=NULL;
  (*r=0);
  int rs[setlen];

  for(j=0;j<setlen;j++){
    if(bmat[rowset[j]][k])
      rs[p++]=j;
  }
  /* cout << "p=" << p << endl; */
  /* cout << "numlists = " << numlists << endl; */
  
  newlists= addintstostrs(r, rs, p, lists, numlists, level);
  /* for(j=0;j<(*r);j++) */
  /*   printvec(newlists[j], level+1); */

  return newlists;
}


//
//for a square bool matrix
//

int **getalllists(int dim, bool **bmat, long *r){
  int **lists=NULL;
  int **newl=NULL;
  int i;
  long j, tot;
  int *xc=(int *)malloc((size_t) (dim*sizeof(int)));
  firstcomb(xc, dim, dim);
  *r=0;
  tot=0;
  for(i=0;i<dim;i++){
    cout << "i = " << i << endl;
    cout << "*r = " << *r << endl;
    newl=growlists(lists, i, tot, xc, dim, i, bmat, r);
    if(lists){
      for(j=0;j<tot;j++)
	free((char *) lists[j]);
      free((char *) lists);
    }
    lists=newl; // reset pointer
    tot=(*r);
  }
  return lists;
}

 

int permpar(int *p, int n){
  int transcount=0;
  int pcpy[n];
  int i, j, p0;
  for(i=0;i<n;i++)
    pcpy[i]=p[i];
  for(i=0;i<n;i++){
    p0=pcpy[i];
    if(p0!=i){
      for(j=i+1;j<n;j++){
	if(pcpy[j]==i){
	  pcpy[i]=i;pcpy[j]=p0;
	  transcount++;
	  break;
	}
      }
    }
  }
  if(transcount%2==0)
    return 1;
  return -1;
}



struct inode
{
  int *iarr;
  int cff;
  struct inode *left;
  struct inode *right;
};

struct inode *talloc(void){
  return (struct inode *)malloc(sizeof(struct inode));
}


int *idup(int *i, int n){
  int j;
  int *k=(int *)malloc((size_t) (n*sizeof(int)));
  for(j=0;j<n;j++)
    k[j]=i[j];
  return k;
}

//
// extracts the coefficient only once
//

struct inode *addtree(struct inode *p, int *w, int n, int *rclist, int dim, int ****J2){
  int cond;
  if(p==NULL){
    p=talloc();
    p->iarr = idup(w,n);
    p->cff=coeff_int(J2, rclist, rclist, dim, w);
    p->left=p->right=NULL;
    //fprintf(stderr, "%d:  ", p->cff);printvec(w,n);
  }
  else if((cond=cmpare(w, p->iarr,n))==0){}// already there
  else if(cond<0)
    p->left=addtree(p->left,w,n,rclist,dim,J2);
  else
    p->right=addtree(p->right,w,n,rclist,dim,J2);
  return p;
}

//
// augments the coefficient
//

struct inode *addtree1(struct inode *p, int cf, int *w, int n){
  int cond;
  if(p==NULL){
    p=talloc();
    p->iarr = idup(w,n);
    p->cff=cf;
    p->left=p->right=NULL;
    //fprintf(stderr, "%d:  ", p->cff);printvec(w,n);
  }
  else if((cond=cmpare(w, p->iarr,n))==0){
    p->cff+=cf;
  }// already there
  else if(cond<0)
    p->left=addtree1(p->left,cf,w,n);
  else
    p->right=addtree1(p->right,cf,w,n);
  return p;
}



void treeprint(struct inode *p, int n){
  if(p!=NULL){
    treeprint(p->left, n);
    fprintf(stderr, "%d,  ", p->cff);
    printvec(p->iarr, n);
    treeprint(p->right, n);
  }
}

long treeread(char *fname, int ***list, int **cff, int *ord){
  FILE *fd;
  char oneline[1000];
  long k=0;
  int j;
  int len=0;
  char *tmp;
  (*ord)=0;
  *list=NULL;
  fd = fopen(fname, "r");
  if(!fd){
    fprintf(stderr, "ERROR in treeread: \"%s\" could not be opened for reading.\n", fname);
    exit(0);
  }
  //cout << "here\n";
  while(getline0(fd, oneline, 1000) > 0){
    if(!(iscomline(oneline))){
      tmp=oneline;
      if(k==0){
	//cout << "k = " << k << endl;
	(*ord)=getnumints(oneline);
	if(!(*ord)){
	  fprintf(stderr, "ERROR in treeread: not a valid tree - no integers found on the first line:\n%s", oneline);
	  exit(0);
	}
	else if((*ord)==1){
	  fprintf(stderr, "ERROR in treeread: not a valid tree - only one integers found on the first line:\n%s", oneline);
	  exit(0);
	}

	(*list)=(int**) malloc(sizeof(int*));
	(*cff)=(int*) malloc(sizeof(int));
	(*list)[k]=(int*) malloc(sizeof(int)*((*ord)-1));
	//cout << "k = " << k << endl;
	(*cff)[k]=getint(&tmp);
	//cout << "ord = " << *ord << endl;
	for(j=0;j<(*ord)-1;j++){
	  //cout << "tmp = \"" << tmp << "\"" << endl;
	  (*list)[k][j]=getint(&tmp);

	}
	//cout << "k = " << k << endl;
      }
      else{
	len=getnumints(oneline);
	if(len!=(*ord)){
	  fprintf(stderr, "ERROR in treeread: not a valid tree - bad line:\n%s", oneline);
	  exit(0);
	}
	(*list)=(int**) realloc((*list), sizeof(int*)*(k+1));
	(*list)[k]=(int*) malloc(sizeof(int)*(len-1));
	(*cff)=(int*) realloc((*cff), sizeof(int)*(k+1));
	(*cff)[k]=getint(&tmp);
	for(j=0;j<len-1;j++)
	  (*list)[k][j]=getint(&tmp);
      }
      k++;
    }
  }
  fclose(fd);
  return k;
}

void listprint(int *cfs, int **mons, int dim, long len){
  long k;
  for(k=0;k<len;k++){
    cout << cfs[k] << ", ";
    printvec1(mons[k],dim);
  }
}




int idiff(int *a, int *b, int n){
  int k=0, k1=0;
  int tot=0;
  while(k<n && k1<n){
    if(a[k]<b[k1]){
      tot++;k++;
    }
    else if(a[k]>b[k1]){
      tot++;k1++;
    }
    else{k++;k1++;}
  }
  return tot;
}

bool idiffis1(int *a, int *b, int n){
  int k=0, k1=0;
  int tot=0;
  while(k<n && k1<n){
    if(a[k]<b[k1]){
      tot++;if(tot>1){return 0;}k++;
    }
    else if(a[k]>b[k1]){
      tot++;if(tot>1){return 0;}k1++;
    }
    else{k++;k1++;}
  }
if(!tot) // same string
    return 0;
  return 1;
}

// compute a determinant given the list of matchings

ex detterm(matrix mmat, int n, int m, int *xc, int *yc, int **lists, int numlists, int dim){
  ex out, tmp;
  int i,j;
  for(i=0;i<numlists;i++){
    tmp=1;
    for(j=0;j<dim;j++)
      tmp*=mmat(xc[lists[i][j]],yc[j]);
    out+=permpar(lists[i], dim)*tmp;
    //out+=tmp;
  }
  return out;
}


//
// identity matrix
//

int **makeidmat(int n){
  int **imat=imatrix(0,n-1,0,n-1);
  int i,j;
  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      if(i==j)
	imat[i][j]=1;
      else
	imat[i][j]=0;
    }
  }
  return imat;
}


// Reorder the columns of a bool matrix putting the sparsest first

bool **reordercolsparse(bool **bmat, int n, int m){
  bool **outmat=bmatrix(0, n-1, 0, m-1);
  int i,j,k,c,r[m],s[m]; // r[m] are column totals, s[m] is the new order
  for(i=0;i<m;i++){ //each col
    r[i]=0;
    for(j=0;j<n;j++){ // each row
      if(bmat[j][i])
	(r[i])++;
    }
  }
  c=0;k=0;
  while(c<n+1){
    for(i=0;i<m;i++){ // each col
      if(r[i]==c) // total=c
	s[k++]=i;
    }
    c++;
  }
  printvec(s, m);
  for(j=0;j<n;j++)
    for(i=0;i<m;i++)
      outmat[j][i]=bmat[j][s[i]];

  return outmat;
}

// Reorder the columns of a bool matrix, an int matrix and 
// an ex matrix putting the sparsest first

void multireordercolsparse(bool **bmatin, bool ***bmatout, bool **bmat2in, bool ***bmat2out, int **imatin, int ***imatout, matrix mmatin, matrix *mmatout, int n, int m){
  int i,j,k,c,r[m],s[m]; // r[m] are column totals, s[m] is the new order
  (*bmatout)=bmatrix(0, n-1, 0, m-1);
  (*bmat2out)=bmatrix(0, n-1, 0, m-1);
  (*imatout)=imatrix(0, n-1, 0, m-1);
  for(i=0;i<m;i++){ //each col
    r[i]=0;
    for(j=0;j<n;j++){ // each row
      if(bmatin[j][i])
	(r[i])++;
    }
  }
  c=0;k=0;
  while(c<n+1){
    for(i=0;i<m;i++){ // each col
      if(r[i]==c) // total=c
	s[k++]=i;
    }
    c++;
  }
  fprintf(stderr, "reordering vector:\n");
  printvec(s, m);
  for(j=0;j<n;j++){
    for(i=0;i<m;i++){
      (*bmatout)[j][i]=bmatin[j][s[i]];
      (*bmat2out)[j][i]=bmat2in[j][s[i]];
      (*imatout)[j][i]=imatin[j][s[i]];
      (*mmatout)(j,i)=mmatin(j,s[i]);
    }
  }
  return;
}



int **imatfromiimat(int ***iimat, int n1, int m1){
  int **imat=imatrix(0, n1-1, 0, m1-1);
  int i,j;
  for(i=0;i<n1;i++){
    for(j=0;j<m1;j++){
      if(iimat[i][j][0]==0)
	imat[i][j]=0;
      else
	imat[i][j]=iimat[i][j][1];
    }
  }
  return imat;
}

// The routine converts an integer matrix to a symbolic matrix
// maxv stores the largest variable subscript

matrix Vmatfromimat(int **imat, int n1, int m1, int *maxv){
  matrix exmat(n1, m1);
  char str1[5];
  int i,j;
  (*maxv)=0;
  for(i=0;i<n1;i++){
    for(j=0;j<m1;j++){
      if(imat[i][j] !=0){
	if(abs(imat[i][j])>(*maxv)) //largest variable subscript
	  (*maxv)=abs(imat[i][j]);

  	if(imat[i][j] >0){
  	  sprintf(str1, "v%d", imat[i][j]);
  	  exmat(i,j)=get_possymbol(str1);
  	}
  	else if(imat[i][j] <0){
  	  sprintf(str1, "v%d", -imat[i][j]);
  	  exmat(i,j)=-get_possymbol(str1);
  	}
      }
    }
  }
  return exmat;
}

matrix Vmatfromiimat(int ***imat, int n1, int m1, int *maxv){
  matrix exmat(n1, m1);
  char str1[5];
  int i,j,k;
  (*maxv)=0;
  for(i=0;i<n1;i++){
    for(j=0;j<m1;j++){
      exmat(i,j)=0;
      for(k=1;k<=imat[i][j][0];k++){
	if(imat[i][j][k] !=0){
	  if(abs(imat[i][j][k])>(*maxv)) //largest variable subscript
	    (*maxv)=abs(imat[i][j][k]);

	  if(imat[i][j][k] >0){
	    sprintf(str1, "v%d", imat[i][j][k]);
	    exmat(i,j)+=get_possymbol(str1);
	  }
	  else if(imat[i][j][k] <0){
	    sprintf(str1, "v%d", -imat[i][j][k]);
	    exmat(i,j)-=get_possymbol(str1);
	  }
	}
      }
    }
  }
  return exmat;
}

int isPmatrix(matrix J, int n){
  int *xc;
  int k;
  int flag;
  ex detex;
  for(k=1;k<=n;k++){
    xc=(int *)malloc((size_t) (k*sizeof(int)));
    firstcomb(xc, n, k);flag=1;
    

    while(flag==1){
//    printvec1(xc,k);
      detex = symbdetsubmat(J, n, n, xc, xc, k);
      printvec1(xc,k);
      if(detex!=0 && !(detex.info(info_flags::positive))){
	cout << detex << endl;
	printexsubmat(J,xc,xc,k,k);
	cout << "An apparently nonpositive principal minor found.\n";
	return 0;
      }
      flag=nextcomb(xc, n, k);
    }
    free((char *) xc);
  }
  return 1;
  
}

bool hassignedpminors(matrix J, int n){
  int *xc;
  int k;
  int flag;
  ex detex,ndetex;
  bool outf=1;
  for(k=1;k<=n;k++){
    xc=(int *)malloc((size_t) (k*sizeof(int)));
    firstcomb(xc, n, k);flag=1;
    

    while(flag==1){
      //printvec1(xc,k);
      detex=symbdetsubmat(J, n, n, xc, xc, k);
      ndetex=-detex;
      printvec1(xc,k);
      if(detex==0)
      	cout << "0" << endl;
      else if(detex.info(info_flags::positive))
      	cout << "+" << endl;
      else if(ndetex.info(info_flags::positive))
	cout << detex << endl;
      //cout << "-" << endl;
      /* if(detex!=0 && !(detex.info(info_flags::positive))&& !(ndetex.info(info_flags::positive))){ */
      /* 	cout << detex << endl; */
      /* 	printexsubmat(J,xc,xc,k,k); */
      /* 	cout << "An apparently unsigned principal minor found.\n"; */
      /* 	outf=0; */
      /* } */


      /* else if(detex==0) */
      /* 	cout << "0" << endl; */
      /* else if(detex.info(info_flags::positive)) */
      /* 	cout << "+" << endl; */
      /* else if(detex.info(info_flags::negative)) */
      /* 	cout << "-" << endl; */

      flag=nextcomb(xc, n, k);
    }
    free((char *) xc);
  }
  return outf;
  
}

bool hassignedminors(matrix J, int n){
  int *xc,*yc;
  int k;
  int flag,flag1;
  ex detex,ndetex;
  bool outf=1;
  for(k=n-1;k<=n-1;k++){
    xc=(int *)malloc((size_t) (k*sizeof(int)));
    yc=(int *)malloc((size_t) (k*sizeof(int)));
    firstcomb(xc, n, k);flag=1;
    while(flag==1){
      firstcomb(yc, n, k);flag1=1;
      while(flag1==1){
	printvec1(xc,k);printvec1(yc,k);
	detex=symbdetsubmat(J, n, n, xc, xc, k);
	ndetex=-detex;
	if(detex!=0 && !(detex.info(info_flags::positive))&& !(ndetex.info(info_flags::positive))){
	  cout << detex << endl;
	  printexsubmat(J,xc,xc,k,k);
	  cout << "An apparently unsigned principal minor found.\n";
	  outf=0;
	}
	else if(detex==0)
	  cout << "0" << endl;
	else if(detex.info(info_flags::positive))
	  cout << "+" << endl;
	else if(detex.info(info_flags::negative))
	  cout << "-" << endl;

	flag1=nextcomb(yc, n, k);
      }
      flag=nextcomb(xc, n, k);
    }
    free((char *) xc);
  }
  return outf;
  
}


matrix multexAB(matrix S, matrix V, int n, int m){
  matrix J(n,n);
  int i,j,k;
  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      J(i,j)=0;
      for(k=0;k<m;k++)
	J(i,j)+=S(i,k)*V(k,j);
    }
  }
  return J;
}




int isQmatrix(matrix J, int n){
  int *xc;
  int k;
  int flag;
  ex detex;
  for(k=1;k<=n;k++){
    xc=(int *)malloc((size_t) (k*sizeof(int)));
    firstcomb(xc, n, k);flag=1;
    detex=0;
    while(flag==1){ // sum all minors of particular size
      //    printvec1(xc,k);
      detex+=symbdetsubmat(J, n, n, xc, xc, k);
      printvec1(xc,k);
      flag=nextcomb(xc, n, k);
    }
    if(detex!=0 && !(detex.info(info_flags::positive))){
      cout << detex << endl;
      cout << "An apparently negative coefficient found.\n";
      return 0;
    }
    free((char *) xc);
  }
  return 1;
}

// The number of different variables in a monomial of degree n

int width(int *k, int n){
  int j, r=1;
  for(j=1;j<n;j++)
    if(k[j]!=k[j-1])
      r++;
  return r;
}

void monovtomonoexp(int *k, int r, int *out, int numv){
  int j;
  for(j=0;j<numv;j++)
    out[j]=0;
  for(j=0;j<r;j++)
    (out[k[j]-1])++;
  return;
}

// inhomogeneous version; assumes the first entry in k is the order

void monovtomonoexp0(int *k, int *out, int numv){
  int j;
  //printvec1(k,k[0]+1);
  for(j=0;j<numv;j++)
    out[j]=0;
  for(j=1;j<k[0]+1;j++)
    (out[k[j]-1])++;
  //printvec1(out,numv);
  return;
}

//
// swtch causes only negative terms to be printed out
//

void treeprint3(int *cfs, int **mons, int dim, long len, int numv, bool swtch){
  long k;
  bool flg;
  int *exps=(int *)malloc((size_t) (numv*sizeof(int)));
  for(k=0;k<len;k++){
    flg=1;
    if(swtch && cfs[k]>=0) // only negative terms
      flg=0;
    if(flg){
      cout << cfs[k] << ", ";
      //printvec1(mons[k],dim);
      monovtomonoexp(mons[k],dim,exps,numv);
      printvec1(exps,numv);
    }
  }
  free((char*)exps);
}

int **monlisttoexplist(int **mons, int dim, long len, int numv, bool srt){
  long k;
  int **expsout=(int **)malloc((size_t) (len*sizeof(int*)));
  for(k=0;k<len;k++){
    expsout[k]=(int *)malloc((size_t) (numv*sizeof(int)));
    monovtomonoexp(mons[k],dim,expsout[k],numv);
  }

  if(srt)
    qsorti(expsout,numv,0,len-1);
  return expsout;
}

int **moncflisttoexplist(int **mons, int dim, int *cfs, long len, int numv, bool srt){
  long k;
  int **expsout=(int **)malloc((size_t) (len*sizeof(int*)));
  for(k=0;k<len;k++){
    expsout[k]=(int *)malloc((size_t) (numv*sizeof(int)));
    monovtomonoexp(mons[k],dim,expsout[k],numv);
  }

  if(srt)
    qsorti2(expsout,cfs,numv,0,len-1);
  return expsout;
}

// don't assume homogeneous polynomials

int **moncflisttoexplist0(int **mons, int *cfs, long len, int numv, bool srt){
  long k;
  int **expsout=(int **)malloc((size_t) (len*sizeof(int*)));
  for(k=0;k<len;k++){
    expsout[k]=(int *)malloc((size_t) (numv*sizeof(int)));
    monovtomonoexp0(mons[k],expsout[k],numv);
  }

  if(srt)
    qsorti2(expsout,cfs,numv,0,len-1);
  return expsout;
}

//
// change all coefficients to have magnitude 1 at the cost of repetition
//

void elongate(int ***explst, int **cflst, long *r, int numv){
  long s,s1=0;
  int tmp[(*r)][numv], cftmp[(*r)];
  int i,j,k,t;
  for(s=0;s<(*r);s++){
    for(i=0;i<numv;i++)
      tmp[s][i]=(*explst)[s][i];
    cftmp[s]=(*cflst)[s];
  }
  ifree_l((*explst),(*r));
  free((char*)(*cflst));

  for(s=0;s<(*r);s++)
    s1+=abs(cftmp[s]);

  (*explst)= (int**) malloc(sizeof(int*)*s1);
  (*cflst)= (int*) malloc(sizeof(int)*s1);
  for(s=0;s<s1;s++)
    (*explst)[s]= (int*) malloc(sizeof(int)*numv);

  s1=0;
  for(s=0;s<(*r);s++){
    t=abs(cftmp[s]);
    for(j=0;j<t;j++){
      for(k=0;k<numv;k++){
	(*explst)[s1][k]=tmp[s][k];
	(*cflst)[s1]=cftmp[s]/t;
      }
      s1++;
    }
  }
  (*r)=s1;
  return;
}

// zeros of vec1 are also zeros of vec2
// i.e. vec2 lies in the same closed face as vec1
// (not a symmetric relationship)

bool sameface(int *vec1, int *vec2, int l){
  int i;
  for(i=0;i<l;i++){
    if(vec1[i]==0 && vec2[i]!=0)
      return 0;
  }
  return 1;
}

long allsameface(int *vec1, int **veclist, long listlen, int veclen, long **vecout){
  long r;
  long tot=0;
  for(r=0;r<listlen;r++){
    if(sameface(vec1, veclist[r], veclen)){
      if(!tot)
	(*vecout)=(long*)malloc((size_t) (sizeof(long)));
      else
	(*vecout)=(long*)realloc((*vecout), sizeof(long)*(tot+1));
      (*vecout)[tot++]=r;
    }

  }
  return tot;
}



int numsameface(int *vec1, int **veclist, int *cflist, long listlen, int veclen, int *totneg){
  long r;
  long tot=0;
  *totneg=0;
  for(r=0;r<listlen;r++){
    if(sameface(vec1, veclist[r], veclen)){
      tot++;
      if(cflist[r]<0)
	(*totneg)++;
    }
  }
  return tot;
}

long allsameface1(int *vec1, int **veclist, int *cflist, long listlen, int veclen, long **vecout, int *totneg){
  long r;
  long tot=0;
  long j=0;
  tot=numsameface(vec1, veclist, cflist, listlen, veclen, totneg);
  (*vecout)=(long*)malloc((size_t) (sizeof(long))*tot);
  for(r=0;r<listlen;r++){//first negative
    if(sameface(vec1, veclist[r], veclen)){
      if(cflist[r]<0)
	(*vecout)[j++]=r;
    }
  }
  for(r=0;r<listlen;r++){//first negative
    if(sameface(vec1, veclist[r], veclen)){
      if(cflist[r]>0)
	(*vecout)[j++]=r;
    }
  }
  return tot;
}

ex getminorsum(matrix M, int n, int k){
  int *xc, *yc;
  int j,flag,flag1;
  ex out=0;
  if(k==0){
    out=1;
    return out;
  }
  xc=(int *)malloc((size_t) (k*sizeof(int)));
  yc=(int *)malloc((size_t) (k*sizeof(int)));
  firstcomb(xc,n,k);

  flag1=1;
  while(flag1==1){//diagonal terms
    out+=expand(pow(symbdetsubmat(M,n,n,xc,xc,k),2));
    flag1=nextcomb(xc,n,k);
  }

  firstcomb(xc,n,k);
  firstcomb(yc,n,k);
  flag1=1;
  while(flag1==1){//diagonal terms
    flag=1;
    for(j=0;j<k;j++)
      yc[j]=xc[j];
    flag=nextcomb(yc,n,k);
    while(flag==1){
      out+=expand(2*symbdetsubmat(M,n,n,xc,yc,k)*symbdetsubmat(M,n,n,yc,xc,k));
      flag=nextcomb(yc,n,k);
    }
    flag1=nextcomb(xc,n,k);
  }
  free((char*)xc);free((char*)yc);
  return out;      
}

// get the sum of principal minors of dimension k

ex getminorsum0(matrix M, int n, int k){
  int *xc;
  int flag1;
  ex out=0;
  if(!k){
    out=1;
    return out;
  }
  xc=(int *)malloc((size_t) (k*sizeof(int)));
  firstcomb(xc,n,k);

  flag1=1;
  while(flag1==1){//diagonal terms
    out+=expand(symbdetsubmat(M,n,n,xc,xc,k));
    flag1=nextcomb(xc,n,k);
  }
  return out;      
}



void ifree(int **imat, int n){
  int j;
  if(n>0){
    for(j=0;j<n;j++)
      free ((char *)(imat[j]));
    free((char *)imat);
  }
}

void ifree_l(int **imat, long n){
  long j;
  if(n>0){
    for(j=0;j<n;j++)
      free ((char *)(imat[j]));
    free((char *)imat);
  }
}

void lfree_i(long **imat, int n){
  int j;
  if(n>0){
    for(j=0;j<n;j++)
      free ((char *)(imat[j]));
    free((char *)imat);
  }
}

void lkeeponly(long **imat, int n, int m){
  int j;
  if(n>0){
    for(j=m;j<n;j++)
      free ((char *)(imat[j]));
  }
}

bool sum0(int *v1, int *v2, int *v3, int k){
  int i;
  for(i=0;i<k;i++){
    if(2*v3[i]!=v1[i]+v2[i])
      return 0;
  }
  return 1;
}

// a perfect square
// working here...

bool psq(int **lst, int *cflst, int i, int j, int k, int numv){
  if(cflst[i]>=0 &&cflst[j]>=0 && cflst[k]>=0)
    return 0;
  if((cflst[i]<0 && cflst[j]<0) || (cflst[i]<0 && cflst[k]<0) || (cflst[j]<0 && cflst[k]<0))
    return 0;
  if(4*cflst[i]*cflst[j]==pow(cflst[k],2) && sum0(lst[i],lst[j],lst[k],numv))
    return 1;
  if(4*cflst[j]*cflst[k]==pow(cflst[i],2) && sum0(lst[j],lst[k],lst[i],numv))
    return 1;
  if(4*cflst[k]*cflst[i]==pow(cflst[j],2) && sum0(lst[k],lst[i],lst[j],numv))
    return 1;
  return 0;
}

// remove perfect squares

void cancel(int **lst, int *cflst, int num, int **inds, int numv){
  int i,j,k;
  for(i=0;i<num;i++){
    for(j=i+1;j<num;j++){
      for(k=j+1;k<num;k++){
	if((*inds)[i]&&(*inds)[j]&&(*inds)[k]&&psq(lst,cflst,i,j,k,numv)){
	  (*inds)[i]=0;(*inds)[j]=0;(*inds)[k]=0;
	}
      }
    }
  }
}

void addnewtoarray(int ***outmat, int sz, int *xc, int len){
  int l;
  if(sz==0) // first entry
    (*outmat) = (int**) malloc(sizeof(int*) * 1);
  else
    (*outmat)=(int**) realloc((*outmat), sizeof(int*) *(sz+1));
  (*outmat)[sz]=(int *)malloc(sizeof(int) *(len));
  for(l=0;l<len;l++)
    (*outmat)[sz][l]=xc[l];
  return;
}

void addnewtolarray(long ***outmat, int sz, long *xc, int len){
  int l;
  if(sz==0) // first entry
    (*outmat) = (long**) malloc(sizeof(long*) * 1);
  else
    (*outmat)=(long**) realloc((*outmat), sizeof(long*) *(sz+1));
  (*outmat)[sz]=(long *)malloc(sizeof(long) *(len));
  for(l=0;l<len;l++)
    (*outmat)[sz][l]=xc[l];
  return;
}

void addnewto1Darray(int **outvec, int sz, int newint){
  if(sz==0) // first entry
    (*outvec) = (int*) malloc(sizeof(int) * 1);
  else
    (*outvec)=(int*) realloc((*outvec), sizeof(int) *(sz+1));
  (*outvec)[sz]=newint;
  return;
}

void addnewto1Dlarray(long **outvec, int sz, long newint){
  if(sz==0) // first entry
    (*outvec) = (long *) malloc(sizeof(long) * 1);
  else
    (*outvec)=(long *) realloc((*outvec), sizeof(long) *(sz+1));
  (*outvec)[sz]=newint;
  return;
}

bool isinintarray(int *v, int numv, int s){
  // is the integer s a member of the array v?
  int i;
  for(i=0;i<numv;i++)
    if(s==v[i])
      return 1;
  return 0;
}

int unionvec(int *a, int na, int *b, int nb, int **outvec){
  int i,k=0;
  for(i=0;i<na;i++){
    addnewto1Darray(outvec,k,a[i]);
    k++;
  }
  for(i=0;i<nb;i++){
    if(!isinintarray(*outvec,k, b[i])){
      addnewto1Darray(outvec,k,b[i]);
      k++;
    }
  }
  return k;
}

int unionlvec(long *a, int na, long *b, int nb, long **outvec){
  int i;
  (*outvec) = (long*) malloc(sizeof(long) * (na+nb));
  for(i=0;i<na;i++)
    (*outvec)[i]=a[i];
  for(i=0;i<nb;i++)
    (*outvec)[na+i]=b[i];
    return na+nb;
}

ex intltopoly_negpart(int **l, int *cfs, long numinds, int n){
  char str1[10];
  ex extot=0;
  ex ex1=1;
  long i;
  int j;
  
  for(i=0;i<numinds;i++){//each monomial
    if(cfs[i]<0){
      ex1=cfs[i];
      for(j=0;j<n;j++){ // each term
	sprintf(str1, "v%d", l[i][j]);
	ex1*=get_possymbol(str1);
      }
      extot+=ex1;
    }
  }
  return extot;
}

// sz is the number of blocks so far
// blocksz is the size of the new block

void addnewblock(int ****outmat, int sz, int **newblock, int blocksz, int len){
  int l;
  if(sz==0) // first entry
    (*outmat) = (int***) malloc(sizeof(int**) * 1);
  else
    (*outmat)=(int***) realloc((*outmat), sizeof(int**) *(sz+1));
  for(l=0;l<blocksz;l++)
    addnewtoarray((*outmat)+sz,l,newblock[l],len);
  return;
}

void freeblocklist(int ***imat, int *blocksz, int n){
  int i,j;
  if(n>0){
    for(j=0;j<n;j++){//each block
      for(i=0;i<blocksz[j];i++)// each vector in block
	free ((char *)(imat[j][i]));
      free((char *)(imat[j]));
    }
    free((char *)imat);
    free((char *)blocksz);//free the array of block sizes
  }
}

void freeblockcflist(int ***imat, int **cf, int *blocksz, int n){
  int i,j;
  if(n>0){
    for(j=0;j<n;j++){//each block
      for(i=0;i<blocksz[j];i++)// each vector in block
	free ((char *)(imat[j][i]));
      free((char *)(imat[j]));
      free((char *)(cf[j]));
    }
    free((char *)imat);free((char *)cf);
    free((char *)blocksz);//free the array of block sizes
  }
}

//
// Take the negative terms only and make an (exponent) list
//

int **monlisttonegexplist(int **mons, int dim, int *cfs, long len, int numv, int *blocksz){
  long k;
  int **expsout=NULL;
  (*blocksz)=0;
  for(k=0;k<len;k++){
    if(cfs[k]<0){
      if((*blocksz)==0) // first entry
	expsout = (int**) malloc(sizeof(int*) * 1);
      else
	expsout=(int**) realloc(expsout, sizeof(int*) *((*blocksz)+1));
      expsout[(*blocksz)]=(int *)malloc((size_t) (numv*sizeof(int)));
      monovtomonoexp(mons[k],dim,expsout[(*blocksz)],numv);
      (*blocksz)++;
    }
  }
  return expsout;
}


void monovtomonoexphom(int *k, int r, int *out, int numv, int deg){
  int j;
  int dg=0;// the degree
  for(j=0;j<numv;j++)
    out[j]=0;
  for(j=0;j<r;j++){
    (out[k[j]-1])++;
    dg++;
  }
  out[numv]=deg-dg;//to make it homogeneous of degree deg
  return;
}
//
// Take the negative terms only and make an (exponent) list
// Also make homogeneous of degree deg
//

int **monlisttohomnegexplist(int **mons, int dim, int *cfs, long len, int numv, int *blocksz, int deg){
  long k;
  int **expsout=NULL;
  (*blocksz)=0;
  for(k=0;k<len;k++){
    if(cfs[k]<0){
      if((*blocksz)==0) // first entry
	expsout = (int**) malloc(sizeof(int*) * 1);
      else
	expsout=(int**) realloc(expsout, sizeof(int*) *((*blocksz)+1));
      expsout[(*blocksz)]=(int *)malloc((size_t) ((numv+1)*sizeof(int)));
      monovtomonoexphom(mons[k],dim,expsout[(*blocksz)],numv,deg);
      (*blocksz)++;
    }
  }
  return expsout;
}

void monoprint(int cf, int *mon, int dim){
  int i;
  bool st=0;
  if(cf==0){
    cout << "0\n";
    return;
  }
  if(cf!=1){
    cout << cf;
    st=1;
  }
  for(i=0;i<dim;i++){
    if(mon[i]){
      if(st){cout << "*";}else{st=1;}
      if(mon[i]>1)
	cout << "v" << i+1 << "^" << mon[i];
      else
	cout << "v" << i+1;
    }
  }
  cout << endl;
  return;
}

//a-b

void vminus(int *a, int *b, int **out, int n){
  int i;
  for(i=0;i<n;i++)
    (*out)[i]=a[i]-b[i];
}

//a+b

void vplus(int *a, int *b, int **out, int n){
  int i;
  for(i=0;i<n;i++)
    (*out)[i]=a[i]+b[i];
}

//2*a-b

void vshift(int *a, int *b, int **out, int n){
  int i;
  for(i=0;i<n;i++)
    (*out)[i]=2*a[i]-b[i];
}

bool eq(int *a, int *b, int n){
  int i;
  for(i=0;i<n;i++){
    if(a[i]!=b[i])
      return 0;
  }
  return 1;
}

bool factorpos(ex ex1){
  ex tmp1;
  unsigned int j1;
  tmp1=factor(ex1);
  // even power
  if(is_a<power>(tmp1) && ((ex_to<numeric>(tmp1.op(1))).to_int())%2==0)
    return 1;
  if(is_a<mul>(tmp1)){
    for (j1=0; j1!=tmp1.nops(); ++j1){
      if(is_a<numeric>(tmp1.op(j1)) && (ex_to<numeric>(tmp1.op(j1))>=0))
	continue;
      if(is_a<symbol>(tmp1.op(j1))) //assume symbols are positive
	continue;
      if(is_a<power>(tmp1.op(j1))){
	if(is_a<symbol>(tmp1.op(j1).op(0)) || ((ex_to<numeric>(tmp1.op(j1).op(1))).to_int())%2==0)
	  continue;
      }
      return 0;
    }
  }
  return 1;
}



void checkfactor(){
  possymbol w("w"),x("x"),z("z"),y("y");
  ex tmp=expand(pow(z,3)*pow((pow(z,2)*(w+x)-w*x*y),2));
  if(factorpos(tmp))
    cout << "nonnegative term\n";
  return;
}



ex intltomono(int cf, int *l, int n){
  char str1[10];
  ex ex1=1;
  int i;
  ex1*=cf;
  for(i=0;i<n;i++){
    sprintf(str1, "v%d", l[i]);
    ex1*=get_possymbol(str1);
  }
  return ex1;
}

ex expltomono(int cf, int *l, int n){
  char str1[10];
  ex ex1=1;
  int i;
  ex1*=cf;
  for(i=0;i<n;i++){
    if(l[i]!=0){
      sprintf(str1, "v%d", i+1);
      ex1*=pow(get_possymbol(str1),l[i]);
    }
  }
  return ex1;
}

void expltomonoprint(int cf, int *l, int n){
  int i;
  if(cf!=1)
    printf("%d ", cf);
  for(i=0;i<n;i++){
    if(l[i]!=0){
      if(l[i]==1)
	printf("v%d ", i+1);
      else
	printf("v%d^%d ", i+1, l[i]);
    }
  }
  printf("\n");
}

// inefficient perhaps?

ex intltopoly(int **l, int *cfs, long *inds, long numinds, int n){
  char str1[10];
  ex extot=0;
  ex ex1=1;
  long i;
  int j;
  
  for(i=0;i<numinds;i++){//each monomial
    ex1=cfs[inds[i]];
    for(j=0;j<n;j++){ // each term
      sprintf(str1, "v%d", l[inds[i]][j]);
      ex1*=get_possymbol(str1);
    }
    extot+=ex1;
  }
  return extot;
}

ex expltopoly(int **l, int *cfs, long *inds, long numinds, int n){
  char str1[10];
  ex extot=0;
  ex ex1=1;
  long i;
  int j;
  for(i=0;i<numinds;i++){//each monomial
    ex1=cfs[inds[i]];
    for(j=0;j<n;j++){ // each symbol
      if(l[inds[i]][j]){
	sprintf(str1, "v%d", j+1);
	ex1*=pow(get_possymbol(str1),l[inds[i]][j]);
      }
    }
    extot+=ex1;
  }
  return extot;
}

ex expltopolysimp(int **l, int *cfs, long numinds, int n){
  char str1[10];
  ex extot=0;
  ex ex1=1;
  long i;
  int j;
  for(i=0;i<numinds;i++){//each monomial
    if(cfs)
      ex1=cfs[i];
    else
      ex1=1;
    for(j=0;j<n;j++){ // each symbol
      if(l[i][j]){
	sprintf(str1, "v%d", j+1);
	ex1*=pow(get_possymbol(str1),l[i][j]);
      }
    }
    extot+=ex1;
  }
  return extot;
}


ex lrexptopolysimp(int **l, int num1, int *cfs, long numinds, int n){
  char str1[10];
  ex extot=0;
  ex ex1=1;
  long i;
  int j;
  for(i=0;i<numinds;i++){//each monomial
    ex1=cfs[i];
    for(j=0;j<n;j++){ // each symbol
      if(l[i][j]){
	sprintf(str1, "v%d", j+1);
	ex1*=pow(get_possymbol(str1),l[i][j]);
      }
    }
    extot+=ex1;
  }
  return extot;
}

//    2*sqrt(cfs[s])*sqrt(cfs[s+s1]);
int splt(int a, int b){
  if(a<0 || b<0){
    fprintf(stderr,"ERROR in splt a=%d,b=%d. Exiting.\n",a,b);
    exit(0);
  }
  if(a==0 || b==0)
    return 0;
  int i=sqrt(a);// integer part of sqrt
  int j=sqrt(b);// integer part of sqrt
  return 2*i*j*min(a/(i*i),b/(j*j));
}



void trysquare(){
  cout << "1,1: " << splt(1,1) << endl;
  cout << "1,2: " << splt(1,2) << endl;
  cout << "2,1: " << splt(2,1) << endl;
  cout << "2,2: " << splt(2,2) << endl;
  cout << "3,2: " << splt(3,2) << endl;
  cout << "1,4: " << splt(1,4) << endl;
  cout << "2,8: " << splt(2,8) << endl;

  /* for(i=0;i<100;i++){ */
  /*   j=sqrt(i); */
  /*   cout << i << ", " << j*j << endl; */
  /* } */
  exit(0);

}

int *midpt(int *vecl, int *vecr, int cfl, int cfr, int n){
  int *vecmid;
  int tmp;
  int i;
  if(cfl<0 || cfr<0) // can't dominate negative with negative
    return NULL;
  vecmid=(int *)malloc((size_t) (n*sizeof(int)));
  for(i=0;i<n;i++){
    tmp=vecl[i]+vecr[i];
    if(tmp%2!=0){// only even allowed
      free((char *)vecmid);
      return NULL;
    }
    else
      vecmid[i]=tmp/2;
  }
  return vecmid;
}

int *intersect(int **l, long *leftinds, int numleft, long *rightinds, int numright, int numv){// intersection of a set of monomials
  int *isect=(int *)malloc((size_t) (numv*sizeof(int)));
  int i,j;
  if(numleft+numright<=0){
    for(i=0;i<numv;i++)
      isect[i]=0;
  }
  for(i=0;i<numv;i++)
    isect[i]=100;

  for(i=0;i<numv;i++){
    for(j=0;j<numleft;j++){
      isect[i]=min(isect[i],l[leftinds[j]][i]);
      if(!isect[i])
	break;
    }
    if(isect[i]){
      for(j=0;j<numright;j++){
	isect[i]=min(isect[i],l[rightinds[j]][i]);
	if(!isect[i])
	  break;
      }
    }
  }
  return isect;
}

int *subtract(int *a, int *b, int numv){
  int i;
  int *out=(int *)malloc((size_t) (numv*sizeof(int)));
  for(i=0;i<numv;i++)
    out[i]=a[i]-b[i];
  return out;
}

int *halve(int *a, int n){
  int i;
  int *out=(int *)malloc((size_t) (n*sizeof(int)));
  for(i=0;i<n;i++)
    out[i]=a[i]/2;
  return out;
}

// Vector vec of length n has only even entries

bool iseven(int *vec, int n){
  int i;
  for(i=0;i<n;i++){
    if(vec[i]%2!=0)
      return 0;
  }
  return 1;
}

ex facsq(int **l, long *leftinds, int *kleft, int numleft, long *rightinds, int *kright, int numright, int numv, int minfac){
  int *exfaci;
  int i,*vc,*vc1;
  ex exout=0,exfac;
  exfaci=intersect(l, leftinds, numleft, rightinds, numright, numv);
  for(i=0;i<numleft;i++){
    vc=subtract(l[leftinds[i]],exfaci,numv);
    if(!iseven(vc,numv)){
      exout=-1;free((char*)vc);return -1;
    }
    vc1=halve(vc,numv);
    /* cout << "hh: " << expltomono(kleft[i],vc1,numv) << endl; */
    exout+=expltomono(kleft[i],vc1,numv);
    free((char*)vc1);free((char*)vc);
  }

  for(i=0;i<numright;i++){
    vc=subtract(l[rightinds[i]],exfaci,numv);
    if(!iseven(vc,numv)){
      exout=-1;free((char*)vc);return -1;
    }
    vc1=halve(vc,numv);
    /* cout << "hh: " << expltomono(kright[i],vc1,numv) << endl; */
    exout-=expltomono(kright[i],vc1,numv);
    free((char*)vc1);free((char*)vc);
  }
  exfac=expltomono(1,exfaci,numv);
  free((char*)exfaci);
  return minfac*exfac*pow(exout,2);
}


// make a perfect square from a list
// takes the terms on the left and those on the right
// computes the midpoints
// assumes that consistency checks have been done earlier

ex makesquare(int **l, int *cfs, long *leftinds, int numleft, long *rightinds, int numright, int n, int maxnegcf, bool simp){
  ex extot=0,extot1=0;
  ex ex1=1;
  int i,j,k,k2;
  int minfac=1000;
  int *vc;
  int *kleft=(int *)malloc((size_t) (numleft*sizeof(int)));
  int *kright=(int *)malloc((size_t) (numright*sizeof(int)));

  if(maxnegcf==2){// simplest case
    minfac=1;
    for(i=0;i<numleft;i++)
      kleft[i]=1;
    for(i=0;i<numright;i++)
      kright[i]=1;
  }
  else if(simp){
    minfac=1;
    for(i=0;i<numleft;i++)
      kleft[i]=sqrt(cfs[leftinds[i]]);
    for(i=0;i<numright;i++)
      kright[i]=sqrt(cfs[rightinds[i]]);//integer square root
  }
  else{
    // get the smallest common factor
    for(i=0;i<numleft;i++){
      k=sqrt(cfs[leftinds[i]]);//integer square root
      kleft[i]=k;
      minfac=min(minfac,cfs[leftinds[i]]/(k*k));//integer division
    }
    for(i=0;i<numright;i++){
      k=sqrt(cfs[rightinds[i]]);//integer square root
      kright[i]=k;
      minfac=min(minfac,cfs[rightinds[i]]/(k*k));//integer division
    }
  }
  extot1=facsq(l, leftinds, kleft, numleft, rightinds, kright, numright, n, minfac);
  if(extot1!=-1){
    free((char*)kleft);free((char*)kright);
    return extot1;
  }


  for(i=0;i<numleft;i++){//each monomial on left
    ex1=expltomono(minfac*kleft[i]*kleft[i], l[leftinds[i]],n);
    extot+=ex1;
  }
  for(i=0;i<numright;i++){//each monomial on right
    ex1=expltomono(minfac*kright[i]*kright[i], l[rightinds[i]],n);
    extot+=ex1;
  }

  // left midpoints (positive)

  for(i=0;i<numleft;i++){
    for(j=i+1;j<numleft;j++){
      vc=midpt(l[leftinds[i]],l[leftinds[j]],cfs[leftinds[i]], cfs[leftinds[j]],n);
      k2=2*kleft[i]*kleft[j]*minfac;
      ex1=expltomono(k2, vc,n);
      if(vc)
	free((char*)vc);
      extot+=ex1;
    }
  }

  // right midpoints (positive)

  for(i=0;i<numright;i++){
    for(j=i+1;j<numright;j++){
      vc=midpt(l[rightinds[i]],l[rightinds[j]],cfs[rightinds[i]], cfs[rightinds[j]],n);
      k2=2*kright[i]*kright[j]*minfac;
      ex1=expltomono(k2, vc,n);
      if(vc)
	free((char*)vc);
      extot+=ex1;
    }
  }

 // mixed midpoints (negative)

  for(i=0;i<numleft;i++){
    for(j=0;j<numright;j++){
      vc=midpt(l[leftinds[i]],l[rightinds[j]],cfs[leftinds[i]], cfs[rightinds[j]],n);
      k2=2*kleft[i]*kright[j]*minfac;
      ex1=expltomono(k2, vc,n);
      if(vc)
	free((char*)vc);
      extot-=ex1;
    }
  }

 


  /* if(extot-expand(extot1)!=0){ */
  /*   fprintf(stderr, "serious problem!\n"); */
  /*   cout << extot << endl; */
  /*   cout << extot1 << endl; */
  /*   exit(0); */
  /* } */
  free((char*)kleft);free((char*)kright);
  /* if(extot1==-1) */
  /*   return extot; */
  return extot;
}




// make a perfect square from a list
// takes the terms on the left and those on the right
// computes the midpoints
// assumes that consistency checks have been done earlier

ex makesquare1(int **l, int *cfs, long r, long *leftinds, int numleft, long *rightinds, int numright, int n){
  ex extot=0;
  ex ex1=1;
  long s1;
  int i,j,k,k2,m;
  int minfac=1000;
  int maxfac=1;
  int *vc;
  int *kleft=(int *)malloc((size_t) (numleft*sizeof(int)));
  int *kright=(int *)malloc((size_t) (numright*sizeof(int)));

  // get the smallest common factor

  for(i=0;i<numleft;i++){
    k=sqrt(cfs[leftinds[i]]);//integer square root
    kleft[i]=k;
    minfac=min(minfac,cfs[leftinds[i]]/(k*k));//integer division
  }
  for(i=0;i<numright;i++){
    k=sqrt(cfs[rightinds[i]]);//integer square root
    kright[i]=k;
    minfac=min(minfac,cfs[rightinds[i]]/(k*k));//integer division
  }

  // check mixed midpoints

  for(i=0;i<numleft;i++){
    for(j=0;j<numright;j++){
      vc=midpt(l[leftinds[i]],l[rightinds[j]],cfs[leftinds[i]], cfs[rightinds[j]],n);
      s1=binisearch(l,r,vc,n);
      if(cfs[s1]<0){
	m=-cfs[s1]/(2*kleft[i]*kright[j]);
	//cout << "  m: " << m << ", cfs[s1]: " << cfs[s1] << endl;
	maxfac=max(maxfac,m);
      }
      if(vc)
	free((char*)vc);
    }
  }
  if(maxfac>minfac)
    cout << "warning: maxfac must be less than minfac\n"<< "  maxfac: " << maxfac << ", minfac: " << minfac << endl;

  for(i=0;i<numleft;i++){//each monomial on left
    ex1=expltomono(maxfac*kleft[i]*kleft[i], l[leftinds[i]],n);
    extot+=ex1;
  }
  for(i=0;i<numright;i++){//each monomial on right
    ex1=expltomono(maxfac*kright[i]*kright[i], l[rightinds[i]],n);
    extot+=ex1;
  }

  // left midpoints (positive)

  for(i=0;i<numleft;i++){
    for(j=i+1;j<numleft;j++){
      vc=midpt(l[leftinds[i]],l[leftinds[j]],cfs[leftinds[i]], cfs[leftinds[j]],n);
      k2=2*kleft[i]*kleft[j]*maxfac;
      ex1=expltomono(k2, vc,n);
      if(vc)
	free((char*)vc);
      extot+=ex1;
    }
  }

  // right midpoints (positive)

  for(i=0;i<numright;i++){
    for(j=i+1;j<numright;j++){
      vc=midpt(l[rightinds[i]],l[rightinds[j]],cfs[rightinds[i]], cfs[rightinds[j]],n);
      k2=2*kright[i]*kright[j]*maxfac;
      ex1=expltomono(k2, vc,n);
      if(vc)
	free((char*)vc);
      extot+=ex1;
    }
  }

 // mixed midpoints (negative)

  for(i=0;i<numleft;i++){
    for(j=0;j<numright;j++){
      vc=midpt(l[leftinds[i]],l[rightinds[j]],cfs[leftinds[i]], cfs[rightinds[j]],n);
      k2=2*kleft[i]*kright[j]*maxfac;
      ex1=expltomono(k2, vc,n);
      if(vc)
	free((char*)vc);
      extot-=ex1;
    }
  }

  free((char*)kleft);free((char*)kright);

  return extot;
}





ex intpairtopoly(int *l1,int cf1,int n1, int *l2, int cf2,int n2){
  char str1[10];
  ex extot=0;
  ex ex1=1;
  int j;
  
  ex1=cf1;
  for(j=0;j<n1;j++){
    sprintf(str1, "v%d", l1[j]);
    ex1*=get_possymbol(str1);
  }
  extot+=ex1;
  ex1=cf2;
  for(j=0;j<n2;j++){ // each term
    sprintf(str1, "v%d", l2[j]);
    ex1*=get_possymbol(str1);
  }
  extot+=ex1;
  return extot;
}

// convert a sublist of l of length m indexed by the entries in t into a polynomial. 

ex intltopoly1(int **l, int *cfs, long *inds, long numinds, int n, int *t, int m){
  char str1[10];
  ex extot=0;
  ex ex1=1;
  long i;
  int j;
  
  for(i=0;i<m;i++){//each monomial
    ex1=cfs[inds[t[i]]];
    for(j=0;j<n;j++){ // each term
      sprintf(str1, "v%d", l[inds[t[i]]][j]);
      ex1*=get_possymbol(str1);
    }
    extot+=ex1;
  }
  return extot;
}

matrix smalltry(int *n, int *numv){
  (*n)=3;(*numv)=3;
  int i;
  ex v[(*numv)+1];
  char str1[10];
  matrix J((*n),(*n));
  for(i=0;i<=(*numv);i++){
    sprintf(str1, "v%d", i);
    v[i]=get_possymbol(str1);
  }
  J(0,0)=v[1];J(0,1)=v[2];J(0,2)=v[3];
  J(1,0)=-v[1];J(1,1)=v[2];J(1,2)=0;
  J(2,0)=-v[1];J(2,1)=0;J(2,2)=v[3];
 
  return J;
}

matrix sixbysix(){
  int n=6;int numv=17;
  int i;
  ex v[numv+1];
  char str1[10];
  matrix J(n,n);
  for(i=0;i<=numv;i++){
    sprintf(str1, "v%d", i);
    v[i]=get_possymbol(str1);
  }
  J(0,0)=v[1]+v[2];J(0,1)=v[3];J(0,2)=-v[8]-v[9];
  J(1,0)=v[1];J(1,1)=v[3];J(1,2)=0;J(1,3)=0;J(1,4)=0;J(1,5)=-v[16];
  J(2,0)=-v[2];J(2,1)=0;J(2,2)=v[8]+v[9];
  J(3,2)=-v[9];J(3,3)=v[10];J(3,4)=v[11];
  J(4,3)=v[10];J(4,4)=v[11]+v[12];J(4,5)=-v[15]-v[16];
  J(5,4)=-v[12];J(5,5)=v[15]+v[16];
  return J;
}



matrix sixbysix1(){
  int n=6;int numv=17;
  int i;
  ex v[numv+1];
  char str1[10];
  matrix J(n,n);
  for(i=0;i<=numv;i++){
    sprintf(str1, "v%d", i);
    v[i]=get_possymbol(str1);
  }
  J(0,0)=v[2];J(0,1)=v[3];J(0,2)=-v[9];
  J(1,0)=0;J(1,1)=v[3];J(1,5)=-v[16];
  J(2,0)=-v[2];J(2,1)=0;J(2,2)=v[9];
  J(3,2)=-v[9];J(3,3)=v[10];J(3,4)=0;
  J(4,3)=v[10];J(4,4)=v[12];J(4,5)=-v[16];
  J(5,4)=-v[12];J(5,5)=v[16];
  return J;
}



//
// Compute the determinant of J^2 + omega^2 I
// from matrix outputting function J
//

void J2pI_f(matrix (*Jcomp)(int*,int*)){
  //, int n,int numv){
  int n,numv;
  matrix J=(*Jcomp)(&n,&numv);
  ex v[numv+1];
  char str1[10];
  int i,k,dim;
  ex tmp,tmp1;
  int **lst;int *cflst;
  long r;

  //  matrix J(n,n);

  for(i=0;i<=numv+1;i++){
    sprintf(str1, "v%d", i);
    v[i]=get_possymbol(str1);
  }


  for(k=0;k<=n;k++){
    tmp1=expand(getminorsum(J,n,k));
    if(tmp1!=0 && !(tmp1.info(info_flags::positive)))
      fprintf(stderr, "could be negative: k=%d\n", k);
    tmp+=expand(pow(v[numv+1], 2*(n-k))*tmp1);
  }
  dim=polytointlist(tmp, &lst, &cflst, &r);//dim=order or monomials
  treeprint3(cflst,lst,dim,r,numv+1,0);

  for(i=0;i<r;i++)
    free ((char *)(lst[i]));
  free((char *)lst);
  free((char *)cflst);
  return;
}

bool ueq2v(int *u, int *v, int n){
  int i;
  for(i=0;i<n;i++){
    if(2*v[i]!=u[i])
      return 0;
  }
  return 1;
}

// vecl --- vecmid --- vecr

int *sympt(int *vecmid, int *vecl, int cf, int n, bool pqswitch){
  int *vecr;
  int i;
  if(pqswitch){
    if(!iseven(vecl,n))
      return NULL;
  }
  if(cf<=0) // can't dominate negative with negative or zero
    return NULL;
  vecr=(int *)malloc((size_t) (n*sizeof(int)));
  for(i=0;i<n;i++){
    vecr[i]=2*vecmid[i]-vecl[i];
    if(vecr[i]<0){ // only positive powers allowed
      free((char *)vecr);
      return NULL;
    }
  }
  return vecr;
}


// take a list of monomials in exponential form. Take a negative monomial mneg which needs to be dominated. Take a positive monomial mpos which can dominate mneg. Find all the other negative monomials which require mpos to dominate them. find the other ends of the dominating pairs. Continue recursively until the set grows no further. Check for all midpoints of left pairs and all midpoints of right pairs. 

// IMPORTANT NOTE: what is gathered
// doesn't have to include all the negative midpoints
// So it may not factorise perfectly

void gather(int **allt, int *cfs, long r, int n, int **neededby,int **needs,int *needexcl,long **lt,long **rt,long **neg,long **mid,long **all,int *numlt,int *numrt,int *numneg,int *nummid,int m0,int startind,int maxsz, int *numsofar, bool onlft, int *failflag, bool tight, bool excflag){
  int j,j1;
  int *vc;
  long s1;
  bool flg1,added;
  //fprintf(stderr, "entering gather...\n");
  //  cout << "numsofar0=" << *numsofar << endl;
  if((*failflag)==1){// failed somewhere
    cout << "failed on entry\n";
    return;
  }
  //  cout << "here1\n";
  if((*numsofar)>maxsz){// too many terms
    (*failflag)=1;
    cout << "entered with too many terms\n";
    return;
  }

  // enter routine; check if a given term has all the necessary midterms. 
  // If so add these midterms and the term itself; if not then return
  // if the term is the first in a left or right list, the check is trivially successful

  if(onlft){
    //first check for/add midpoints with existing left terms

    for(j1=0;j1<(*numlt);j1++){ // all must exist
      vc=midpt(allt[m0],allt[(*lt)[j1]],cfs[m0],cfs[(*lt)[j1]],n);
      if(!vc || (s1=binisearch(allt,r,vc,n))<0 || cfs[s1]<=0 || isinllist(s1, (*lt), (*numlt)) || isinllist(s1, (*rt), (*numrt)) || (abs(cfs[s1]) < splt(cfs[m0], cfs[(*lt)[j1]]))){// midpoint not found or negative or already a left or right term or has insufficient weight
    	if(vc)
    	  free((char*)vc);
    	if(tight || excflag){(*failflag)=1;/* cout << "major fail\n"; */}else{(*failflag)=-1;}
    	//cout << "failing to find midterm\n";printvec1(allt[m0],n);printvec1(allt[(*lt)[j1]],n);
    	return;
      }
      if(vc)
    	free((char*)vc);
    }

    // all midterms good

    for(j1=0;j1<(*numlt);j1++){
      vc=midpt(allt[m0],allt[(*lt)[j1]],cfs[m0],cfs[(*lt)[j1]],n);
      s1=binisearch(allt,r,vc,n);
      if(!isinllist(s1, (*mid), (*nummid))){
    	(*numsofar)++;
    	if((*numsofar)>maxsz){(*failflag)=1;fprintf(stderr, "Clique got too large.\n");return;}
    	addnewto1Dlarray(mid,(*nummid),s1);(*nummid)++;
    	addnewto1Dlarray(all, (*numsofar)-1,s1);
    	//cout << "adding midterm\n   "; printvec1(allt[s1],n);
      }
      free((char*)vc);
    }

    // now add the actual new term

   //cout << "here7\n";
    (*numsofar)++;
    if((*numsofar)>maxsz){(*failflag)=1;fprintf(stderr, "Clique got too large.\n");return;}

    addnewto1Dlarray(lt, (*numlt),m0);(*numlt)++;
    addnewto1Dlarray(all, (*numsofar)-1,m0);
    //cout << "adding on left\n   "; printvec1(allt[m0],n);

  }
  else if(!onlft){
    //cout << "here8\n";

    // first check for/add midpoints with existing rightterms
    for(j1=0;j1<(*numrt);j1++){
      vc=midpt(allt[m0],allt[(*rt)[j1]],cfs[m0],cfs[(*rt)[j1]],n);
      if(!vc || (s1=binisearch(allt,r,vc,n))<0 || cfs[s1]<=0 || isinllist(s1, (*lt), (*numlt)) || isinllist(s1, (*rt), (*numrt)) || (abs(cfs[s1]) < splt(cfs[m0], cfs[(*rt)[j1]]))){// midpoint not found or negative or already a left or right term or has insufficient weight
    	if(vc)
    	  free((char*)vc);
    	if(tight || excflag){(*failflag)=1;/* cout << "major fail\n"; */}else{(*failflag)=-1;}
	//cout << "failing to find midterm\n";printvec1(allt[m0],n);printvec1(allt[(*rt)[j1]],n);
    	return;
      }
      if(vc)
    	free((char*)vc);
    }

    // all midterms good

    for(j1=0;j1<(*numrt);j1++){
      vc=midpt(allt[m0],allt[(*rt)[j1]],cfs[m0],cfs[(*rt)[j1]],n);
      s1=binisearch(allt,r,vc,n);
      if(!isinllist(s1, (*mid), (*nummid))){
    	(*numsofar)++;
    	if((*numsofar)>maxsz){(*failflag)=1;fprintf(stderr, "Clique got too large.\n");return;}

    	addnewto1Dlarray(mid,(*nummid),s1);(*nummid)++;
    	addnewto1Dlarray(all, (*numsofar)-1,s1);
    	//cout << "adding midterm\n   "; printvec1(allt[s1],n);
      }
      free((char*)vc);
    }


    // now add the actual new term

    (*numsofar)++;
    if((*numsofar)>maxsz){(*failflag)=1;fprintf(stderr, "Clique got too large.\n");return;}

    addnewto1Dlarray(rt, (*numrt),m0);(*numrt)++;
    addnewto1Dlarray(all, (*numsofar)-1,m0);
    //cout << "adding on right\n   "; printvec1(allt[m0],n);
  }

  //cout << "here2\n";  
  //
  //
  // Now branch out from the current term
  //
  //

  for(j=startind;j<neededby[m0][0]+1;j+=2){
    //cout << "herek\n";
    //(*failflag)=0;// shouldn't be needed
    added=0;
    flg1=0;if(!excflag ||needexcl[neededby[m0][j]]){flg1=1;}
    if(onlft && !isinllist(neededby[m0][j+1], (*rt), (*numrt)) && flg1){// new right term
      /* if(needs[neededby[m0][j]][0]==2)//has to be tight - negative term allows only one domination */
      /* 	gather(allt,cfs,r,n,neededby,needs,needexcl,lt,rt,neg,mid,all,numlt,numrt,numneg,nummid,neededby[m0][j+1],1,maxsz,numsofar,0,failflag,1, excflag); */
      /* else */
      gather(allt,cfs,r,n,neededby,needs,needexcl,lt,rt,neg,mid,all,numlt,numrt,numneg,nummid,neededby[m0][j+1],1,maxsz,numsofar,0,failflag,tight, excflag);
      added=1;
      if((*failflag)==1)
	return;
      /* if((*failflag) && needexcl[neededby[m0][j]]){ */
      /* 	/\* printvec(allt[m0],n); *\/ */
      /* 	/\* printvec(allt[neededby[m0][j]],n); *\/ */
      /* 	/\* printvec(allt[neededby[m0][j+1]],n); *\/ */
      /* 	/\* printvec(needs[neededby[m0][j]],needs[neededby[m0][j]][0]+1); *\/ */
      /* 	/\* //exit(0); *\/ */
      /* 	(*failflag)==1; return; */
      /* } */
    }
    else if(!onlft && !isinllist(neededby[m0][j+1], (*lt), (*numlt)) && flg1){// add a left term
      /* if(needs[neededby[m0][j]][0]==2)//has to be tight - negative term allows only one domination */
      /* 	gather(allt,cfs,r,n,neededby,needs,needexcl,lt,rt,neg,mid,all,numlt,numrt,numneg,nummid,neededby[m0][j+1],maxsz,numsofar,1,failflag,1, excflag); */
      /* else */
      gather(allt,cfs,r,n,neededby,needs,needexcl,lt,rt,neg,mid,all,numlt,numrt,numneg,nummid,neededby[m0][j+1],1,maxsz,numsofar,1,failflag,tight, excflag);
      added=1;
      if((*failflag)==1)
	return;
      /* if((*failflag) && needexcl[neededby[m0][j]]){ */
      /* 	/\* printvec(allt[m0],n); *\/ */
      /* 	/\* printvec(allt[neededby[m0][j]],n); *\/ */
      /* 	/\* printvec(allt[neededby[m0][j+1]],n); *\/ */
      /* 	/\* printvec(needs[neededby[m0][j]],needs[neededby[m0][j]][0]+1); *\/ */
      /* 	/\* //exit(0); *\/ */
      /* 	(*failflag)==1; return; */
      /* } */
    }
    // should only do this if previous steps were fully successful
    if(added && !(*failflag) && !isinllist(neededby[m0][j], (*neg), (*numneg))){ // negative midpoints if endpoints added
      (*numsofar)++;
      if((*numsofar)>maxsz){(*failflag)=1;fprintf(stderr, "Clique got too large.\n");return;}

      addnewto1Dlarray(neg, (*numneg),neededby[m0][j]);(*numneg)++;
      addnewto1Dlarray(all, (*numsofar)-1,neededby[m0][j]);


      //cout << "adding negative term\n   "; printvec1(allt[neededby[m0][j]],n);
    }

  }
  //cout << "numsofar5=" << *numsofar << endl;
  //fprintf(stderr, "exiting gather...\n");
  return;

}

//
//what is the size of the intersection of a set of monomials?
//

int numshare(int **allt, long r, long *inds, int num, int n){
  int i,j;
  int share=0;
  bool flg;
  for(j=0;j<n;j++){
    flg=1;
    for(i=1;i<num;i++){
      if(allt[inds[i]][j]!=allt[inds[0]][j]){
	flg=0;break;
      }
    }
    if(flg)
      share++;
  }
  return share;
}

int hasonlyfullblocks(long *all, int numall, long **blocks, int *numinblock, int numblocks){
  int i,j;
  int numfull=0;
  bool flag;
  for(i=0;i<numblocks;i++){
    flag=1;
    for(j=0;j<numinblock[i];j++){
      if(!(isinllist(blocks[i][j], all,numall))){//absent
	if(flag==1 && j>0)// half a block
	  return -1;
	flag=0;
      }
      else{//present
	if(flag==0)// half a block
	  return -1;
      }
    }
    if(flag)//fullblock
      numfull++;
  }
  return numfull;
}

//
// Give this access to the negative blocks so it can check if each clique found has complete negative blocks
// long **negindblocks, int *numinnegblock, int numnegblocks
//
// swtch=1: largest clique
// swtch=0: smallest clique
// swtch=2: first clique
// swtch=3: smallest width
// swtch=4: smallest width which still contains full blocks
// swtch=5: largest full block



int maxminclique(int **allt, int *cfs, long r, int n, int **neededby,int **needs,int *needexcl,long **all,int *numsofar,long **lt,long **rt,int *numlt, int *numrt, long **negindblocks, int *numinnegblock, int numnegblocks, bool tight, int swtch){
  long *neg=NULL,*mid=NULL;
  int numneg=0,nummid=0;
  int maxsz=1000;
  int failflag=0;
  long m0,mspec;
  int k,maxval=0;
  int minval=100;
  int k1,k2,sharebig=0;
  bool flag=0;
  int absnegmax=0;
  int maxfullblocks=-1;
  *numlt=0;*numrt=0;
  *all=NULL;*lt=NULL;*rt=NULL;
  cout << "pppp\n";
  for(m0=0;m0<r;m0++){
    //cout << "m0="<<m0<<endl;
    if(cfs[m0]>0 && neededby[m0][0]>0){
      (*numlt)=0;(*numrt)=0;numneg=0;nummid=0;
      (*numsofar)=0;failflag=0;
      gather(allt,cfs,r,n,neededby,needs,needexcl,lt,rt,&neg,&mid,all,numlt,numrt,&numneg,&nummid,m0,1,maxsz,numsofar,1,&failflag,tight,0);
      //cout << "gathered\n" << failflag <<" " << (*numsofar) << endl;

      if((*numsofar)>=3 && (failflag!=1)){
	k2=hasonlyfullblocks(*all,*numsofar,negindblocks,numinnegblock,numnegblocks);

 
	/* if(!failflag) */
	/* 	cout << (*numsofar) << endl; */
	if(swtch==1){// largest clique
	  if((*numsofar)>maxval){ 
	    maxval=(*numsofar);
	    mspec=m0;flag=1;// some clique found
	  }
	}
	else if(swtch==0){// smallest clique
	  if((*numsofar)<minval){
	    minval=(*numsofar);
	    mspec=m0;flag=1;// some clique found
	  }
	}
	else if(swtch==2){ // first real clique
	  mspec=m0;flag=1;
	}
	else if(swtch==3){
	  if((k1=numshare(allt, r, (*all), *numsofar,n))>sharebig){// smallest width
	    sharebig=k1;
	    fprintf(stderr, "sharebig = %d\n", sharebig);
	    mspec=m0;flag=1;
	  }
	}
	else if(swtch==4){
	  if(k2>0 && (k1=numshare(allt, r, (*all), *numsofar,n))>sharebig){
	    sharebig=k1;
	    fprintf(stderr, "sharebig = %d\n", sharebig);
	    mspec=m0;flag=1;
	  }
	}
	else if(swtch==5){
	  if(k2>maxfullblocks){
	    maxfullblocks=k2;
	    mspec=m0;flag=1;
	  }
	}

	if((*lt)){free((char*)(*lt));(*lt)=NULL;}if((*rt)){free((char*)(*rt));(*rt)=NULL;}if(neg){free((char*)neg);neg=NULL;}if(mid){free((char*)mid);mid=NULL;}if(*all){free((char*)(*all));(*all)=NULL;}
      }
    }
    if(flag && swtch==2)// first clique
      break;
  }

  if(flag){
    if(swtch==5)
      fprintf(stderr, "choosing clique with %d complete blocks.\n", maxfullblocks);
    m0=mspec;
    if(cfs[m0]>0 && neededby[m0][0]>0){// should definitely succeed!
      (*numlt)=0;(*numrt)=0;numneg=0;nummid=0;
      (*numsofar)=0;failflag=0;
      gather(allt,cfs,r,n,neededby,needs,needexcl,lt,rt,&neg,&mid,all,numlt,numrt,&numneg,&nummid,m0,1,maxsz,numsofar,1,&failflag,tight,0);
      //cout << "gathered\n" << failflag <<" " << (*numsofar) << endl;
      if((*numsofar) && (failflag!=1)){// should definitely succeed
	if(tight)
	  cout << "found a tight clique of size " << (*numsofar) << endl;
	else
	  cout << "found a loose clique of size " << (*numsofar) << endl;
	cout << "numlt = " << (*numlt) << endl;
	for(k=0;k<(*numlt);k++){
	  cout << "  " << cfs[(*lt)[k]] << ", "; printvec1(allt[(*lt)[k]],n);
	}
	cout << "numrt = " << (*numrt) << endl;
	for(k=0;k<(*numrt);k++){
	  cout << "  " << cfs[(*rt)[k]] << ", "; printvec1(allt[(*rt)[k]],n);
	}
	cout << "numneg = " << numneg << endl;
	
	for(k=0;k<numneg;k++){
	  cout << "  " << cfs[neg[k]] << ", "; printvec1(allt[neg[k]],n);
	  absnegmax=max(absnegmax,abs(cfs[neg[k]]));
	}
	cout << "nummid = " << nummid << endl;
	for(k=0;k<nummid;k++){
	  cout << "  " << cfs[mid[k]] << ", "; printvec1(allt[mid[k]],n);
	}
	if(neg){free((char*)neg);neg=NULL;}if(mid){free((char*)mid);mid=NULL;}
	return absnegmax;
	
	//exit(0);
      }
      /* fprintf(stderr, "numsofar=%d\n",*numsofar); */
      /* printlvec(*all,*numsofar); */
      
      if((*lt)){free((char*)(*lt));(*lt)=NULL;}if((*rt)){free((char*)(*rt));(*rt)=NULL;}if(neg){free((char*)neg);neg=NULL;}if(mid){free((char*)mid);mid=NULL;}if(*all){free((char*)(*all));(*all)=NULL;}
    }
  }
  if(tight){cout << "no tight cliques found\n";}else{cout << "no loose cliques found\n";}
  (*numsofar)=0;
  return absnegmax;
}

void allcliques(int **allt, int *cfs, long r, int n, int **neededby,int **needs, int *needexcl){
  long *lt=NULL,*rt=NULL,*neg=NULL,*mid=NULL,*all=NULL;
  int numlt=0,numrt=0,numneg=0,nummid=0;
  int maxsz=50;
  int numsofar=0;
  int failflag=0;
  long m0;
  int k;
  for(m0=0;m0<r;m0++){
    if(cfs[m0]>0 && neededby[m0][0]>0){
      numlt=0;numrt=0;numneg=0;nummid=0;
      numsofar=0;failflag=0;
      gather(allt,cfs,r,n,neededby,needs,needexcl,&lt,&rt,&neg,&mid,&all,&numlt,&numrt,&numneg,&nummid,m0,1,maxsz,&numsofar,1,&failflag,1,0);
      //cout << "gathered\n" << failflag <<" " << numsofar << endl;
      if(numsofar && (failflag!=1)){
	cout << "found a clique of size " << numsofar << endl;
	cout << "numlt = " << numlt << endl;
	cout << "numrt = " << numrt << endl;
	cout << "numneg = " << numneg << endl;
	cout << "nummid = " << nummid << endl;
	for(k=0;k<numlt;k++){
	  cout << cfs[lt[k]] << ", "; printvec1(allt[lt[k]],n);
	}
	for(k=0;k<numrt;k++){
	  cout << cfs[rt[k]] << ", "; printvec1(allt[rt[k]],n);
	}
	for(k=0;k<numneg;k++){
	  cout << cfs[neg[k]] << ", "; printvec1(allt[neg[k]],n);
	}
	for(k=0;k<nummid;k++){
	  cout << cfs[mid[k]] << ", "; printvec1(allt[mid[k]],n);
	}
	
	
	//exit(0);
      }
      if(lt){free((char*)lt);lt=NULL;}if(rt){free((char*)rt);rt=NULL;}if(neg){free((char*)neg);neg=NULL;}if(mid){free((char*)mid);mid=NULL;}if(all){free((char*)all);all=NULL;}
    }
  }
  return;
}

void ipswap(int **p, int i, int j){
  int *temp=p[i];
  p[i]=p[j];
  p[j]=temp;
}

void lswap(long *p, int i, int j){
  long temp=p[i];
  p[i]=p[j];
  p[j]=temp;
}

void lpswap(long **p, int i, int j){
  long *temp=p[i];
  p[i]=p[j];
  p[j]=temp;
}

//

void sizesort4(int *cliquesz, long **lt, long **rt, int *numlt, int *numrt, int *absnegmax, int left, int right)
{
  int i,last;

  if (left >= right)
    return;
  iswap(cliquesz,left,(left+right)/2);
  iswap(absnegmax,left,(left+right)/2);
  lpswap(lt,left,(left+right)/2);
  lpswap(rt,left,(left+right)/2);
  iswap(numlt,left,(left+right)/2);
  iswap(numrt,left,(left+right)/2);
  last = left;

  for (i=left+1; i<=right;i++){
    if (cliquesz[i]>cliquesz[left]){// or <
      last++;
      iswap(cliquesz,last,i);
      iswap(absnegmax,last,i);
      lpswap(lt,last,i);
      lpswap(rt,last,i);
      iswap(numlt,last,i);
      iswap(numrt,last,i);
    }
  }

  iswap(cliquesz,left,last);
  iswap(absnegmax,left,last);
  lpswap(lt,left,last);
  lpswap(rt,left,last);
  iswap(numlt,left,last);
  iswap(numrt,left,last);

  sizesort4(cliquesz, lt, rt, numlt,numrt, absnegmax, left, last-1);
  sizesort4(cliquesz, lt, rt, numlt,numrt, absnegmax, last+1,right);
  return;
}



int allcliques1(int **allt, int *cfs, long r, int n, int **neededby, int **needs, int *needexcl, long ***left, long ***right, int **numleft, int **numright, int **absnegmax, int maxsz, int maxgetcliques, bool tight, bool onlycmplt, bool excflag){
  long *lt=NULL,*rt=NULL,*neg=NULL,*mid=NULL,*all=NULL;
  int numlt=0,numrt=0,numneg=0,nummid=0;
  int *cliquequal=NULL;// "quality" of a clique
  int numsofar=0;
  int failflag=0;
  int *vc;
  long m0,s1;
  int k,k1;
  long *outvec;
  int numcliques=0;
  bool flg,cmplt;
  long **cliquearray=NULL;
  int *cliquesz=NULL;
  int tmpabsneg=0;
  int st,mx;
  int maxcliques=100, cliqueqtmp;
  for(m0=0;m0<r;m0++){
    if(excflag)
      mx=2;
    else
      mx=neededby[m0][0];
    if(cfs[m0]>0 && neededby[m0][0]>0){
      //cout << "checking: (mx="<<mx<<"): "; printvec1(allt[m0],n);
      for(st=1;st<mx;st+=2){
	numlt=0;numrt=0;numneg=0;nummid=0;
	numsofar=0;failflag=0;
	//fprintf(stderr, "entering gather.\n");
	gather(allt,cfs,r,n,neededby,needs,needexcl,&lt,&rt,&neg,&mid,&all,&numlt,&numrt,&numneg,&nummid,m0,st,maxsz,&numsofar,1,&failflag,tight,excflag);
	//fprintf(stderr, "exiting gather.\n");
	//cout << "gathered\n" << failflag <<" " << numsofar << endl;
	if(numsofar>=3 && numneg && (failflag!=1)){
	  //cout << "found clique\n";
	  unionlvec(lt,numlt,rt,numrt, &outvec);
	  flg=0;
	  for(k=0;k<numcliques;k++){
	    if(unordllistareeq(outvec,numlt+numrt,cliquearray[k],cliquesz[k])){// only check left and right halves
	      flg=1;
	      break;
	    }
	  }
	  if(!flg){// not already in array of cliques; add + build output info

	    if(onlycmplt){//only accept complete cliques
	      //is the clique complete?
	      cmplt=1;k=0;k1=0;
	      while(cmplt && k<numlt){
		while(cmplt && k1<numrt){
		  vc=midpt(allt[lt[k]],allt[rt[k1]],cfs[lt[k]],cfs[rt[k1]],n);
		  s1=binisearch(allt,r,vc,n);
		  if(s1<0 || !(isinllist(s1,neg,numneg)))
		    cmplt=0;
		  k1++;
		  if(vc){free((char*)vc);}
		}
		k++;
	      }
	    }

	    if(!onlycmplt || cmplt){
	      //fprintf(stderr, "found clique number %d (size: %d, neg: %d)\n", numcliques+1, numsofar, numneg);
	      cout << "exiting with l,r,mid,n, tot " << numlt << ", " << numrt << ", " <<  nummid << ", " <<  numneg << ", " <<  numsofar << endl;
	      cout << "found a clique of size " << numsofar << endl;
	      for(k=0;k<numlt;k++){
		cout << cfs[lt[k]] << ", "; printvec1(allt[lt[k]],n);
	      }
	      for(k=0;k<numrt;k++){
		cout << cfs[rt[k]] << ", "; printvec1(allt[rt[k]],n);
	      }
	      for(k=0;k<numneg;k++){
		cout << cfs[neg[k]] << ", "; printvec1(allt[neg[k]],n);
	      }
	      for(k=0;k<nummid;k++){
		cout << cfs[mid[k]] << ", "; printvec1(allt[mid[k]],n);
	      }

	      //	      fprintf(stderr, "found a clique\n");
	      addnewto1Darray(&cliquesz,numcliques,numlt+numrt);
	      addnewtolarray(&cliquearray,numcliques,outvec,numlt+numrt);
	      addnewtolarray(left,numcliques,lt,numlt);
	      addnewto1Darray(numleft,numcliques,numlt);
	      addnewtolarray(right,numcliques,rt,numrt);
	      addnewto1Darray(numright,numcliques,numrt);
	      cliqueqtmp=(10*numneg)/(numlt+numrt) + numneg/5;
	      //cliqueqtmp=(10*numneg)/(numlt+numrt) - numneg/10;
	      cout << "clique quality: " << cliqueqtmp << endl;
	      addnewto1Darray(&cliquequal,numcliques,cliqueqtmp);

	      for(k=0;k<numneg;k++)
		tmpabsneg=max(tmpabsneg,abs(cfs[neg[k]]));
	      addnewto1Darray(absnegmax,numcliques,tmpabsneg);

	      numcliques++;
	      if(maxgetcliques && numcliques>=maxgetcliques){
	      	fprintf(stderr, "got %d cliques\n", maxgetcliques);
		sizesort4(cliquequal, *left, *right, *numleft, *numright, *absnegmax, 0, numcliques-1);
		//sizesort4(cliquesz, *left, *right, *numleft, *numright, *absnegmax, 0, numcliques-1);
		if(numcliques>maxcliques){// free the ones which won't be returned
		  lkeeponly(*left, numcliques, maxcliques);
		  lkeeponly(*right, numcliques, maxcliques);
		}

	      	free((char*)outvec);
	      	if(lt){free((char*)lt);lt=NULL;}if(rt){free((char*)rt);rt=NULL;}if(neg){free((char*)neg);neg=NULL;}if(mid){free((char*)mid);mid=NULL;}if(all){free((char*)all);all=NULL;}
	      	if(cliquearray){lfree_i(cliquearray, numcliques);}
	      	if(cliquesz){free((char *)cliquesz);}
	      	if(cliquequal){free((char *)cliquequal);}

		if(numcliques>maxcliques)
		  return maxcliques;
		return numcliques;
	      }
	    }
	  }
	  free((char*)outvec);
	}
	if(lt){free((char*)lt);lt=NULL;}if(rt){free((char*)rt);rt=NULL;}if(neg){free((char*)neg);neg=NULL;}if(mid){free((char*)mid);mid=NULL;}if(all){free((char*)all);all=NULL;}
      }
    }
  }
  fprintf(stderr, "got %d cliques\n", numcliques);
  cout << "total cliques found = "<< numcliques << endl;

  sizesort4(cliquequal, *left, *right, *numleft, *numright, *absnegmax, 0, numcliques-1);
  //sizesort4(cliquesz, *left, *right, *numleft, *numright, *absnegmax, 0, numcliques-1);
  
  if(numcliques>maxcliques){// free the ones which won't be returned
    lkeeponly(*left, numcliques, maxcliques);
    lkeeponly(*right, numcliques, maxcliques);
  }
  if(cliquearray){lfree_i(cliquearray, numcliques);}
  if(cliquesz){free((char *)cliquesz);}
  if(cliquequal){free((char *)cliquequal);}

  if(numcliques>maxcliques)
    return maxcliques;
  return numcliques;
}



// Find all the ways that a negative monomial vec may be dominated 
// symmetrically by pairs of vectors from the list allt
// assume that allt is ordered
// needs is a vector of pairs which can dominate vec

float numsplits(long indx, int n, int **allt, int *cfs, long r, int **needcount, int **needexcl, int **needs, int **neededby, bool q){
  long s,s1;
  int *vc;
  int j,j1,j2;
  int spl=0;
  int powsplit=0;
  int cf=cfs[indx];
  //fprintf(stderr, "entering numsplits.\n");
  if(!q){
    cout << "---------------------" << endl;
    /* cout << cf << ", "; printvec1(allt[indx],n); */expltomonoprint(cfs[indx],allt[indx],n);
  }
  for(s=0;s<r;s++){
    vc=sympt(allt[indx],allt[s],cfs[s],n,0);
    if(vc && (s1=binisearch(allt+s,r-s,vc,n))>=0 && cfs[s+s1]>0){
      (needs[indx][0])+=2;
      j=(needs[indx][0])+1;
      needs[indx]=(int*) realloc((needs[indx]), sizeof(int)*(j)); 
      
      (neededby[s][0])+=2;
      (neededby[s+s1][0])+=2;
      j1=(neededby[s][0])+1;
      j2=(neededby[s+s1][0])+1;
      neededby[s]=(int*) realloc((neededby[s]), sizeof(int)*(j1));
      neededby[s+s1]=(int*) realloc((neededby[s+s1]), sizeof(int)*(j2));

      powsplit+=splt(cfs[s],cfs[s+s1]);
      needs[indx][j-2]=s;needs[indx][j-1]=s+s1;
      neededby[s][j1-2]=indx;neededby[s][j1-1]=s+s1;
      neededby[s+s1][j2-2]=indx;neededby[s+s1][j2-1]=s;
      ((*needcount)[s])++;
      ((*needcount)[s+s1])++;
      if(!q){
	/* cout << cfs[s]<< ", ";  printvec1(allt[s],n);  */expltomonoprint(cfs[s],allt[s],n);
	/* cout << cfs[s+s1]<< ", ";  printvec1(allt[s+s1],n);  */expltomonoprint(cfs[s+s1],allt[s+s1],n);
      }
      spl++;

    }
    if(vc)
      free((char*)vc);
  }
  if(powsplit==abs(cf)){ //no choice
    (*needexcl)[indx]=1;
  }
  if(!q){
    cout << "powsplit = " << powsplit << endl;
    if(powsplit<abs(cf)){
      cout << "Warning: powers don't add up.\n";
      fprintf(stderr, "failure to split.\n");
      //exit(0);
    }
    else if(powsplit==abs(cf))
      cout << "All splits needed.\n";
  }
  if(powsplit<abs(cf)) // basically a failure to split in this way
    return -1.0;
  //fprintf(stderr, "exiting numsplits.\n");
  return sqrt(((float)powsplit)/((float)(abs(cf)))-0.99);
}

// Do all negative monomials in allt split?

bool allhavesplits(int **allt, int *cfs, long r, int n, bool q){
  long s,s1,s2;
  int *vc;
  int powsplit, cf;
  for(s2=0;s2<r;s2++){
    cf=cfs[s2];
    if(cf<0){
      powsplit=0;
      if(!q){
	cout << "---------------------" << endl;
	/* cout << cf << ", "; printvec1(allt[s2],n); */expltomonoprint(cfs[s2],allt[s2],n);
      }
      for(s=0;s<r;s++){
	vc=sympt(allt[s2],allt[s],cfs[s],n,0);
	if(vc && (s1=binisearch(allt+s,r-s,vc,n))>=0 && cfs[s+s1]>0){
	  powsplit+=splt(cfs[s],cfs[s+s1]);
	  if(!q){
	    /* cout << cfs[s]<< ", ";  printvec1(allt[s],n);  */expltomonoprint(cfs[s],allt[s],n);
	    /* cout << cfs[s+s1]<< ", ";  printvec1(allt[s+s1],n);  */expltomonoprint(cfs[s+s1],allt[s+s1],n);
	  }
	}
	if(vc){free((char*)vc);}
      }
 
      if(!q){
	cout << "powsplit = " << powsplit << endl;
	if(powsplit<abs(cf)){
	  cout << "Warning: powers don't add up.\n";
	  fprintf(stderr, "failure to split.\n");
	}
	else if(powsplit==abs(cf))
	  cout << "All splits needed.\n";
      }
      if(powsplit<abs(cf)) // failure to split
	return 0;
    }
  }
  return 1;
}


ex remsqs(ex tmp, int **explst1, int *cflst1, long r1, int **explst2, int *cflst2, long r2, int numv, int *used1, int *used2, int clq1, int clq2, bool q){
  ex tmp2=tmp;
  ex tmp1,extot=0;
  bool swtch=0;
  long s, m;
  if(!q){cout << endl << "current state of the squares:\n" << endl;}
  for(s=0;s<r1;s++){
    if(used1[s]==-1){
      if(!swtch){
	extot+=expltomono(1,explst1[s],numv);swtch=1;
      }
      else
	extot-=expltomono(1,explst1[s],numv);
    }
  }
  tmp1=expand(pow(extot,2));tmp2-=tmp1;
  if(extot!=0 && !q){cout << "(incompatible square):  " << pow(extot,2) << endl;}

  extot=0;swtch=0;
  for(s=0;s<r2;s++){
    if(used2[s]==-2){
      if(!swtch){
	extot+=expltomono(1,explst2[s],numv);swtch=1;
      }
      else
	extot-=expltomono(1,explst2[s],numv);
    }
  }
  tmp1=expand(pow(extot,2));tmp2-=tmp1;
  if(extot!=0 && !q){cout << "(incompatible square):  " << pow(extot,2) << endl;}

  for(m=1;m<=clq1;m++){// each clique from Re
    extot=0;
    for(s=0;s<r1;s++){// or should we use coefficients 1?
      if(used1[s]==m)
	extot+=expltomono(cflst1[s],explst1[s],numv);
    }
    tmp1=expand(pow(extot,2));tmp2-=tmp1;
    if(!q){cout << pow(extot,2) << endl;}
  }

  for(m=1;m<=clq2;m++){// each clique from Im
    extot=0;
    for(s=0;s<r2;s++){// or should we use coefficients 1?
      if(used2[s]==m)
	extot+=expltomono(cflst2[s],explst2[s],numv);
    }
    tmp1=expand(pow(extot,2));tmp2-=tmp1;
    if(!q){cout << pow(extot,2) << endl;}
  }
  if(!q){cout << endl;}
  return tmp2;
}


void extolst(ex tmp, int numv, int ***explst, int **cflst, long *r){
  int **lst;

  polytointlist0(tmp, &lst, cflst, r);
  (*explst)=moncflisttoexplist0(lst,*cflst,*r,numv+1,1);// redefines cflst
  //elongate(explst1,cflst1,r1,numv);
  ifree_l(lst, *r);
  return;
}


// count the negative terms in an expression

int negcnt(ex tmp){
  int **lst,*cflst;
  long r, s;
  int negcount=0;
  polytointlist0(tmp, &lst, &cflst, &r);
  //extolst(tmp, numv, &explst, &cflst, &r);
  for(s=0;s<r;s++){
    if(cflst[s]<0)
      negcount++;
  }
  ifree_l(lst,r); free((char*)cflst); 
  return negcount;
}

void cliquediff(int *used1, long r1, int *used2, long r2, int *used1tmp, int *used2tmp, int *used1t, int *used2t){
  long s;
  for(s=0;s<r1;s++){
    if(!used1tmp[s])
      used1t[s]=used1[s];
    else
      used1t[s]=0;
  }

  for(s=0;s<r2;s++){
    if(!used2tmp[s])
      used2t[s]=used2[s];
    else
      used2t[s]=0;
  }

  return;

}

int biggestclique(int *used1, long r1, int *used2, long r2, int *used1tmp, int *used2tmp, int *clqnum){
  long s,s1;
  int val,k=0;
  int sz,szold=0;
  for(s=0;s<r1;s++){used1tmp[s]=used1[s];}for(s=0;s<r2;s++){used2tmp[s]=used2[s];}
  //count cliques

  for(s=0;s<r1;s++){
    if((val=used1tmp[s])){
      sz=0;
      for(s1=s;s1<r1;s1++){
	if(used1tmp[s1]==val){
	  used1tmp[s1]=0;sz++;
	}
      }
      if(sz>szold){szold=sz;(*clqnum)=val;k=1;}
    }
  }

  for(s=0;s<r2;s++){
    if((val=used2tmp[s])){
      sz=0;
      for(s1=s;s1<r2;s1++){
	if(used2tmp[s1]==val){
	  used2tmp[s1]=0;sz++;
	}
      }
      if(sz>szold){szold=sz;(*clqnum)=val;k=2;}
    }
  }
  fprintf(stderr, "szold=%d\n", szold);
  return k;
}

bool clqsremleavessplits(ex tmp, int **explst1, int *cflst1, int r1, int *used1, int **explst2, int *cflst2, int r2, int *used2, int numv){
  ex tmp2;
  bool retval=0;
  long s;
  int k=0,cq1,cq2;
  int clqnum;
  int *used1tmp=(int *)malloc((size_t) (r1*sizeof(int)));
  int *used2tmp=(int *)malloc((size_t) (r2*sizeof(int)));
  int *used1t=(int *)malloc((size_t) (r1*sizeof(int)));
  int *used2t=(int *)malloc((size_t) (r2*sizeof(int)));
  int **explst, *cflst;
  long r;

  k=biggestclique(used1, r1, used2, r2, used1tmp, used2tmp,&clqnum);
  if(k==1){
    for(s=0;s<r1;s++){if(used1[s]==clqnum){used1tmp[s]=1;}}
    cq1=1;cq2=0;
  }
  else if(k==2){
    for(s=0;s<r2;s++){if(used2[s]==clqnum){used2tmp[s]=1;}}
    cq1=0;cq2=1;
  }
  fprintf(stderr, "clqnum=%d", clqnum);
  tmp2=remsqs(tmp, explst1, cflst1, r1, explst2, cflst2, r2, numv, used1tmp, used2tmp, cq1, cq2, 0);
  /* extolst(tmp2, numv, &explst, &cflst, &r); */
  /* fprintf(stderr, "negative count = %d\n", negcnt(tmp2)); */
  /* if(allhavesplits(explst, cflst, r, numv, 1)) */
  /*   retval=1; */
  /* ifree_l(explst,r); free((char*)cflst);  */

  cliquediff(used1, r1, used2, r2, used1tmp, used2tmp, used1t, used2t);
  k=biggestclique(used1t, r1, used2t, r2, used1tmp, used2tmp,&clqnum);
  if(k==1){
    for(s=0;s<r1;s++){if(used1t[s]==clqnum){used1tmp[s]=1;}}
    cq1=1;cq2=0;
  }
  else if(k==2){
    for(s=0;s<r2;s++){if(used2t[s]==clqnum){used2tmp[s]=1;}}
    cq1=0;cq2=1;
  }
  fprintf(stderr, "clqnum=%d", clqnum);
  tmp2=remsqs(tmp, explst1, cflst1, r1, explst2, cflst2, r2, numv, used1tmp, used2tmp, cq1, cq2, 0);

  cliquediff(used1t, r1, used2t, r2, used1tmp, used2tmp, used1t, used2t);
  k=biggestclique(used1t, r1, used2t, r2, used1tmp, used2tmp,&clqnum);
  if(k==1){
    for(s=0;s<r1;s++){if(used1t[s]==clqnum){used1tmp[s]=1;}}
    cq1=1;cq2=0;
  }
  else if(k==2){
    for(s=0;s<r2;s++){if(used2t[s]==clqnum){used2tmp[s]=1;}}
    cq1=0;cq2=1;
  }
  fprintf(stderr, "clqnum=%d", clqnum);
  tmp2=remsqs(tmp, explst1, cflst1, r1, explst2, cflst2, r2, numv, used1tmp, used2tmp, cq1, cq2, 0);

  extolst(tmp2, numv, &explst, &cflst, &r);
  fprintf(stderr, "negative count = %d\n", negcnt(tmp2));
  if(allhavesplits(explst, cflst, r, numv, 1))
    retval=1;
  ifree_l(explst,r); free((char*)cflst); 

  free((char*)used1tmp);free((char*)used2tmp);
  free((char*)used1t);free((char*)used2t);
  return retval;
}

int *sum(int *a, int *b, int n){
  int i;
  int *out=(int *)malloc((size_t) (n*sizeof(int)));
  for(i=0;i<n;i++)
    out[i]=a[i]+b[i];
  return out;

}

int honestsplits(long indx, int n, int **allt, int *cfs, long r, int **lst1, int *cfs1, int r1, int **lst2, int *cfs2, int r2, int *used1, int *used2, int **needcount1, int **needs1, int **neededby1, int **needcount2, int **needs2, int **neededby2, int **needexcl, bool q){
  long s,s1;
  int *vc;
  int j,j1,j2;
  int spl=0;
  int powsplit=0;
  int cf=cfs[indx];
  //fprintf(stderr, "entering honestsplits.\n");
  if(!q){
    cout << "---------------------" << endl;
    /* cout << cf << ", "; printvec1(allt[indx],n); */expltomonoprint(cfs[indx],allt[indx],n);
  }
  for(s=0;s<r1;s++){
    for(s1=s+1;s1<r1;s1++){
      if(cfs1[s]*cfs1[s1]<0 && !used1[s] && !used1[s1]){
	vc=sum(lst1[s],lst1[s1],n);
	if(areequal(vc, allt[indx], n)){
	  (needs1[indx][0])+=2;j=(needs1[indx][0])+1;
	  needs1[indx]=(int*) realloc((needs1[indx]), sizeof(int)*(j)); 
      
	  (neededby1[s][0])+=2;(neededby1[s1][0])+=2;
	  j1=(neededby1[s][0])+1;j2=(neededby1[s1][0])+1;
	  neededby1[s]=(int*) realloc((neededby1[s]), sizeof(int)*(j1));
	  neededby1[s1]=(int*) realloc((neededby1[s1]), sizeof(int)*(j2));
	  powsplit+=2*abs(cfs1[s]*cfs1[s1]);;
	  needs1[indx][j-2]=s;needs1[indx][j-1]=s1;
	  neededby1[s][j1-2]=indx;neededby1[s][j1-1]=s1;
	  neededby1[s1][j2-2]=indx;neededby1[s1][j2-1]=s;
	  ((*needcount1)[s])++;((*needcount1)[s1])++;
	  if(!q){
	    /* cout << cfs[s]<< ", ";  printvec1(allt[s],n);  */expltomonoprint(1,lst1[s],n);
	    /* cout << cfs[s+s1]<< ", ";  printvec1(allt[s+s1],n);  */expltomonoprint(1,lst1[s1],n);
	  }
	  spl++;

	}
	if(vc)
	  free((char*)vc);
      }
    }
  }

  for(s=0;s<r2;s++){
    for(s1=s+1;s1<r2;s1++){
      if(cfs2[s]*cfs2[s1]<0 && !used2[s] && !used2[s1]){
	vc=sum(lst2[s],lst2[s1],n);
	if(areequal(vc, allt[indx], n)){
	  (needs2[indx][0])+=2;j=(needs2[indx][0])+1;
	  needs2[indx]=(int*) realloc((needs2[indx]), sizeof(int)*(j)); 
      
	  (neededby2[s][0])+=2;(neededby2[s1][0])+=2;
	  j1=(neededby2[s][0])+1;j2=(neededby2[s1][0])+1;
	  neededby2[s]=(int*) realloc((neededby2[s]), sizeof(int)*(j1));
	  neededby2[s1]=(int*) realloc((neededby2[s1]), sizeof(int)*(j2));
	  powsplit+=2*abs(cfs2[s]*cfs2[s1]);
	  needs2[indx][j-2]=s;needs2[indx][j-1]=s1;
	  neededby2[s][j1-2]=indx;neededby2[s][j1-1]=s1;
	  neededby2[s1][j2-2]=indx;neededby2[s1][j2-1]=s;
	  ((*needcount2)[s])++;((*needcount2)[s1])++;
	  if(!q){
	    /* cout << cfs[s]<< ", ";  printvec1(allt[s],n);  */expltomonoprint(1,lst2[s],n);
	    /* cout << cfs[s+s1]<< ", ";  printvec1(allt[s+s1],n);  */expltomonoprint(1,lst2[s1],n);
	  }
	  spl++;
	}
	if(vc)
	  free((char*)vc);
      }
    }
  }
  if(!q){cout << "powsplit = " << powsplit << endl;}
  if(powsplit==abs(cf)){ //no choice
    if(!q){cout << "All splits needed.\n";}
    (*needexcl)[indx]=1;
  }

  if(powsplit<abs(cf)){
    if(!q){
      cout << "Warning: powers don't add up.\n";
      fprintf(stderr, "failure to split (honestsplits).\n");
    }
    return 0;
    //exit(0);
  }

 //fprintf(stderr, "exiting honestsplits.\n");
 return spl;
}




// 
// returns true if the negative term allt[indx] is not the midpoint of a pair
//

bool isorphan(long indx, int **allt, int *cfs, long r, int n){
  long s,s1;
  int *vc;
  int powsplit=0;
  int cf=cfs[indx];
  for(s=0;s<r;s++){
    vc=sympt(allt[indx],allt[s],cfs[s],n,0);
    if(vc && (s1=binisearch(allt+s,r-s,vc,n))>=0 && cfs[s+s1]>0){
      powsplit+=splt(cfs[s],cfs[s+s1]);
      //cout << indx << ": " << cf << ", "; printvec1(allt[indx],n);
      //cout << "split successfully " << powsplit << endl;
    }
    if(vc)
      free((char*)vc);
  }
  if(powsplit<abs(cf)){
    //cout << indx << ": " << cf << ", "; printvec1(allt[indx],n);cout << "split failed " << powsplit << endl;
    return 1;
  }
  return 0;
}

bool hasorphans(int **allt, int *cfs, long r, int n){
  long s;
  //cout << "entering hasorphans\n";
  for(s=0;s<r;s++){
    if(cfs[s]<0 && isorphan(s, allt, cfs, r, n)){
      //cout << "orphan: "<< s << ": "; printvec1(allt[s],n); 
      //cout << "exiting hasorphans\n";
      return 1;
    }
  }
  //cout << "exiting hasorphans\n";
  return 0;
}






//
// does removing a particular perfect square leave "orphaned" negative terms? 
// assume that the square exists in the list
//

bool sqremleavesorphans(int **allt, int *cfs, long r, int n, int **lt, int numlt, int **rt, int numrt, bool loose){
  int i,j;
  long s,s1;
  int cfs1[r];
  int *vc;
  for(s=0;s<r;s++)
    cfs1[s]=cfs[s];

  /* if(hasorphans(allt,cfs1,r,n)){ */
  /*   fprintf(stderr, "already has orphans!\n"); */
  /*   exit (0); */
  /* } */

  for(i=0;i<numlt;i++){
    for(j=0;j<numlt;j++){
      vc=sum(lt[i],lt[j],n);
      //cout << "removing \n"; printvec1(vc,n);
      if((s1=binisearch(allt,r,vc,n))<0 || cfs1[s1]<=0){// assume a proper clique
	//fprintf(stderr, "not a clique\n");printvec(vc,n);
	free((char*)vc);
	return 1;
	//exit(0);
      }
      free((char*)vc);
      cfs1[s1]--;
      //     cout << cfs[s1] << ", ";printvec1(allt[s1],n);
    }
  }
  for(i=0;i<numrt;i++){
    for(j=0;j<numrt;j++){
      vc=sum(rt[i],rt[j],n);
      //printvec1(vc,n);
      if((s1=binisearch(allt,r,vc,n))<0 || cfs1[s1]<=0){// assume a proper clique
	//fprintf(stderr, "not a clique\n");printvec(vc,n);
	free((char*)vc);
	return 1;
	//exit(0);
      }
      free((char*)vc);
      //cout << cfs[s1] << ", ";printvec1(allt[s1],n);
      cfs1[s1]--;
    }
  }
  if(loose){// need not be there
    for(i=0;i<numlt;i++){
      for(j=0;j<numrt;j++){
	vc=sum(lt[i],rt[j],n);
	//printvec1(vc,n);
	if((s1=binisearch(allt,r,vc,n))>=0)// need not be there
	  cfs1[s1]+=2;
	free((char*)vc);
      }
    }
  }
  else{ // only squares which actually exist in the full expansion
    for(i=0;i<numlt;i++){
      for(j=0;j<numrt;j++){
	vc=sum(lt[i],rt[j],n);
	//printvec1(vc,n);
	if((s1=binisearch(allt,r,vc,n))<0 || cfs1[s1]>=0)// not there or not needed
	  return 1;
	else
	  cfs1[s1]+=2;
	free((char*)vc);
      }
    }

  }
  //cout << "from \n"; listprint(cfs,allt,n,r);

  if(hasorphans(allt,cfs1,r,n)){
    //fprintf(stderr, "leaves orphans\n");
    return 1;
  }
  return 0;

}

void veccp(int *cflst,long r,int *cflst1){
  long s;
  for(s=0;s<r;s++)
    cflst1[s]=cflst[s];
  return;
}

void veccp0(int *vecin,int r,int *vecout){
  int s;
  for(s=0;s<r;s++)
    vecout[s]=vecin[s];
  return;
}


void addtocliquelist(int **lt,int k1,int **rt, int k2, int numv, int ****leftlist, int **leftsz, int ****rightlist, int **rightsz, int *numclqs){

  //fprintf(stderr, "entering atcl\n");
  addnewblock(leftlist, *numclqs, lt, k1, numv);
  addnewto1Darray(leftsz,*numclqs,k1);
  addnewblock(rightlist, *numclqs, rt, k2, numv);
  addnewto1Darray(rightsz,*numclqs,k2);
  (*numclqs)++;
  //fprintf(stderr, "exiting atcl\n");
  return;
}
bool issublist(int **biglist, long r, int **smalllist, int sz, int n){
  int i;
  for(i=0;i<sz;i++)
    if(binisearch(biglist,r,smalllist[i],n)<0)
      return 0;
  return 1;
}

// assume ordered

bool listseq(int **a, int **b, int n, int k){
  int j;
  for(j=0;j<k;j++){
    if(!areequal(a[j], b[j], n))
      return 0;
  }
  return 1;
}

//assume all are ordered
bool incliquelist(int **lt,int k1,int **rt, int k2, int numv, int ***leftlist, int *leftsz, int ***rightlist, int *rightsz, int numclqs){
  int i;
  for(i=0;i<numclqs;i++){// each existing clique
    if(k1==leftsz[i] && k2==rightsz[i]){// correctsz
      if(listseq(lt, leftlist[i],numv,k1) && listseq(rt, rightlist[i],numv,k2))
	return 1;
    }
  }
  return 0;
}

bool sqmidsneeded(int **allt, int *cfs, long r, int numv,int **lt, int numlt, int **rt, int numrt, bool allswtch){
  int i,j;
  int *vc;
  long s1;
  bool someflg=0;// some mid-point is needed
  //fprintf(stderr, "entering sqmid...\n");
  for(i=0;i<numlt;i++){
    for(j=0;j<numrt;j++){
      vc=sum(lt[i],rt[j],numv);
      if(vc && ((s1=binisearch(allt,r,vc,numv))<0 || cfs[s1]>0)){// not needed
	if(allswtch){
	  free((char*)vc);//fprintf(stderr, "exiting sqmid 1...\n");
	  return 0;
	}
      }
      else
	someflg=1;
      if(vc){free((char*)vc);}
    }
  }
  //fprintf(stderr, "exiting sqmid 2...\n");
  return someflg;
}



bool disj(int **lt,int **rt, int k1, int k2, int numv){
  int i;
  //  fprintf(stderr, "entering disj\n");
  for(i=0;i<k2;i++){
    if(isinarray3(lt,k1,rt[i],numv)>=0){
      //fprintf(stderr, "exiting disj1\n");
      return 0;
    }
  }
  //fprintf(stderr, "exiting disj2\n");
  return 1;
}
// To make a square from a set, choose disjoint left-set and right-set
// Left and right are not uniquely defined: two squares are equal 


void goodsqs(int **allt, int *cfs, long r, int n, int **ltall, int numltall, int **rtall, int numrtall, int ****leftlist,int ****rightlist,int **leftsz,int **rightsz,int *numclqs,int k1max, int k2max, bool allmids){
  int i,j;
  int **lt;
  int **rt;
  int k1,k2;
  int *xc, *yc;
  int combflag,combflag1;

  for(k1=1;k1<=k1max;k1++){
    fprintf(stderr, "k1=%d\n", k1);
    lt = (int **)malloc((size_t) (k1*sizeof(int*)));
    for(i=0;i<k1;i++)
      lt[i]=(int *)malloc((size_t) (n*sizeof(int)));
    xc=(int *)malloc((size_t) (k1*sizeof(int)));

    for(k2=1;k2<=k2max;k2++){
      fprintf(stderr, "k2=%d\n", k2);
      rt = (int **)malloc((size_t) (k2*sizeof(int*)));
      for(i=0;i<k2;i++)
	rt[i]=(int *)malloc((size_t) (n*sizeof(int)));
      yc=(int *)malloc((size_t) (k2*sizeof(int)));

      firstcomb(xc,numltall,k1);combflag=1;
      while(combflag==1){
	for(j=0;j<k1;j++)
	  veccp0(ltall[xc[j]],n,lt[j]);
	combflag1=firstcombfrom(yc,numltall,k2,xc[0]+1);
	while(combflag1==1){
	  for(j=0;j<k2;j++)
	    veccp0(ltall[yc[j]],n,rt[j]);
	  if(disj(lt,rt,k1,k2,n) && !incliquelist(lt,k1,rt,k2,n,*leftlist,*leftsz,*rightlist,*rightsz,*numclqs) && !sqremleavesorphans(allt, cfs,r,n,lt,k1,rt,k2,1) && sqmidsneeded(allt,cfs,r,n,lt,k1,rt,k2,allmids)){//all midpoints needed
	    addtocliquelist(lt,k1,rt,k2,n,leftlist,leftsz,rightlist,rightsz,numclqs);
	    fprintf(stderr, "squares found: %d\n",*numclqs);
	    /* cout << "good pair:\n"; printvec1(lt[0],n); printvec1(rt[0],n); cout << endl; */
	  }
	  combflag1=nextcomb(yc,numltall,k2);
	}
	combflag=nextcomb(xc,numltall,k1);
      }
      free((char*)yc);ifree(rt,k2);
    }
    free((char*)xc);ifree(lt,k1);
  }


  for(k1=1;k1<=k1max;k1++){
    fprintf(stderr, "k1(2)=%d\n", k1);
    lt = (int **)malloc((size_t) (k1*sizeof(int*)));
    for(i=0;i<k1;i++)
      lt[i]=(int *)malloc((size_t) (n*sizeof(int)));
    xc=(int *)malloc((size_t) (k1*sizeof(int)));

    for(k2=1;k2<=k2max;k2++){
      fprintf(stderr, "k2(2)=%d\n", k2);
      rt = (int **)malloc((size_t) (k2*sizeof(int*)));
      for(i=0;i<k2;i++)
	rt[i]=(int *)malloc((size_t) (n*sizeof(int)));
      yc=(int *)malloc((size_t) (k2*sizeof(int)));
      firstcomb(xc,numrtall,k1);combflag=1;
      while(combflag==1){
	for(j=0;j<k1;j++)
	  veccp0(rtall[xc[j]],n,lt[j]);
	combflag1=firstcombfrom(yc,numrtall,k2,xc[0]+1);
	//printvec(yc,k2);fprintf(stderr, "combflag1=%d\n", combflag1);
	while(combflag1==1){
	  for(j=0;j<k2;j++)
	    veccp0(rtall[yc[j]],n,rt[j]);
	  //
	  if(disj(lt,rt,k1,k2,n) && !incliquelist(lt,k1,rt,k2,n,*leftlist,*leftsz,*rightlist,*rightsz,*numclqs) && !sqremleavesorphans(allt, cfs,r,n,lt,k1,rt,k2,1) && sqmidsneeded(allt,cfs,r,n,lt,k1,rt,k2,allmids)){//all midpoints needed?
	    fprintf(stderr, "here\n"); exit(0);	    
	    addtocliquelist(lt,k1,rt,k2,n,leftlist,leftsz,rightlist,rightsz,numclqs);
	    fprintf(stderr, "squares found: %d\n",*numclqs);
	    /* cout << "good pair:\n";  printvec1(lt[0],n);  printvec1(rt[0],n); cout << endl; */
	  }
	  combflag1=nextcomb(yc,numrtall,k2);
	}
	combflag=nextcomb(xc,numrtall,k1);
      }
      free((char*)yc);ifree(rt,k2);
    }
    free((char*)xc);ifree(lt,k1);
  }
  fprintf(stderr, "exiting\n");
  return;
}


long trimid(int **allt, int n, int *cfs, long r, long j1, long k1){
  int *vc;
  long s1;
  vc=midpt(allt[j1],allt[k1],cfs[j1],cfs[k1],n);
  if(vc){
    if((s1=binisearch(allt,r,vc,n))<0 || cfs[s1]<=0){// midpoint not found or negative
      free((char*)vc);
      return -1;
    }
    free((char*)vc);
    return s1;
  }

  return -1;
}

//check if needed mid-points are available

bool hasalltrimid(int **allt, int n, int *cfs, long r, int **neededby, long s){
  long j1,k1;
  int j,k,m;
  m=neededby[s][0];//number which need allt[s]
  for(j=2;j<m-1;j+=2){// go through opposite pairs
    j1=neededby[s][j];
    for(k=j+2;k<m+1;k+=2){
      k1=neededby[s][k];
      if(trimid(allt, n, cfs, r, j1, k1)<0)
	return 0;
    }
  }
  return 1;
}

void getalltrimids(int **allt, int n, int *cfs, long r, int **neededby, bool *trimids){
  long s;
  for(s=0;s<r;s++)
    trimids[s]=hasalltrimid(allt, n, cfs, r, neededby, s);
  return;
}

//Get the Real and imaginary part of det(J+iI)

void ReImJpiI(matrix J, int n, int numv, int ***explst1, int **cflst1, long *r1, int ***explst2, int **cflst2, long *r2, bool hom){
  ex v[numv+2];
  char str1[10];
  int i,k;
  ex tmp,tmp1=0,tmp2=0;
  int **lst;
  bool fl;
  long rtmp;

  //  matrix J(n,n);

  for(i=0;i<=numv+1;i++){
    sprintf(str1, "v%d", i);
    v[i]=get_possymbol(str1);
  }

  fl=1;
  for(k=n;k>=0;k-=2){// starts with det(J) and alternates sign
    if(!hom){tmp=getminorsum0(J,n,k);}else{tmp=expand(pow(v[numv+1],n-k)*getminorsum0(J,n,k));}
    if(fl){tmp1+=tmp;fl=0;}
    else{tmp1-=tmp;fl=1;}
  }
  cout << tmp1 << endl;

  //  cout << endl << tmp << endl;
  polytointlist0(tmp1, &lst, cflst1, r1);
  if(!hom){
    (*explst1)=moncflisttoexplist0(lst,*cflst1,*r1,numv,1);// redefines cflst
    rtmp=*r1;
    elongate(explst1,cflst1,r1,numv);
  }
  else{
    (*explst1)=moncflisttoexplist0(lst,*cflst1,*r1,numv+1,1);// redefines cflst
    rtmp=*r1;
    elongate(explst1,cflst1,r1,numv+1);
  }
  ifree_l(lst, rtmp);
  //cout << "sanity check = " <<  expand(tmp1-expltopolysimp(*explst1,*cflst1,*r1,numv)) << endl;

  fl=1;
  for(k=n-1;k>=0;k-=2){// starts with J_{n-1} and alternates sign
    if(!hom){tmp=getminorsum0(J,n,k);}else{tmp=expand(pow(v[numv+1],n-k)*getminorsum0(J,n,k));}
    if(fl){tmp2+=tmp;fl=0;}
    else{tmp2-=tmp;fl=1;}
  }
  cout << tmp2 << endl;

  //  cout << endl << tmp << endl;
  polytointlist0(tmp2, &lst, cflst2, r2);
  if(!hom){
    (*explst2)=moncflisttoexplist0(lst,*cflst2,*r2,numv,1);// redefines cflst
    rtmp=*r2;
    elongate(explst2,cflst2,r2,numv);
  }
  else{
    (*explst2)=moncflisttoexplist0(lst,*cflst2,*r2,numv+1,1);// redefines cflst
    rtmp=*r2;
    elongate(explst2,cflst2,r2,numv+1);
  }
  ifree_l(lst, rtmp);
  //cout << "sanity check = " <<  expand(tmp2-expltopolysimp(*explst2,*cflst2,*r2,numv)) << endl;

  return;
}

//homogeneous version

ex ex_J2pI_simphom(matrix J, int n, int numv){
  ex v[numv+2];
  char str1[10];
  int i,k;
  ex tmp,tmp1;

  //  matrix J(n,n);

  for(i=0;i<=numv+1;i++){
    sprintf(str1, "v%d", i);
    v[i]=get_possymbol(str1);
  }

  for(k=0;k<=n;k++){
    tmp1=expand(getminorsum(J,n,k));
    tmp+=expand(pow(v[numv+1], 2*(n-k))*tmp1);
  }
  return tmp;
}



ex GetJn(matrix J, int n, int numv, bool hom, ex Jn[]){
  ex v[numv+2];
  char str1[10];
  int i,k;
  ex tmp;

  //  matrix J(n,n);

  for(i=0;i<=numv+1;i++){
    sprintf(str1, "v%d", i);
    v[i]=get_possymbol(str1);
  }

  for(k=n;k>=0;k-=1){// starts with det(J) and alternates sign
    if(!hom){tmp=getminorsum0(J,n,k);}else{tmp=expand(pow(v[numv+1],n-k)*getminorsum0(J,n,k));}
    Jn[k]=tmp;
    cout << "J" << k << ": " << tmp << ";" << endl;
  }

  //tmp=factor(Jn[6]);
  //tmp=expand(Jn[4]*Jn[4] + Jn[6]*Jn[2]);
  //tmp=expand((Jn[5]-Jn[3]+Jn[0]) - Jn[1]*(Jn[4]-Jn[2]+Jn[0]));

  //cout << "Jn6: " << tmp << ";" << endl;

  tmp=ex_J2pI_simphom(J, n, numv);
  //tmp=tmp-tmp.subs(v[8]==0);
  //cout << "qqq: " << tmp << endl;
  //tmp=expand(pow(Jn[0]-Jn[2]+Jn[4],2)+pow(Jn[1]-Jn[3],2));
  //tmp=expand((Jn[0]-Jn[2]+Jn[4])*(Jn[0].subs(lst(v[2]==0,v[3]==0))-Jn[2].subs(lst(v[2]==0,v[3]==0))+Jn[4].subs(lst(v[2]==0,v[3]==0)))+(Jn[1]-Jn[3])*(Jn[1].subs(lst(v[2]==0,v[3]==0))-Jn[3].subs(lst(v[2]==0,v[3]==0))));

/*   tmp-=expand(pow(v[15]*v[1]*pow(v[17],2)*v[13]+v[15]*pow(v[17],2)*v[13]*v[5]-v[12]*v[1]*v[6]*v[13]*v[5]+pow(v[17],2)*v[13]*v[7]*v[5]+v[1]*v[11]*pow(v[17],2)*v[7]-v[10]*v[15]*v[1]*v[6]*v[5]+v[15]*v[1]*v[11]*pow(v[17],2)-v[15]*v[1]*v[6]*v[11]*v[5]+v[11]*pow(v[17],2)*v[7]*v[4]+v[10]*pow(v[17],2)*v[7]*v[4]+v[15]*v[11]*pow(v[17],2)*v[5]-v[10]*v[12]*v[1]*v[6]*v[5]+v[12]*pow(v[17],2)*v[13]*v[5]+v[10]*v[12]*v[1]*pow(v[17],2)+v[12]*v[1]*pow(v[17],2)*v[13]+v[1]*pow(v[17],2)*v[13]*v[7]+v[10]*v[12]*pow(v[17],2)*v[4]-v[12]*v[13]*v[2]*v[3]*v[5]-v[13]*v[2]*v[7]*v[3]*v[5]+v[10]*v[15]*pow(v[17],2)*v[5]-v[10]*v[12]*v[2]*v[3]*v[5]-v[15]*v[11]*v[2]*v[3]*v[5]+v[11]*pow(v[17],2)*v[7]*v[5]+v[10]*pow(v[17],2)*v[7]*v[5]+v[15]*pow(v[17],2)*v[13]*v[4]+pow(v[17],2)*v[13]*v[7]*v[4]-v[15]*v[1]*v[6]*v[13]*v[5]+v[10]*v[12]*pow(v[17],2)*v[5]-v[11]*v[2]*v[7]*v[3]*v[5]+v[10]*v[15]*pow(v[17],2)*v[4]-v[10]*v[2]*v[7]*v[3]*v[5]-v[15]*v[13]*v[2]*v[3]*v[5]-v[10]*v[15]*v[2]*v[3]*v[5]+v[10]*v[1]*pow(v[17],2)*v[7]+v[15]*v[11]*pow(v[17],2)*v[4]+v[12]*pow(v[17],2)*v[13]*v[4]+v[10]*v[15]*v[1]*pow(v[17],2),2)*pow(v[17],8)); */


/* tmp-=expand(pow(v[6]*v[11]*pow(v[17],2)*v[3]+v[8]*pow(v[17],2)*v[14]*v[3]-v[12]*v[8]*v[4]*v[14]*v[10]-v[12]*v[6]*v[14]*v[3]*v[10]+pow(v[17],2)*v[2]*v[14]*v[3]-v[12]*v[8]*v[14]*v[3]*v[10]+v[6]*v[11]*pow(v[17],2)*v[4]+v[1]*v[6]*pow(v[17],2)*v[14]-v[12]*v[6]*v[4]*v[14]*v[10]+v[1]*v[8]*pow(v[17],2)*v[14]-v[11]*v[2]*v[7]*v[14]*v[3]+v[6]*pow(v[17],2)*v[14]*v[3]-v[11]*v[8]*v[7]*v[4]*v[14]-v[1]*v[11]*v[8]*v[7]*v[14]+v[1]*v[8]*pow(v[17],2)*v[13]+v[1]*v[11]*v[8]*pow(v[17],2)-v[12]*v[1]*v[6]*v[14]*v[10]+v[8]*pow(v[17],2)*v[13]*v[4]+v[11]*v[8]*pow(v[17],2)*v[4]-v[12]*v[2]*v[4]*v[14]*v[10]+v[6]*pow(v[17],2)*v[13]*v[3]+v[11]*pow(v[17],2)*v[2]*v[4]-v[12]*v[1]*v[8]*v[14]*v[10]+v[8]*pow(v[17],2)*v[4]*v[14]+pow(v[17],2)*v[13]*v[2]*v[3]+v[1]*v[6]*pow(v[17],2)*v[13]+v[11]*v[8]*pow(v[17],2)*v[3]-v[11]*v[2]*v[7]*v[4]*v[14]+v[6]*pow(v[17],2)*v[4]*v[14]+v[8]*pow(v[17],2)*v[13]*v[3]-v[11]*v[8]*v[7]*v[14]*v[3]+pow(v[17],2)*v[2]*v[4]*v[14]-v[12]*v[2]*v[14]*v[3]*v[10]+pow(v[17],2)*v[13]*v[2]*v[4]+v[1]*v[6]*v[11]*pow(v[17],2)+v[11]*pow(v[17],2)*v[2]*v[3]+v[6]*pow(v[17],2)*v[13]*v[4],2)*pow(v[17],8)); */




/*  tmp-=expand(2*pow(v[17],8)*v[13]*pow(v[12]*pow(v[17],2)*v[4]+v[12]*v[1]*pow(v[17],2)-v[2]*v[3]*v[5]*v[15]+pow(v[17],2)*v[7]*v[5]-v[12]*v[2]*v[3]*v[5]-v[1]*v[6]*v[5]*v[15]-v[2]*v[7]*v[3]*v[5]+pow(v[17],2)*v[5]*v[15]+v[12]*pow(v[17],2)*v[5]-v[12]*v[1]*v[6]*v[5]+v[1]*pow(v[17],2)*v[7]+v[1]*pow(v[17],2)*v[15]+pow(v[17],2)*v[4]*v[15]+pow(v[17],2)*v[7]*v[4],2)*v[14]); */

/*  tmp-=expand(pow(v[5],2)*pow(v[17],8)*pow(v[10]*v[12]*v[2]*v[14]-v[8]*pow(v[17],2)*v[13]+v[10]*v[12]*v[6]*v[14]+v[11]*v[2]*v[7]*v[14]-v[6]*v[11]*pow(v[17],2)+v[11]*v[8]*v[7]*v[14]-v[6]*pow(v[17],2)*v[13]-v[11]*v[8]*pow(v[17],2)-v[6]*pow(v[17],2)*v[14]+v[10]*v[12]*v[8]*v[14]-pow(v[17],2)*v[13]*v[2]-pow(v[17],2)*v[2]*v[14]-v[11]*pow(v[17],2)*v[2]-v[8]*pow(v[17],2)*v[14],2)); */

/*  tmp-=expand(pow(v[5]*v[2]*v[7]*v[3]+v[5]*v[15]*v[2]*v[3]-v[5]*v[12]*pow(v[17],2)-v[12]*pow(v[17],2)*v[4]+v[5]*v[12]*v[1]*v[6]-v[12]*v[1]*pow(v[17],2)-v[5]*v[15]*pow(v[17],2)-pow(v[17],2)*v[7]*v[4]+v[5]*v[12]*v[2]*v[3]-v[15]*v[1]*pow(v[17],2)+v[5]*v[15]*v[1]*v[6]-v[5]*pow(v[17],2)*v[7]-v[15]*pow(v[17],2)*v[4]-v[1]*pow(v[17],2)*v[7],2)*pow(v[17],8)*pow(v[14],2)); */

/*  tmp-=expand(pow(v[17],8)*pow(v[10]*v[12]*pow(v[17],2)*v[2]+v[15]*v[11]*pow(v[17],2)*v[2]+v[15]*pow(v[17],2)*v[13]*v[2]-v[15]*v[1]*v[6]*v[13]*v[9]+pow(v[17],2)*v[13]*v[2]*v[7]+v[10]*v[15]*pow(v[17],2)*v[2]+v[12]*pow(v[17],2)*v[13]*v[2]+v[10]*pow(v[17],2)*v[2]*v[7]+v[11]*pow(v[17],2)*v[2]*v[7]-v[15]*v[1]*v[6]*v[11]*v[9]-v[12]*v[1]*v[6]*v[13]*v[9]-v[10]*v[15]*v[1]*v[6]*v[9],2)); */

/*  tmp-=expand(pow(v[17],8)*v[4]*v[5]*pow(v[2]*v[7]*v[14]*v[11]+v[8]*v[14]*v[10]*v[12]-pow(v[17],2)*v[13]*v[6]-v[8]*pow(v[17],2)*v[13]+v[2]*v[14]*v[10]*v[12]-pow(v[17],2)*v[13]*v[2]+v[8]*v[7]*v[14]*v[11]+v[14]*v[10]*v[12]*v[6]-pow(v[17],2)*v[14]*v[6]-pow(v[17],2)*v[2]*v[11]-v[8]*pow(v[17],2)*v[14]-pow(v[17],2)*v[6]*v[11]-pow(v[17],2)*v[2]*v[14]-v[8]*pow(v[17],2)*v[11],2)); */

/*  tmp-=expand(pow(v[17],8)*pow(v[8]*pow(v[17],2)*v[3]*v[12]+pow(v[17],2)*v[2]*v[3]*v[12]-v[8]*v[7]*v[4]*v[16]*v[11]+pow(v[17],2)*v[12]*v[1]*v[6]+v[8]*pow(v[17],2)*v[12]*v[1]+v[8]*pow(v[17],2)*v[4]*v[12]-v[8]*v[7]*v[16]*v[1]*v[11]-v[2]*v[7]*v[4]*v[16]*v[11]+pow(v[17],2)*v[4]*v[12]*v[6]-v[8]*v[7]*v[3]*v[16]*v[11]+pow(v[17],2)*v[3]*v[12]*v[6]+pow(v[17],2)*v[2]*v[4]*v[12],2)); */

/*  tmp-=expand(pow(v[17],8)*pow(pow(v[17],2)*v[13]*v[2]+v[8]*pow(v[17],2)*v[14]+v[11]*pow(v[17],2)*v[2]-v[11]*v[8]*v[7]*v[14]-v[8]*v[14]*v[10]*v[12]+v[11]*pow(v[17],2)*v[6]+pow(v[17],2)*v[2]*v[14]+v[8]*pow(v[17],2)*v[13]-v[14]*v[10]*v[12]*v[6]+pow(v[17],2)*v[13]*v[6]-v[11]*v[2]*v[7]*v[14]+v[11]*v[8]*pow(v[17],2)+pow(v[17],2)*v[14]*v[6]-v[2]*v[14]*v[10]*v[12],2)*v[4]*v[5]); */

/*  // constructed manually */

/*  tmp-=expand(pow(v[17],8)*pow(v[11]*pow(v[17],2)*v[7]*v[16]+v[10]*v[12]*v[6]*v[9]*v[16]-v[7]*v[9]*v[10]*v[14]*v[16]-v[7]*v[8]*v[10]*v[14]*v[16]-v[2]*v[7]*v[10]*v[14]*v[16]-v[12]*pow(v[17],4)-pow(v[17],4)*v[16]-v[15]*pow(v[17],4),2)); */

 // a subset of J6*J2


 /* tmp-=expand(2*Jn[6]*pow(v[17],7)*(v[2]*v[10]+v[2]*v[12]+v[1]*v[12]+v[2]*v[11]+v[3]*v[12]+v[2]*v[14]+v[2]*v[13]+v[2]*v[3]+v[1]*v[9]+v[1]*v[8]+v[3]*v[8]+v[3]*v[9]+v[2]*v[15]+v[2]*v[16]+v[1]*v[15]+v[1]*v[16]+v[3]*v[16]+v[3]*v[15]+(v[4]+v[5])*(v[15]+v[16]+v[12])+(v[8]+v[9])*(v[11]+v[12]+v[13]+v[14]+v[15]+v[16])+v[10]*(v[12]+v[15]+v[16])+v[11]*(v[15]+v[16]))); */

 //tmp=expand(Jn[5]*Jn[5]);


  //tmp=expand(Jn[2]*Jn[2]+2*Jn[0]*Jn[4]-2*Jn[1]*Jn[3]+Jn[1]*Jn[1]+Jn[3]*Jn[3]- 2*Jn[0]*Jn[2] + 2*Jn[1]*Jn[5] + Jn[5]*Jn[5]);

  //tmp=expand(Jn[1]*Jn[1]+ 2*Jn[1]*Jn[5] + Jn[5]*Jn[5]);

  /* tmp-=expand(pow(v[17],12)*pow(pow(v[17],2)*v[1]-v[2]*v[3]*v[5]+pow(v[17],2)*v[4]+pow(v[17],2)*v[5]-v[6]*v[5]*v[1],2) + pow(v[17],12)*pow(v[14]*v[10]*v[12]+v[11]*v[7]*v[14]-pow(v[17],2)*v[13]-v[11]*pow(v[17],2)-pow(v[17],2)*v[14],2) + pow(pow(v[17],2)*v[12]-v[11]*v[7]*v[16],2)*pow(v[17],12) + pow(v[17],12)*pow(pow(v[17],2)*v[2]-v[6]*v[9]*v[1],2) + pow(pow(v[17],2)*v[9]-v[6]*v[14]*v[10]-v[6]*v[10]*v[12],2)*pow(v[17],12) + pow(pow(v[17],2)*v[7]-v[2]*v[3]*v[16]-v[3]*v[16]*v[5],2)*pow(v[17],12) + pow(v[17],14)*pow(v[11]*v[7]-v[14]*v[10],2) + pow(v[17],14)*pow(v[7]*v[14]-v[11]*v[16],2) + pow(v[6]*v[5]-pow(v[17],2),2)*pow(v[17],12)*v[3]*v[1] + pow(v[6]*v[1]-v[9]*v[5],2)*pow(v[17],14)); */

 //tmp=expand(Jn[5]*Jn[5]-2*Jn[6]*Jn[4] + 2*v[3]*pow(v[5],2)*v[6]*v[7]*v[10]*v[11]*v[12]*v[14]*v[16]*pow(v[17],8) + 2*v[3]*v[4]*v[5]*v[6]*v[7]*v[10]*v[11]*v[12]*v[14]*v[16]*pow(v[17],8) + 2*pow(v[3],2)*v[5]*v[6]*v[7]*v[10]*v[11]*v[12]*v[14]*v[16]*pow(v[17],8) + 2*pow(v[3],2)*pow(v[5],2)*pow(v[9],2)*v[10]*v[11]*v[12]*v[14]*pow(v[17],8) + 4*pow(v[3],2)*pow(v[5],2)*v[8]*v[9]*v[10]*v[11]*v[12]*v[14]*pow(v[17],8) + 2*pow(v[3],2)*pow(v[5],2)*pow(v[8],2)*v[10]*v[11]*v[12]*v[14]*pow(v[17],8));

 //v5

  return tmp;
}

ex J2pI_simp(matrix J, int n, int numv, int ***explst, int **cflst, long *r){
  ex v[numv+2];
  char str1[10];
  int i,k;
  ex tmp,tmp1,tmp2;
  int **lst;

  //  matrix J(n,n);

  for(i=0;i<=numv+1;i++){
    sprintf(str1, "v%d", i);
    v[i]=get_possymbol(str1);
  }

  for(k=0;k<=n;k++)
    tmp+=expand(getminorsum(J,n,k));

  polytointlist0(tmp, &lst, cflst, r);
  (*explst)=moncflisttoexplist0(lst,*cflst,*r,numv,1);// redefines cflst
  //elongate(explst1,cflst1,r1,numv);
  ifree_l(lst, *r);
  return tmp;
}

//homogeneous version

ex J2pI_simphom(matrix J, int n, int numv, int ***explst, int **cflst, long *r){
  ex v[numv+2];
  char str1[10];
  int i,k;
  ex tmp,tmp1;
  int **lst;

  //  matrix J(n,n);

  for(i=0;i<=numv+1;i++){
    sprintf(str1, "v%d", i);
    v[i]=get_possymbol(str1);
  }

  for(k=0;k<=n;k++){
    tmp1=expand(getminorsum(J,n,k));
    tmp+=expand(pow(v[numv+1], 2*(n-k))*tmp1);
  }

  polytointlist0(tmp, &lst, cflst, r);
  (*explst)=moncflisttoexplist0(lst,*cflst,*r,numv+1,1);// redefines cflst
  //elongate(explst1,cflst1,r1,numv);
  ifree_l(lst, *r);
  return tmp;
}





//
// Compute the determinant of J^2 + omega^2 I
// for matrix J
//

void J2pI(matrix J, int n, int numv){
  ex v[numv+2];
  char str1[10];
  int i,k,dim;
  ex tmp,tmp1,tmp2;
  int **lst,**explst,*cflst;
  long r,s,s1,t;
  float ns;
  int *needcount,*needexcl;
  int *vc;

  //  matrix J(n,n);

  for(i=0;i<=numv+1;i++){
    sprintf(str1, "v%d", i);
    v[i]=get_possymbol(str1);
  }


  for(k=0;k<=n;k++){
    tmp2=getminorsum(J,n,k);
    //tmp2=getminorsum(J,n,k);
    //cout << endl <<endl << tmp2 << endl << endl;
    tmp1=expand(tmp2);
    //cout << "k=" << k << ":" << tmp2 << endl;
    if(tmp1!=0 && !(tmp1.info(info_flags::positive))&& !(tmp2.info(info_flags::positive))){
      fprintf(stderr, "could be negative: k=%d\n", k);
      //cout << tmp1 << endl;
    }
    tmp+=expand(pow(v[numv+1], 2*(n-k))*tmp1);
  }
  //  cout << endl << tmp << endl;
  dim=polytointlist(tmp, &lst, &cflst, &r);//dim=order or monomials
  explst=moncflisttoexplist(lst,dim,cflst,r,numv+1,1);// redefines cflst
  listprint(cflst, explst, numv+1, r);
  //treeprint3(cflst,lst,dim,r,numv+1,0);
  //exit(0);
  needcount = (int *)malloc((size_t) (r*sizeof(int))); 
  needexcl = (int *)malloc((size_t) (r*sizeof(int))); 
  for(s=0;s<r;s++){
    needcount[s]=0;
    needexcl[s]=0;
  }

  int **needs=(int **) malloc((size_t)(r*sizeof(int*)));
  int **neededby=(int **) malloc((size_t)(r*sizeof(int*)));
  for(s=0;s<r;s++){
    neededby[s]=(int *) malloc((size_t)(sizeof(int)));
    neededby[s][0]=0;
    needs[s] = (int*) malloc(sizeof(int));
    needs[s][0]=0; //number needed
  }

  for(s=0;s<r;s++){
    if(cflst[s]<0){
      ns=numsplits(s,numv+1,explst,cflst,r,&needcount,&needexcl,needs,neededby,0);
      cout << "numsplits[" <<  s << "]: " << ns << endl;
    }
  }

  cout << "can be used by\n" << endl;
  for(s=0;s<r;s++){
    if(cflst[s]>0){
      cout << s << ": " << cflst[s] << ", "; printvec1(explst[s],numv+1);
      printvec1(neededby[s],neededby[s][0]+1);
      cout << endl;
    }
  }


  cout << "simple squares\n" << endl;
  t=0;
  ex tmp3;
  long inds[3];
  cout << tmp << endl;
  for(s=0;s<r;s++){
    if(neededby[s][0]==2 && neededby[s][2]>s && neededby[neededby[s][2]][0]==2 && (splt(cflst[s],cflst[neededby[s][2]]) >=abs(cflst[neededby[s][1]]))){ // both halves needed exactly once, powers add up
      t++;
      cout << t << endl;
      cout << "   " << cflst[neededby[s][1]] << ", "; printvec1(explst[neededby[s][1]],numv+1);
      cout << "   " << cflst[s] << ", "; printvec1(explst[s],numv+1);
      cout << "   " << cflst[neededby[s][2]] << ", "; printvec1(explst[neededby[s][2]],numv+1);
      inds[0]=s;inds[1]=neededby[s][1];inds[2]=neededby[s][2];
      tmp3=expltopoly(explst, cflst, inds, 3, numv+1);
      cout << endl << tmp3 << endl;

      tmp-=tmp3;
    }

  }
  cout << endl << tmp;
  exit(0);


  bool *trimids=(bool *) malloc((size_t)(r*sizeof(bool)));
  getalltrimids(explst, numv+1, cflst, r, neededby,trimids);
  cout << "nontrivial trimids\n" << endl;
  for(s=0;s<r;s++){
    if(neededby[s][0] > 2 && trimids[s]){// nontrivial
      printvec1(neededby[s],neededby[s][0]+1);
      cout << endl;
    }
  }
  
 


  cout << "simple triangles\n" << endl;
  t=0;
  for(s=0;s<r;s++){
    if(neededby[s][0]==4 && neededby[neededby[s][2]][0]==2 && neededby[neededby[s][4]][0]==2 && (splt(cflst[s],cflst[neededby[s][2]]) >=abs(cflst[neededby[s][1]])) && (splt(cflst[s],cflst[neededby[s][4]]) >=abs(cflst[neededby[s][3]])) && trimids[s]){

      vc=midpt(explst[neededby[s][2]],explst[neededby[s][4]],cflst[neededby[s][2]],cflst[neededby[s][4]],numv+1);
      s1=binisearch(explst,r,vc,n);
      if(vc)
    	free((char*)vc);
      if(splt(cflst[neededby[s][4]],cflst[neededby[s][4]]) <=abs(cflst[s1])){

    	t++;
    	cout << t << endl;
    	cout << "   " << cflst[neededby[s][1]] << ", "; printvec1(explst[neededby[s][1]],numv+1);
    	cout << "   " << cflst[neededby[s][3]] << ", "; printvec1(explst[neededby[s][3]],numv+1);
    	cout << "   " << cflst[neededby[s][2]] << ", "; printvec1(explst[neededby[s][2]],numv+1);
    	cout << "   " << cflst[neededby[s][4]] << ", "; printvec1(explst[neededby[s][4]],numv+1);
    	cout << "   " << cflst[s1] << ", "; printvec1(explst[s1],numv+1);

      }
    }
  }


  
  //treeprint3(cflst,lst,dim,r,numv+1,0);
  
  ifree_l(needs, r);
  ifree_l(neededby, r);
  ifree_l(lst, r);
  ifree_l(explst, r);
  free((char *)cflst); free((char *)needcount);free((char *)needexcl);free((char*)trimids);
  return;
}


void printblock(int **block, int numinblock, int n){
  int i;
  for(i=0;i<numinblock;i++){
    cout << "   "; printvec1(block[i],n);
  }
  cout << endl;
  return;
}

void printcfblock(int **block, int *cfs, int numinblock, int n){
  int i;
  for(i=0;i<numinblock;i++){
    printf("%d,  ", cfs[i]); printvec1(block[i],n);
  }
  printf("\n");
  return;
}

void printblock1(int **block, int *cfs, int numinblock, int n){
  int i;
  for(i=0;i<numinblock;i++){
    fprintf(stderr, "%d,  ", cfs[i]); printvec(block[i],n);
  }
  fprintf(stderr, "\n");
  return;
}

//
// get the indices of a set of monomials in a list from a larger list
// Assume memory has been allocated for out
//

void getsublist(int **biglist, long r, int **smalllist, int sz, int n, long *out){
  int i;
  for(i=0;i<sz;i++)
    out[i]=binisearch(biglist,r,smalllist[i],n);
  return;
}


long **blockstoinds(int **biglist, long r, int ***blocks, int *numinblock, int numblocks, int n){
  int k;
  long **out;
  out=(long **)malloc((size_t) (numblocks*sizeof(long*)));
  for(k=0;k<numblocks;k++){
    out[k]=(long *)malloc((size_t) (numinblock[k]*sizeof(long))); 
    getsublist(biglist, r, blocks[k], numinblock[k], n, out[k]);
  }
  return out;
}

/* int triangle(int n){ */
/*   if(n<=0 || n > 10000){ */
/*     cout << "ERROR" << endl; */
/*     return -1; */
/*   } */
/*   if(n==1) */
/*     return 1; */
/*   return n+triangle(n-1); */
/* } */

/* int fibonacci(int n){ */
/*   if(n<1 || n>100000){ */
/*     cout << "ERROR" << endl; */
/*     exit(0); */
/*   } */
/*   if (n==1) */
/*     return 1; */
/*   return fibonacci(n-1)+fibonacci(n-2); */
/* } */

/* ex power111(ex n){ */
/*   if(n<2 || n>10){ */
/*     cout << "ERROR" << endl; */
/*     return -1; */
/*   } */
/*   if(n==2) */
/*     return 2; */
/*   return pow(power111(n-1),n); */
/* } */

/* long factorial1(int n){ */
/*   if (n<0 || n>1000){ */
/*     cout << "ERROR" << endl; */
/*     return -1; */
/*   } */
/*   if (n==0) */
/*     return 1; */
/*   return n*factorial1(n-1); */
/* } */

ex lststobrac(int **explst, long *lst1, int num1, long *lst2, int num2, int numv){
  int i;
  ex extot;
  extot=0;
  for(i=0;i<num1;i++)
    extot+=expltomono(1,explst[lst1[i]],numv);
  for(i=0;i<num2;i++)
    extot-=expltomono(1,explst[lst2[i]],numv);
  return extot;
}

void forcedclq(long indx, long r, int *needexcl, int **needs, int **neededby, int *used, int *tot){
  long m;
  int j,k,flg;
  if(!needexcl[indx]) // not forced
    return;
  if(used[indx]==-1)// already dealt with
    return;
  // flg=1 means odds on left; flg=2 means odds on right
  used[indx]=-1;
  flg=0;
  for(j=1;j<=needs[indx][0];j++){// step through pairs needed by indx
    if(used[needs[indx][j]]){
      if(j%2==used[needs[indx][j]]%2){
	flg=1;break;
      }
      else{
	flg=2;break;
      }
    }
  }
  // set used flags
  if(!flg || flg==1){ //odd,even...
    for(j=1;j<=needs[indx][0];j++){
      if(used[needs[indx][j]] && used[needs[indx][j]]%2!=j%2){
	fprintf(stderr, "impossible situation(1).\n");/* exit(0); */used[indx]=-1;(*tot)=0;return;
      }
      else{
	if(j%2==0 && !used[needs[indx][j]]){
	  used[needs[indx][j]]=2;(*tot)++;
	}
	else if(!used[needs[indx][j]]){
	  used[needs[indx][j]]=1;(*tot)++;
	}
      }
    }
  }
  else if(flg==2){ //even,odd...
    for(j=1;j<=needs[indx][0];j++){
      if(used[needs[indx][j]] && used[needs[indx][j]]!=j%2+1){
	fprintf(stderr, "impossible situation(2).\n");/* exit(0); */used[indx]=-1;(*tot)=0;return;
      }
      else{
	if(j%2==0 && !used[needs[indx][j]]){
	  used[needs[indx][j]]=1;(*tot)++;
	}
	else if(!used[needs[indx][j]]){
	  used[needs[indx][j]]=2;(*tot)++;
	}
      }
    }
  }

  for(j=1;j<=needs[indx][0];j++){//step through positive terms surrounding indx
    for(k=1;k<=neededby[needs[indx][j]][0];k+=2){// step through negative terms dominated by each positive term
      m=neededby[needs[indx][j]][k];
      forcedclq(m,r,needexcl,needs,neededby,used,tot);
    }
  }
  return;
}

int isfreehonest(int **explst, int *cflst, int numv, long r, int **explst1, int *cflst1, int r1, int **explst2, int *cflst2, int r2, int *used, int *used1, int *used2, int *used1tmp, int *used2tmp){
  long s,s1;
  int *vc;
  int flg=0, sgn=0;
  //fprintf(stderr, "entering isfreehonest\n");
  for(s=0;s<r;s++){
    if(used[s]>0 && !iseven(explst[s],numv)){
      //fprintf(stderr, "couldn't be halved.\n");
      return 0;
    }
  }
  for(s=0;s<r;s++){
    if(used[s]>0){// left or right in clique
      vc=halve(explst[s], numv);
      if((s1=binisearch(explst1,r1,vc,numv))>=0){
	if(flg==2){
	  //fprintf(stderr, "mixed(1).\n");
	  return 0;
	}
	if(used1[s1]){// not available?
	  //fprintf(stderr, "failing s1 = %ld\n", s1);
	  return 0;
	}
	else{
	  //fprintf(stderr, "s1 = %ld\n", s1);
	  used1tmp[s1]=1;
	}
	flg=1;
	if(!sgn){
	  if(used[s]==1)
	    sgn=cflst1[s1];
	  else
	    sgn=-cflst1[s1];
	}
	if((used[s]==1 && sgn*cflst1[s1]<0) || (used[s]==2 && sgn*cflst1[s1]>0)){
	  //fprintf(stderr, "sign change(1).\n");
	  return 0;
	}
      }
      else if((s1=binisearch(explst2,r2,vc,numv))>=0){
	if(flg==1){
	  //fprintf(stderr, "mixed(2).\n");
	  return 0;
	}
	if(used2[s1]){//available?
	  //fprintf(stderr, "failing s1 = %ld\n", s1);
	  return 0;
	}
	else{
	  //fprintf(stderr, "s1 = %ld\n", s1);
	  used2tmp[s1]=1;
	}
	flg=2;
	if(!sgn){
	  if(used[s]==1)
	    sgn=cflst2[s1];
	  else
	    sgn=-cflst2[s1];
	}
	if((used[s]==1 && sgn*cflst2[s1]<0) || (used[s]==2 && sgn*cflst2[s1]>0)){
	  //fprintf(stderr, "sign change(2).\n");
	  return 0;
	}

      }
      else // not found
	return 0;
    }
  }
  //fprintf(stderr, "returning %d\n", flg);
  return flg;
}


//
//assume allocation already done
//

void lstcp(int **explst,int *cflst,long r,int numv,int **explst1, int *cflst1){
  long s;
  int k;
  for(s=0;s<r;s++){
    cflst1[s]=cflst[s];
    for(k=0;k<numv;k++)
      explst1[s][k]=explst[s][k];
  }
  return;
}



// alters cflst

bool removelst(int **biglist, int *cflst, long r, int **smalllist, int *cfs, int sz, int n, bool rmall){
  int i;
  long out[sz];
  bool flg=1;
  for(i=0;i<sz;i++){
    out[i]=binisearch(biglist,r,smalllist[i],n);
    if(out[i]<0 || (cfs[i]>0 && (cflst[out[i]]<0 || cfs[i]>cflst[out[i]]))){
      if(rmall){
	fprintf(stderr, "WARNING: Unable to remove all squares\n");
	return 0;
      }
      else{
	flg=0;break;
      }
    }
  }
  if(flg){
    for(i=0;i<sz;i++){
      cflst[out[i]]-=cfs[i];
      cout << cfs[i] << ", "; printvec1(biglist[out[i]],n);
    }
    return 1;
  }
  return 0;
}



// remove a square
// alters cflst

bool rmsq(int **l, int *cfs, int *changelst, int r, long *leftinds, int numleft, long *rightinds, int numright, int n, int maxnegcf, ex *rmt, /* int ***newnegmid, int *numnnegmid,  */bool simp){
  int i,j,k,k2;
  int minfac=1000;
  int *vc;
  long s1;
  int *kleft=(int *)malloc((size_t) (numleft*sizeof(int)));
  int *kright=(int *)malloc((size_t) (numright*sizeof(int)));
  ex trm=0;

  if(maxnegcf==2){// simplest case
    minfac=1;
    for(i=0;i<numleft;i++)
      kleft[i]=1;
    for(i=0;i<numright;i++)
      kright[i]=1;
  }
  else if(simp){
    minfac=1;
    for(i=0;i<numleft;i++)
      kleft[i]=sqrt(cfs[leftinds[i]]);
    for(i=0;i<numright;i++)
      kright[i]=sqrt(cfs[rightinds[i]]);//integer square root
  }
  else{
    // get the smallest common factor
    for(i=0;i<numleft;i++){
      k=sqrt(cfs[leftinds[i]]);//integer square root
      kleft[i]=k;
      minfac=min(minfac,cfs[leftinds[i]]/(k*k));//integer division
    }
    for(i=0;i<numright;i++){
      k=sqrt(cfs[rightinds[i]]);//integer square root
      kright[i]=k;
      minfac=min(minfac,cfs[rightinds[i]]/(k*k));//integer division
    }
  }

  for(i=0;i<numleft;i++){//each monomial on left
    k2=minfac*kleft[i]*kleft[i];
    cfs[leftinds[i]]-=k2;changelst[leftinds[i]]+=k2;
    trm=expltomono(k2,l[leftinds[i]],n);(*rmt)+=trm;
  }
  for(i=0;i<numright;i++){//each monomial on right
    k2=minfac*kright[i]*kright[i];
    cfs[rightinds[i]]-=k2;changelst[rightinds[i]]+=k2;
    trm=expltomono(k2,l[rightinds[i]],n);(*rmt)+=trm;
  }

  // left midpoints (positive)

  for(i=0;i<numleft;i++){
    for(j=i+1;j<numleft;j++){
      vc=midpt(l[leftinds[i]],l[leftinds[j]],cfs[leftinds[i]], cfs[leftinds[j]],n);
      if(vc && (s1=binisearch(l,r,vc,n))>=0){
	k2=2*kleft[i]*kleft[j]*minfac;
	trm=expltomono(k2,vc,n);(*rmt)+=trm;
	cfs[s1]-=k2;changelst[s1]+=k2;
	free((char*)vc);
      }
      else{
	fprintf(stderr, "The gathering process failed.\n");
	exit(0);
      }
    }
  }

  // right midpoints (positive)

  for(i=0;i<numright;i++){
    for(j=i+1;j<numright;j++){
      vc=midpt(l[rightinds[i]],l[rightinds[j]],cfs[rightinds[i]], cfs[rightinds[j]],n);
      if(vc && (s1=binisearch(l,r,vc,n))>=0){
	k2=2*kright[i]*kright[j]*minfac;
	trm=expltomono(k2,vc,n);(*rmt)+=trm;
	cfs[s1]-=k2;changelst[s1]+=k2;
	free((char*)vc);
      }
      else{//failure
	fprintf(stderr, "The gathering process failed.\n");
	exit(0);
      }
    }
  }

 // mixed midpoints (negative)

  for(i=0;i<numleft;i++){
    for(j=0;j<numright;j++){
      vc=midpt(l[leftinds[i]],l[rightinds[j]],cfs[leftinds[i]], cfs[rightinds[j]],n);
      k2=2*kleft[i]*kright[j]*minfac;
      if(vc){
	s1=binisearch(l,r,vc,n);
	/* if(s1<0){ */
	/*   addnewtoarray((*newnegmid),*numnegmid,vc, n); */
	/*   fprintf(stderr, "WARNING: not a perfect square.\n"); */
	/*   (*numnegmid)++; */
	/* } */
	if(s1>=0){
	  cfs[s1]+=k2;changelst[s1]-=k2;
	}
	trm=expltomono(k2,vc,n);(*rmt)-=trm;
	free((char*)vc);
      }
      else
	return 0;
    }
  }

  free((char*)kleft);free((char*)kright);
  return 1;

}

// remove a perfect square which we know to be in the list
// all coefficients are 1. The square is provided in unexpanded form

bool rmsq1(int **l, int *cfs, int r, int n, int **lt, int numleft, int **rt, int numright, ex *rmt){
  int i,j;
  int *vc;
  long s1;
  ex trm;
  for(i=0;i<numleft;i++){//each monomial on left
    for(j=0;j<numleft;j++){
      vc=sum(lt[i],lt[j],n);
      s1=binisearch(l,r,vc,n);
      if(s1<0){fprintf(stderr, "not a true square 1\n");exit(0);}
      if(cfs[s1]<=0){ free((char*)vc);return 0;}
      cfs[s1]--;
      free((char*)vc);
    }
  }
  for(i=0;i<numright;i++){//each monomial on left
    for(j=0;j<numright;j++){
      vc=sum(rt[i],rt[j],n);
      s1=binisearch(l,r,vc,n);
      if(s1<0){fprintf(stderr, "not a true square 2\n");exit(0);}
      if(cfs[s1]<=0){free((char*)vc);return 0;}
      cfs[s1]--;
      free((char*)vc);
    }
  }

  for(i=0;i<numleft;i++){//each monomial on left
    for(j=0;j<numright;j++){
      vc=sum(lt[i],rt[j],n);
      s1=binisearch(l,r,vc,n);
      if(s1>=0)//{fprintf(stderr, "not a true square 3\n");exit(0);}
	cfs[s1]+=2;
      free((char*)vc);
    }
  }

  for(i=0;i<numleft;i++)
    trm+=expltomono(1,lt[i],n);
  for(i=0;i<numright;i++)
    trm-=expltomono(1,rt[i],n);

  (*rmt)=pow(trm,2);

  return 1;
}




int **cfchangetoblock(int **explst, int *cflst,long r,int *newcflst, int *blocksz, int **blockcfs, int numv){
  long k;
  int i,cf;
  int **expsout=NULL;  
  (*blocksz)=0;
  for(k=0;k<r;k++){
    if((cf=(newcflst[k]-cflst[k]))){// change
      if((*blocksz)==0){ // first entry
	expsout = (int**) malloc(sizeof(int*) * 1);
	(*blockcfs) = (int*) malloc(sizeof(int) * 1);
      }
      else{
	expsout=(int**) realloc(expsout, sizeof(int*) *((*blocksz)+1));
	(*blockcfs) = (int*) realloc((*blockcfs), sizeof(int) * ((*blocksz)+1));
      }
      expsout[(*blocksz)]=(int *)malloc((size_t) ((numv)*sizeof(int)));
      (*blockcfs)[(*blocksz)]=cf;
      for(i=0;i<numv;i++)
	expsout[(*blocksz)][i]=explst[k][i];
      (*blocksz)++;
    }
  }
  return expsout;


}

// recursive routine
//


bool sqrem(int **explst, int *cflst, long r, int numv, int *level, int *changelst, int *totiter, int ****sqout, int ***sqcfs, int **sqoutlen, int *numsqout){
  int j;
  long s;
  float ns;
  int *needcount,*needexcl;
  int **needs, **neededby;
  int negcount=0;
  long **lt=NULL,**rt=NULL;
  int *numlt=NULL,*numrt=NULL,*absnegmax=NULL;
  int cflst1[r];
  int changelst1[r];
  int totcliques;
  int maxsz=500;
  int **newblock;
  int blocksz;
  int *blockcfs;
  ex tmp2;

  if(*level>=100)
    return 0;

  (*level)++;(*totiter)++;

  needcount = (int *)malloc((size_t) (r*sizeof(int))); 
  needexcl = (int *)malloc((size_t) (r*sizeof(int))); 
  for(s=0;s<r;s++){
    if(cflst[s]<0){negcount++;}
    needcount[s]=0;needexcl[s]=0;
  }
  fprintf(stderr, "level = %d, negative count: %d\n", *level, negcount);
  //cout << negcount << " negative terms\n";
  //listprint(cflst, explst, numv, r);
  //cout << endl;

  if(!negcount){// what to output?
    (*level)--;free((char *)needcount);free((char *)needexcl);
    return 1;
  }

  needs=(int **) malloc((size_t)(r*sizeof(int*)));
  neededby=(int **) malloc((size_t)(r*sizeof(int*)));
  for(s=0;s<r;s++){
    neededby[s]=(int *) malloc((size_t)(sizeof(int)));
    neededby[s][0]=0;
    needs[s] = (int*) malloc(sizeof(int));
    needs[s][0]=0; //number needed
  }
  for(s=0;s<r;s++){ 
    if(cflst[s]<0){
      ns=numsplits(s,numv,explst,cflst,r,&needcount,&needexcl,needs,neededby,0);
      //fprintf(stderr, "cflst[%d] = %d, ns=%.4f\n", s, cflst[s],ns);
      if(ns<0.0001){// not able to split any more
	cout << "failure to split:\n   " << cflst[s] << ", "; printvec1(explst[s],numv);
	(*level)--;
	ifree_l(needs, r);ifree_l(neededby, r);
	free((char *)needcount);free((char *)needexcl);
	return 0; // exit, but not fatal
      }
    }
  }
  totcliques=allcliques1(explst,cflst,r,numv,neededby,needs,needexcl,&lt,&rt,&numlt,&numrt,&absnegmax,maxsz,0,0,1,0);
  printf("totcliques=%d\n", totcliques);
  ifree_l(needs, r);ifree_l(neededby, r);
  free((char *)needcount);free((char *)needexcl);
      
  for(j=0;j<totcliques;j++){// remove positive quantities
    veccp(cflst,r,cflst1);veccp(changelst,r,changelst1);
    tmp2=0;
    if(!rmsq(explst,cflst1,changelst1,r,lt[j],numlt[j],rt[j],numrt[j],numv,absnegmax[j],&tmp2,1)){//remove a square
      cout << "ERROR: an apparent square wasn't found.\n";
      exit(0);
    }
    if(sqrem(explst, cflst1, r, numv, level, changelst1,totiter,sqout,sqcfs,sqoutlen, numsqout)){//success

      cout << "removing at least: " << tmp2 << endl;
      // as we exit successfully add to the outsquares

      //printvec1(changelst,r);cout << endl;
      //printvec1(changelst1,r);cout << endl;
      newblock=cfchangetoblock(explst, changelst,r,changelst1,&blocksz,&blockcfs,numv);
      if(blocksz>0){
	//printblock1(newblock,blockcfs,blocksz,numv);
	addnewblock(sqout, *numsqout, newblock, blocksz, numv);
	addnewto1Darray(sqoutlen,*numsqout,blocksz);
	addnewtoarray(sqcfs,*numsqout,blockcfs,blocksz);
	(*numsqout)++;
	ifree(newblock,blocksz);
	free((char*)blockcfs);
      }
 
      //veccp(cflst1,r,cflst);veccp(changelst1,r,changelst);
      if(lt){lfree_i(lt,totcliques);}if(rt){lfree_i(rt,totcliques);}
      if(numlt){free((char*)numlt);}if(numrt){free((char*)numrt);}
      if(absnegmax){free((char*)absnegmax);}
      //fprintf(stderr, "exiting sqrem\n");

      return 1;
    }
    //    fprintf(stderr, "looping.\n");
  }
  if(lt){lfree_i(lt,totcliques);}if(rt){lfree_i(rt,totcliques);}
  if(numlt){free((char*)numlt);}if(numrt){free((char*)numrt);}if(absnegmax){free((char*)absnegmax);}
  //fprintf(stderr, "exiting without success.\n");
  (*level)--;//fprintf(stderr, "level = %d\n", *level);
  return 0; 
}

bool getpsq(int **explst, int *cflst, long r, int numv, int ***sqs, int *sz, int **cfs, int numsqs, int *changelst, int ****sqout, int ***sqcfs, int **sqoutlen, int *numsqout){
  long s;
  int negtot=0;
  int level=0;
  int totiter=0;

  for(s=0;s<r;s++){
    if(cflst[s]<0)
      negtot+=abs(cflst[s]);
  }

  if(negtot==0)// no negative terms
    return 1;

  /* for(i=0;i<numsqs;i++){ // remove all existing squares (alters cflst) */
  /*   if(sz[i]>0 && removelst(explst, cflst, r, sqs[i], cfs[i], sz[i], numv, 0)){ */
  /*     numrem++; */
  /*   } */
  /* } */
  /* if(numrem){ */
  /*   cout << "Removed " << numrem << " terms\n"; */
  /*   for(s=0;s<r;s++){ */
  /*     if(cflst[s]<0) */
  /* 	negtot+=abs(cflst[s]); */
  /*   } */
  /*   if(negtot==0)// no negative terms */
  /*     return 1; */
  /* } */

  // try to remove all squares
  // if successful output what is removed as a new list
  if(!sqrem(explst, cflst, r, numv, &level,changelst,&totiter, sqout, sqcfs, sqoutlen, numsqout)){
    fprintf(stderr, "failed\n");
    return 0;
  }

  fprintf(stderr, "successfully removed all squares\n");
  return 1;
}

void addonblocks(int ****sqs, int ***sqcfs, int **sqlen, int *numsqs, int ***sqnew, int **sqcfsnew, int *sqlennew, int numsqnew, int numv){
  int i;
  for(i=0;i<numsqnew;i++){
    //printblock1(newblock,blockcfs,blocksz,numv);
    addnewblock(sqs,*numsqs, sqnew[i], sqlennew[i], numv);
    addnewtoarray(sqcfs,*numsqs,sqcfsnew[i],sqlennew[i]);
    addnewto1Darray(sqlen,*numsqs,sqlennew[i]);
    (*numsqs)++;
  }
  return;
}

matrix Jsubmat(matrix J, int n, int m, int *xc, int *yc, int xlen, int ylen){
  matrix J1(xlen,ylen);
  int i,j;
  if(xlen>n || ylen > m){
    fprintf(stderr, "Can't make a submatrix of dimensions %d X %d from a matrix of dimension %d X %d\n", xlen, ylen, n, m);
    exit(0);
  }
  for(i=0;i<xlen;i++){// submatrix
    for(j=0;j<ylen;j++)
      J1(i,j)=J(xc[i],xc[j]);
  }
  return J1;
}

void J2pI_trial(matrix J, int n, int numv){
  ex v[numv+2];
  char str1[10];
  ex tmp,tmp2,tmp3;
  int i,j,k,dim,combflag;
  int *xc;
  ex out=0;
  // global list of useful squares
  int ***sqs=NULL,**sqcfs=NULL, *sqlen=NULL, numsqs=0;
  int **lst, *cflst, **explst;
  int *changelst=NULL;
  int negcount=0;
  long r, s;
  int ***sqout=NULL, **sqoutcfs=NULL, *sqoutlen=NULL, numsqout=0;
  matrix J1;

  for(k=0;k<=numv+1;k++){
    sprintf(str1, "v%d", k);
    v[k]=get_possymbol(str1);
  }

  printexmat(J,n,n);

  for(k=1;k<n;k++){ // one dimension at a time
    xc=(int *)malloc((size_t) (k*sizeof(int)));
    firstcomb(xc,n,k);
    combflag=1;

    // as we move up a level increase the powers of omega by two

    for(i=0;i<numsqs;i++){
      for(j=0;j<sqlen[i];j++){
	sqs[i][j][numv]+=2;
      } 
    }
    printf("numsqs = %d\n", numsqs);
    cout << "Current list of stored squares:\n";
    for(i=0;i<numsqs;i++){
      printf("sqlen[%d] = %d\n", i, sqlen[i]);
      for(j=0;j<sqlen[i];j++){
	cout << sqcfs[i][j] << ", ";
	printvec1(sqs[i][j],numv+1);
      }
      tmp3=expltopolysimp(sqs[i], sqcfs[i], sqlen[i], numv+1);
      cout << tmp3 << endl << endl;
    }

    while(combflag==1){// each principal submatrix

      negcount=0;numsqout=0;
      printvec1(xc,k);
      J1=Jsubmat(J,n,n,xc,xc,k,k);
      printexmat(J1,k,k);
      tmp=0;
      for(i=0;i<=k;i++)
	tmp+=expand(pow(v[numv+1], 2*(k-i))*getminorsum(J1,k,i));

      cout << tmp << endl;

      dim=polytointlist(tmp, &lst, &cflst, &r);//dim=order of monomials
      explst=moncflisttoexplist(lst,dim,cflst,r,numv+1,1);// redefines cflst
      changelst=(int *)malloc((size_t) (r*sizeof(int)));
      for(s=0;s<r;s++){
	if(cflst[s]<0)
	  negcount+=abs(cflst[s]);
	changelst[s]=0; // initialise the list of changes
      }
      if(negcount!=0){ // try to remove 

	if(!getpsq(explst, cflst, r, numv+1, sqs, sqlen, sqcfs, numsqs,changelst,&sqout,&sqoutcfs,&sqoutlen,&numsqout)){
	  fprintf(stderr, "Failed at dimension %d. Couldn't remove negative terms.\n", k);
	  exit(0);
	}
	//printvec(changelst,r);
	//sanity check
	tmp2=tmp;
	for(i=0;i<numsqout;i++){
	  //fprintf(stderr, "sqoutlen[%d]=%d\n",i,sqoutlen[i]);
	  cout << "Removing...\n";
	  tmp3=expltopolysimp(sqout[i], sqoutcfs[i], sqoutlen[i], numv+1);
	  cout << tmp3 << endl << endl;
	  tmp2-=tmp3;
	}

	cout<< "final polynomial:\n";
	cout << tmp2 << endl;

	addonblocks(&sqs, &sqcfs, &sqlen, &numsqs, sqout, sqoutcfs, sqoutlen, numsqout,numv+1);
	printf("numsqs = %d\n", numsqs);
	cout << "Full list of stored squares:\n";
	for(i=0;i<numsqs;i++){
	  printf("sqlen[%d] = %d\n", i, sqlen[i]);
	  for(j=0;j<sqlen[i];j++){
	    cout << sqcfs[i][j] << ", ";
	    printvec1(sqs[i][j],numv+1);
	  }
	  tmp3=expltopolysimp(sqs[i], sqcfs[i], sqlen[i], numv+1);
	  cout << tmp3 << endl << endl;
	}
	freeblockcflist(sqout, sqoutcfs, sqoutlen, numsqout);
      }

      combflag=nextcomb(xc,n,k);
      ifree_l(lst, r);
      ifree_l(explst, r);
      free((char *)cflst);
      free((char *)changelst);
    }

    free((char *)xc);
  }


  return;

}






matrix multAB(int **S, matrix V, int n, int m){
  matrix J(n,n);
  int i,j,k;
  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      J(i,j)=0;
      for(k=0;k<m;k++)
	J(i,j)+=S[i][k]*V(j,k);
    }
  }
  return J;
}

bool flgset(int newf, int a, int b){
  if(newf==2){
    if(!a && !b)
      return 1;
  }
  else if(newf==1){
    if((!a || !b) && a>=0 && b>=0)
      return 1;
  }
  else{
    if(a>=0 && b>=0 && (a!=b || (a==0 && b==0)))
      return 1;
  }
  return 0;
}

// new=2 is a new pair; new=1 has at least one new member; new=0 allows both members to be in use

int validpair1(int **explst, long t1, long t2, long *u1, long *u2, int **explst1, int **explst2, int *cflst1, int *cflst2, long r1, long r2, int *used1, int *used2, int numv, int newf){
  long s1, s2;

  //try list 1
  for(s1=0;s1<r1;s1++){
    if(ueq2v(explst[t1],explst1[s1], numv)){
      for(s2=s1+1;s2<r1;s2++){
	//cout << used1[s1]<< ", " << used1[s2] << endl;
	//cout << "flgset = " << flgset(newf, used1[s1], used1[s2]) << endl;
	if(ueq2v(explst[t2],explst1[s2], numv) &&  flgset(newf, used1[s1], used1[s2])){
	  if(cflst1[s1]*cflst1[s2]<0){
	    if(cflst1[s1]>0){*u1=s1;*u2=s2;}else{*u1=s2;*u2=s1;}
	    return 1;
	  }
	}
      }
    }
    else if(ueq2v(explst[t2],explst1[s1], numv)){
     for(s2=s1+1;s2<r1;s2++){
       //cout << used1[s1]<< ", " << used1[s2] << endl;
       //cout << "flgset = " << flgset(newf, used1[s1], used1[s2]) << endl;
	if(ueq2v(explst[t1],explst1[s2], numv) && flgset(newf, used1[s1], used1[s2])){
	  if(cflst1[s1]*cflst1[s2]<0){
	    if(cflst1[s1]>0){*u1=s1;*u2=s2;}else{*u1=s2;*u2=s1;}
	    return 1;
	  }
	}
      }
    }
  }
  // try list 2
  for(s1=0;s1<r2;s1++){
    if(ueq2v(explst[t1],explst2[s1], numv)){
      for(s2=s1+1;s2<r2;s2++){
	//cout << used2[s1]<< ", " << used2[s2] << endl;
	//cout << "flgset = " << flgset(newf, used2[s1], used2[s2]) << endl;
	if(ueq2v(explst[t2],explst2[s2], numv) && flgset(newf, used2[s1], used2[s2])){
	  if(cflst2[s1]*cflst2[s2]<0){
	    if(cflst2[s1]>0){*u1=s1;*u2=s2;}else{*u1=s2;*u2=s1;}
	    return 2;
	  }
	}
      }
    }
    else if(ueq2v(explst[t2],explst2[s1], numv)){
      for(s2=s1+1;s2<r2;s2++){
	//cout << used2[s1]<< ", " << used2[s2] << endl;
	//cout << "flgset = " << flgset(newf, used2[s1], used2[s2]) << endl;
	if(ueq2v(explst[t1],explst2[s2], numv) && flgset(newf, used2[s1], used2[s2])){
	  if(cflst2[s1]*cflst2[s2]<0){
	    if(cflst2[s1]>0){*u1=s1;*u2=s2;}else{*u1=s2;*u2=s1;}
	    return 2;
	  }
	}
      }
    }
  }
  return 0;
}


bool isprod(int *a, int *b, int *prod, int numv){
  int i;
  for(i=0;i<numv;i++){
    if((a[i]+b[i])!=prod[i])
      return 0;
  }
  return 1;
}







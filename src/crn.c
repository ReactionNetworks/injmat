/* Copyright (C) 2010-2013, Murad Banaji
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

int sgn(int i){
  if(!i)
    return 0;
  if(i>0)
    return 1;
  return -1;
}

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

/* Check if a minor is sign singular only */
/* Everything else returns failure */

int minorisSS(int **imat, int n, int m, int *vec1, int *vec2, int k, unsigned long fk, int **pms){
  unsigned long i;
  /* worth checking for a row or column of zeros first */
  if(minorhas0rc(imat, n, m, vec1, vec2, k))
    return 1;
  for(i=0;i<fk;i++){
    if(unsterm1(imat,n,m,vec1,pms[i],k)) /* nonzero term */
      return 0;
  }
  return 1;
}

/* like minorisSS except doesn't first compute */
/* all the permutations, so should avoid memory failure */

int minorisSS2(int **imat, int n, int m, int *vec1, int *vec2, int k, unsigned long fk){
  int *vec;
  int *veclr;
  int i, par=1, flag=1;
  unsigned long t=0;

  /* worth checking for a row or column of zeros first */
  if(minorhas0rc(imat, n, m, vec1, vec2, k))
    return 1;

  vec=(int *)malloc((size_t) ((k)*sizeof(int)));
  veclr=(int *)malloc((size_t) ((k)*sizeof(int)));

  for(i=0;i<k;i++){vec[i]=vec2[i];veclr[i]=-1;}

  while(flag){
    if(unsterm1(imat, n, m, vec1, vec,k)) /* nonzero term */
      return 0;
    flag=nextperm(&vec, &veclr, &par, k);
    if(t%1000000==0){
      fprintf(stderr, "%d percent complete\n", (int)(((double)t)/((double)fk)*100));
    }
    t++;
  }
  free((char *) vec);free((char *) veclr);
  return 1;
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

int minorisSNSSS_rec(int **imat, int n, int m, int *vec1, int *vec2, int k){

  /* worth checking for a column or row of zeros first */

  if(minorhas0rc(imat, n, m, vec1, vec2, k))
    return 1;
  if(qualdetsubmat(imat, n, m, vec1, vec2, k)!=2)
    return 1;
  return 0;

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

int minorisSNSSS2(int **imat, int n, int m, int *vec1, int *vec2, int k, unsigned long fk){
  int tmp, tm1=0; 
  int *vec;
  int *veclr;
  int i, par=1, flag=1;
  unsigned long t=0;

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
    /* if(t%1000000==0){ */
    /*   fprintf(stderr, "%d percent complete\n", (int)(((double)t)/((double)fk)*100)); */
    /* } */
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

// is column k of matrix mat (with n rows) equal or plus/minus an earlier column?

bool eqorminus(int **mat, int n, int k){
  int i=0,j;
  int flgminus=1;
  int flg=1;
  if(k==0)
    return 0;

  for(j=0;j<k;j++){ //each col
    i=0;flgminus=1;flg=1;
    while(i<n && mat[i][j]==0 && mat[i][k]==0)//initial zeros
      i++;
    if(mat[i][j]==mat[i][k])
      i++;
    else if (mat[i][j]==-mat[i][k]){
      flgminus=-1;i++;
    }
    else
      flg=0;
    while(flg && i<n){
      if(mat[i][j]!=flgminus*mat[i][k])
	flg=0;
      i++;
    }
    if(flg)
      return 1;
  }
  return 0;
}

// list repeated columns (possibly after a sign change) 
// of an n X m matrix mat. Store the information in a binary vector
// of what to keep and what to discard

int reps(int **mat, int n, int m, bool **keeps){
  int k, tot=1;
  (*keeps)[0]=1;
  for(k=1;k<m;k++){// keep the first column
    if(eqorminus(mat, n, k))
      (*keeps)[k]=0;
    else{
      (*keeps)[k]=1;tot++;
    }
  }
  return tot;
}

// remove redundant columns from a stoichiometric matrix
// i.e. columns which (upto sign change) appear previously 

int **redmat(int **mat, int n, int m, int *m1){
  bool *keeps=(bool *)malloc((size_t) ((m)*sizeof(bool)));
  int **tmp;
  int i,j,tot=0;

  (*m1)=reps(mat,n,m,&keeps);
  tmp=imatrix(0, n-1, 0, (*m1)-1);
  for(i=0;i<m;i++){//column
    if(keeps[i]){
      for(j=0;j<n;j++)//row
	tmp[j][tot]=mat[j][i];
      tot++;
    }
  }
  return tmp;
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
  /* unsigned long fk=factorial(k); */

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
    /* if(t%1000000==0){ */
    /*   fprintf(stderr, "%d percent complete\n", (int)(((double)t)/((double)fk)*100)); */
    /* } */
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
  int **pms,j;
  int dt=0,qd=0, tmp;
  unsigned long i, fk;
  int tm1=0;
  int par=1,flag=1;
  unsigned long t=0;
  int *vec, *veclr;
  int rtrig=8;
 
  /* worth checking if either is SS first */

  if(minorhas0rc(imat2, n, m, vec1, vec2, k) || minorhas0rc(imat1, n, m, vec1, vec2, k))
    return 1;

  if((dt=detsubmat(imat1, n, m, vec1, vec2, k))==0)
    return 1;
  if((qd=qualdetsubmat(imat2, n, m, vec1, vec2, k))==0)
    return 1;
  if(qd==2)
    return 0;
  if(dt*qd<0) // opp signs: keep track of this info
    return -1;
  if(dt*qd>0) // a positive product: keep track of this info
    return 2;

  fk=factorial(k);
  if(k<=rtrig){
    /* fprintf(stderr, "entering here.\n"); */
    pms=allperms1(vec2, k);
    for(i=0;i<fk;i++){
      tmp=unsterm1(imat2, n, m, vec1, pms[i],k);
      if(tmp){
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
    free_imatrix(pms,0, fk-1, 0, k);
  }
  else{ // large submatrix
    vec=(int *)malloc((size_t) ((k)*sizeof(int)));
    veclr=(int *)malloc((size_t) ((k)*sizeof(int)));

    for(j=0;j<k;i++){vec[j]=vec2[j];veclr[j]=-1;}

    while(flag){
      tmp=unsterm1(imat2, n, m, vec1, vec, k);
      if(tmp){
	if(abs(tmp)>1){ // unsigned term
	  free((char *) vec);free((char *) veclr);
	  return 0;
	}
	// fails to be SNS
	if(tm1==0)
	  tm1=par*tmp;
	else if(par*tmp*tm1 < 0){
	  free((char *) vec);free((char *) veclr);
	  return 0;
	}
      }
      flag=nextperm(&vec, &veclr, &par, k);
      /* if(t%1000000==0){ */
      /* 	fprintf(stderr, "%d percent complete\n", (int)(((double)t)/((double)fk)*100)); */
      /* } */
      t++;
    }

    free((char *) vec);free((char *) veclr);
  }
  if(dt*tm1<0) // opp signs: keep track of this info
    return -1;
  else if(dt*tm1>0) // a positive product: keep track of this info
    return 2;

  return 1;
}

/* imat1 is an n X m matrix, while imat2 is an n X m sign-pattern */
/* check if they are compatible. */
/* flag = 3 means the matrices are compatible and r-strongly compatible */
/* flag = 2 means the matrices are compatible but not r-strongly compatible */
/* flag = 1 means the matrices are r-strongly compatible but not compatible */
/* flag = 0 means the matrices are none of the above */

int mat_signpat_compat0(int **imat1, int imatrank, int **imat2, int n, int m, int q){
  int k;
  long r1, r2, cnk, cmk;
  int **xcombs;
  int **ycombs;
  int flag=3, flg;
  //  int imatrank=matrank(imat1,n,m);
  bool posprod=0;

  for(k=imatrank;k>=1;k--){
    if(k==imatrank-1 && !posprod) // not r-strongly compatible
      flag=2;
    xcombs=allcombsgen(n,k);
    ycombs=allcombsgen(m,k);
    cnk=comb(n, k);cmk=comb(m, k);
    if(n+m>15){
      fprintf(stderr, "\nchecking %.0f minors of size %d...\n", ((double) cnk)*((double) (cmk)), k);
    }
    for(r1=0;r1<cnk;r1++){
      if(r1%100==0 && (n+m>15))
	fprintf(stderr, ".");
      for(r2=0;r2<cmk;r2++){
	flg=submat_signpat_compat(imat1, imat2, n, m, xcombs[r1], ycombs[r2], k);
	if(flg<=0){ // fail to be compatible: return 1 or 0
	  fprintf(stderr, "\nnumerical/pattern submatrices which fail to be compatible:\n\n");
	  printsubmat(imat1, xcombs[r1], ycombs[r2], k, k);
	  fprintf(stderr, "*** and ***\n");
	  printsubmat(imat2, xcombs[r1], ycombs[r2], k, k);
	  if(k==imatrank){flag=0;}else if(flag>=2){flag=flag-2;}
	  if(q){
	    free_imatrix(xcombs, 0, cnk-1, 0, k-1);
	    free_imatrix(ycombs, 0, cmk-1, 0, k-1);
	    return flag;
	  }
	}
	else if(flg==2 && k==imatrank) // a strictly positive product
	  posprod=1;

      }
    }
    free_imatrix(xcombs, 0, cnk-1, 0, k-1);
    free_imatrix(ycombs, 0, cmk-1, 0, k-1);
  }

  return flag;
}

// are a matrix and sign pattern negative r-compatible

int mat_signpat_negr_compat(int **imat1, int imatrank, int **imat2, int n, int m, int q){
  int k=imatrank;
  long r1, r2, cnk, cmk;
  int **xcombs;
  int **ycombs;
  int flag=1, flg;
  int posprod=0;
  fprintf(stderr, "\nChecking negative compatibility...\n");
  xcombs=allcombsgen(n,k);
  ycombs=allcombsgen(m,k);
  cnk=comb(n, k);cmk=comb(m, k);
  fprintf(stderr, "\nchecking %.0f minors of size %d...\n", ((double) cnk)*((double) (cmk)), k);
  for(r1=0;r1<cnk;r1++){
    if(r1%100==0)
      fprintf(stderr, ".");
    for(r2=0;r2<cmk;r2++){
      flg=submat_signpat_compat(imat1, imat2, n, m, xcombs[r1], ycombs[r2], k);
      if(flg!=-1 && flg!=1){// fail to be negatively compatible
	fprintf(stderr, "\nnumerical/pattern submatrices which fail to be negatively compatible:\n\n");
	printsubmat(imat1, xcombs[r1], ycombs[r2], k, k);
	fprintf(stderr, "*** and ***\n");
	printsubmat(imat2, xcombs[r1], ycombs[r2], k, k);
	flag=0;
	if(q){
	  free_imatrix(xcombs, 0, cnk-1, 0, k-1);
	  free_imatrix(ycombs, 0, cmk-1, 0, k-1);
	  return flag;
	}
      }
      else if(flg==-1)
	posprod=1;
    }
  }
  free_imatrix(xcombs, 0, cnk-1, 0, k-1);
  free_imatrix(ycombs, 0, cmk-1, 0, k-1);

  if(flag && posprod)
    return 1;
  return 0;
}

// wrapper for mat_signpat_compat0

int mat_signpat_compat(int **imat1, int imatrank, int **imat2, int n, int m, int q){
  int flg;
  if(n+m>15)
    fprintf(stderr, "Depending on the speed of your computer, this could take a while.\n");

  flg=mat_signpat_compat0(imat1, imatrank, imat2, n, m, q);
  if(flg)
    return flg;
  else if(mat_signpat_negr_compat(imat1, imatrank, imat2, n, m, q))
    return -1;
  return 0;
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
      if(submat_signpat_compat(imat1, imat2, n, m, xcombs[r1], ycombs[r2], k)<=0){
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




/* Check if n X m matrix imat1 and n X m sign pattern imat2 */
/* are compatible. Remove empty rows/columns first. This is */
/* done twice because in theory each removal can create a new */
/* row/column of zeros in the other matrix */

int arecompat(int **imat1, int **imat2, int n, int m, int q){
  int **imat1a, **imat2a, **imat1b, **imat2b;
  int n1, m1, n2, m2;
  int flag=0;
  int imatrank=matrank(imat1, n, m);
  // simplify first
  fprintf(stderr, "Checking compatibility...\n");
  fprintf(stderr, "Removing rows and columns of zeros...\n");
  simppair(imat1, imat2, n, m, &imat1a, &imat2a, &n1, &m1);
  simppair(imat1a, imat2a, n1, m1, &imat1b, &imat2b, &n2, &m2);
  printmat(imat1b, n2, m2);
  printmat(imat2b, n2, m2);
  fprintf(stderr, "Checking if matrix and sign-pattern are compatible, %d-strongly compatible or %d-strongly negatively compatible...\n", imatrank, imatrank);
  flag=mat_signpat_compat(imat1b, imatrank, imat2b, n2, m2, q);

  if(flag==3)
    fprintf(stderr, "Finished checking compatibility: matrix and sign-pattern are compatible and %d-strongly compatible.\n---------------------------------------\n",imatrank);
  else if(flag==2)
    fprintf(stderr, "Finished checking compatibility: matrix and sign-pattern are compatible but not %d-strongly compatible.\n---------------------------------------\n", imatrank);
  else if(flag==1)
    fprintf(stderr, "Finished checking compatibility: matrix and sign-pattern are %d-strongly compatible but not compatible.\n---------------------------------------\n", imatrank);
  else if(flag==-1)
    fprintf(stderr, "Finished checking negative compatibility: matrix and sign-pattern are %d-strongly negatively compatible.\n---------------------------------------\n", imatrank);
  else
    fprintf(stderr, "Finished checking compatibility: matrix and sign-pattern are not compatible or %d-strongly compatible or %d-strongly negatively compatible.\n---------------------------------------\n", imatrank, imatrank);

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
  int flag=0;
  int **imat, n1, m1;
  int imatrank=matrank(mat,n,m);
  int rkbad=0; // rank of bad minors found
  fprintf(stderr, "Checking if the matrix is SSD or %d-SSD...\n",imatrank);
  if(n+m > 15)
    fprintf(stderr, "...please be patient. This could take time.\n");
  imat=simpmat(mat, n, m, &n1, &m1);
  flag=isSSD1(imat, n1, m1, 1, &rkbad, q);
  free_imatrix(imat,0, n1-1, 0, m1-1);
  if(flag){
    fprintf(stderr, "Finished checking SSD: the matrix is SSD.\n---------------------------------------\n");
    return 2;
  }
  else if(rkbad==imatrank){
    fprintf(stderr, "Finished checking SSD: the matrix is not SSD or %d-SSD.\n---------------------------------------\n",imatrank);
    return 0;
  }
  flag=isSSD1(mat, n, m, 0, &rkbad, q);

  if(flag){
    fprintf(stderr, "Finished checking SSD: the matrix is %d-SSD but not SSD.\n---------------------------------------\n",imatrank);
    return 1;
  }
  else{
    fprintf(stderr, "Finished checking SSD: the matrix is not SSD or %d-SSD.\n---------------------------------------\n",imatrank);
    return 0;
  }
}



/* Check if n X m matrix mat is CSD (All square */
/* submatrices are sign nonsingular or sign singular) */

/* flag = 2 means CSD */
/* flag = 1 means r-CSD but not CSD */
/* flag = 0 means neither */


int isCSD1(int **imat, int n, int m, bool allm, int *rkbad, int q){
  int k,r;
  long r1, r2, cnk, cmk;
  int **xcombs;
  int **ycombs;
  int flag=1;
  unsigned long fk;
  int **pms;
  int ret;
  int rtrig=4;
  int rmin=2;
  int imatrank=matrank(imat,n,m);

  if(allm){
    r=min(n,m);rmin=2;
  }
  else{
    r=imatrank;rmin=r;
  }

  for(k=r;k>=rmin;k--){
    //fprintf(stderr, "k=%d\n", k);
    xcombs=allcombsgen(n,k);ycombs=allcombsgen(m,k);
    cnk=comb(n, k);cmk=comb(m, k);
    fk=factorial(k);
    if(k<=rtrig){// smallish matrix
      for(r2=0;r2<cmk;r2++){
	pms=allperms1(ycombs[r2], k);
	for(r1=0;r1<cnk;r1++){
	  if(k>imatrank)// only check for sign singularity
	    ret=minorisSS(imat, n, m, xcombs[r1], ycombs[r2], k, fk, pms);
	  else
	    ret=minorisSNSSS(imat, n, m, xcombs[r1], ycombs[r2], k, fk, pms);

	  if(!ret){
	    fprintf(stderr, "submatrix which fails to be sign nonsingular or sign singular:\n");
	    printsubmat(imat, xcombs[r1], ycombs[r2], k, k);
	    flag=0;if(!(*rkbad)){(*rkbad)=k;}// rank of bad submatrix
	    if(q){
	      free_imatrix(pms,0, fk-1, 0, k);
	      free_imatrix(xcombs, 0, cnk-1, 0, k-1);
	      free_imatrix(ycombs, 0, cmk-1, 0, k-1);
	      return flag;
	    }
	  }
	}
      }
    }
    else{// large matrix
      for(r2=0;r2<cmk;r2++){
	for(r1=0;r1<cnk;r1++){
	  if(k>imatrank)// only check for sign singularity
	    ret=minorisSS2(imat, n, m, xcombs[r1], ycombs[r2], k, fk);
	  else
	    ret=minorisSNSSS_rec(imat, n, m, xcombs[r1], ycombs[r2], k);//recursive version
	  //	    ret=minorisSNSSS2(imat, n, m, xcombs[r1], ycombs[r2], k,fk);

	  if(!ret){
	    fprintf(stderr, "submatrix which fails to be sign nonsingular or sign singular:\n");
	    printsubmat(imat, xcombs[r1], ycombs[r2], k, k);
	    flag=0;if(!(*rkbad)){(*rkbad)=k;}// rank of bad submatrix
	    if(q){
	      free_imatrix(xcombs, 0, cnk-1, 0, k-1);
	      free_imatrix(ycombs, 0, cmk-1, 0, k-1);
	      return flag;
	    }
	  }

	}
      }
    }
    free_imatrix(xcombs, 0, cnk-1, 0, k-1);
    free_imatrix(ycombs, 0, cmk-1, 0, k-1);
  }

  return flag;
}

/* Check if n X m matrix imat is SSD */

int isSSD1(int **imat, int n, int m, bool allm, int *rkbad, int q){
  int k, ret=0;
  int r=100;
  long r1, r2, cnk, cmk;
  int **xcombs;
  int **ycombs;
  int flag=1;
  unsigned long fk,t=0;
  int **pms;
  int rmin=2,szswtch=0;
  int imatrank=matrank(imat,n,m);


  if(allm){
    r=imatrank;rmin=2;
  }
  else{
    r=imatrank;rmin=r;
  }

  // all minors of size greater than imatrank are definitely singular, so nothing to check
  
  for(k=r;k>=rmin;k--){
    xcombs=allcombsgen(n,k);
    ycombs=allcombsgen(m,k);
    cnk=comb(n, k);cmk=comb(m, k);
    fprintf(stderr, "\nchecking %.0f minors of size %d...\n", ((double) cnk)*((double) (cmk)), k);
    fk=factorial(k);

    t=0;
    for(r2=0;r2<cmk;r2++){
      /* if(t%100==0) */
      /* 	fprintf(stderr, "t=%.0f\n", (double) t); */
      /* t++; */

      if(k<=szswtch)
	pms=allperms1(ycombs[r2], k);
      for(r1=0;r1<cnk;r1++){
	if(k>szswtch) // large submatrix
	  ret=minorisSNSsing2(imat, n, m, xcombs[r1], ycombs[r2], k);
	else
	  ret=minorisSNSsing(imat, n, m, xcombs[r1], ycombs[r2], k, fk, pms);

	if(!ret){
	  fprintf(stderr, "submatrix which fails to be sign-nonsingular or singular (t=%.0f):\n", (double) t);
	  //fprintf(stderr, "det = %d\n", detsubmat(imat, n, m, xcombs[r1], ycombs[r2], k));
	  //fprintf(stderr, "qualdet = %d\n", qualdetsubmat(imat, n, m, xcombs[r1], ycombs[r2], k-1));
	  printsubmat(imat, xcombs[r1], ycombs[r2], k, k);
	  flag=0;if(!(*rkbad)){(*rkbad)=k;}// rank of bad submatrix
	  if(q){ // exit here for speed
	    if(k<=szswtch)
	      free_imatrix(pms,0, fk-1, 0, k);
	    free_imatrix(xcombs, 0, cnk-1, 0, k-1);
	    free_imatrix(ycombs, 0, cmk-1, 0, k-1);
	    return flag;
	  }
	}
      }
      if(k<=szswtch)
	free_imatrix(pms,0, fk-1, 0, k);
    }
    free_imatrix(xcombs, 0, cnk-1, 0, k-1);
    free_imatrix(ycombs, 0, cmk-1, 0, k-1);
  }

  return flag; // SSD (and hence automatically r-SSD)
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
  int flag;
  int **imat, n1, m1;
  int imatrank=matrank(mat,n,m);
  int rkbad=0; // rank of bad minors found
  fprintf(stderr, "Checking if the matrix is CSD or %d-CSD...\n", imatrank);
  printmat(mat, n, m);
  imat=simpmat(mat, n, m, &n1, &m1);
  flag=isCSD1(imat, n1, m1, 1, &rkbad, q);
  free_imatrix(imat,0, n1-1, 0, m1-1);

  if(flag){
    fprintf(stderr, "Finished checking CSD: the matrix is CSD.\n---------------------------------------\n");
    return 2;
  }
  else if(rkbad==imatrank){
    fprintf(stderr, "Finished checking CSD: the matrix is not CSD or %d-CSD.\n---------------------------------------\n", imatrank);
    return 0;
  }

  flag=isCSD1(mat, n, m, 0, &rkbad, q);

  if(flag){
    fprintf(stderr, "Finished checking CSD: the matrix is %d-CSD but not CSD.\n---------------------------------------\n", imatrank);
    return 1;
  }
  else{
    fprintf(stderr, "Finished checking CSD: the matrix is not CSD or %d-CSD.\n---------------------------------------\n", imatrank);
    return 0;
  }
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
  int n1,m1,flag1=0,flag2=0;
  int rkbad=0;
  int imatrank=matrank(mat,n,m);
  imat=simpmat(mat, n, m, &n1, &m1); // first simplify and then double
  tmp=doublemat(imat, n1, m1);
  //fprintf(stderr, "Checking if this matrix is WSD...\n");
  //printmat(tmp, n1, 2*m1);
  flag1=isWSD1(tmp, n1, 2*m1, 1, &rkbad, q); //all
  free_imatrix(tmp, 0, n1-1, 0, 2*m1-1);
  free_imatrix(imat,0, n1-1, 0, m1-1);
  if(flag1==0 && rkbad==imatrank){// no need to check rWSD
    fprintf(stderr, "Finished checking WSD: the matrix is neither r-strongly WSD nor WSD.\n---------------------------------------\n");
    return 0;
  }
  // for r-strong WSD, don't simplify
  tmp=doublemat(mat, n, m);
  fprintf(stderr, "Checking if this matrix is WSD or r-strongly WSD...\n");
  printmat(tmp, n, 2*m);
  flag2=isWSD1(tmp, n, 2*m, 0, &rkbad, q); // only maximal
  free_imatrix(tmp, 0, n-1, 0, 2*m-1);

  if(flag1 && flag2){
    fprintf(stderr, "Finished checking WSD: the matrix is WSD and r-strongly WSD.\n---------------------------------------\n");
    return 3;
  }
  else if (flag1){
    fprintf(stderr, "Finished checking WSD: the matrix is WSD but not r-strongly WSD.\n---------------------------------------\n");
    return 2;
  }
  else if(flag2){
    fprintf(stderr, "Finished checking WSD: the matrix is r-strongly WSD but not WSD.\n---------------------------------------\n");
    return 1;
  }
  else{
    fprintf(stderr, "Finished checking WSD: the matrix is neither r-strongly WSD nor WSD.\n---------------------------------------\n");
    return 0;
  }


}



/* A wrapper for isWSD1, first removing rows/columns of zeros */

int isWSD(int **mat, int n, int m, int q){
  int **imat;
  int n1,m1,flag=0;
  int rkbad=0;
  imat=simpmat(mat, n, m, &n1, &m1);
  flag=isWSD1(imat, n1, 2*m1, 1, &rkbad, q);
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

/* flag = 3 means the matrices are compatible and r-strongly compatible */
/* flag = 2 means the matrices are compatible but not r-strongly compatible */
/* flag = 1 means the matrices are r-strongly compatible but not compatible */
/* flag = 0 means the matrices are none of the above */

int isWSD1(int **imat, int n, int m, bool allm, int *rkbad, int q){
  int k,rmin,rmax,tt;
  long r1,r2,cnk,cmk;
  int **xcombs;
  int **ycombs;
  int flag=1;
  int **imat2, **imat3, n1, m1;
  int mr1, mr2;
  bool posprod=0;

  WSpair(imat, n, m, &imat2, &imat3, &n1, &m1); // imat3 is sparser

  if(n1!=n || m1!=m){
    fprintf(stderr, "The simplified matrix pair...\n");
    printmat(imat2, n1, m1);
    printmat(imat3, n1, m1);
  }

  mr1=matrank(imat2,n1,m1);mr2=matrank(imat3,n1,m1);
  rmax=min(mr1,mr2);

  if(!allm && rmax<mr1){ // can't be r-strongly compatible
    free_imatrix(imat2, 0, n1-1, 0, m1-1);
    free_imatrix(imat3, 0, n1-1, 0, m1-1);
    return 0;
  }

  if(allm)
    rmin=2;
  else
    rmin=rmax;

  for(k=rmax;k>=rmin;k--){
    xcombs=allcombsgen(n1,k);ycombs=allcombsgen(m1,k);
    cnk=comb(n1,k);cmk=comb(m1,k);
    for(r1=0;r1<cnk;r1++){
      for(r2=0;r2<cmk;r2++){
	if((tt=fixedminorcompat(imat2, imat3, n1, m1, xcombs[r1], ycombs[r2],k))<0){
	  fprintf(stderr, "submatrix which fails det(S)det(S-) >=0:\n");
	  printsubmat(imat, xcombs[r1], ycombs[r2], k, k);
	  flag=0;if((*rkbad)==0){(*rkbad)=k;}
	  if(q){
	    free_imatrix(xcombs, 0, cnk-1, 0, k-1);
	    free_imatrix(ycombs, 0, cmk-1, 0, k-1);
	    free_imatrix(imat2, 0, n1-1, 0, m1-1);
	    free_imatrix(imat3, 0, n1-1, 0, m1-1);
	    return flag;
	  }
	}
	else if(tt==1)// a positive product
	  posprod=1;
      }
    }

    free_imatrix(xcombs, 0, cnk-1, 0, k-1);
    free_imatrix(ycombs, 0, cmk-1, 0, k-1);
  }
  free_imatrix(imat2, 0, n1-1, 0, m1-1);
  free_imatrix(imat3, 0, n1-1, 0, m1-1);

  if(!allm && !posprod)// not r-strongly WSD
    return 0;

  return flag;
}

// working here: this needs to check r-strong compatibility...

/* Check if n X m matrices imat1 and imat2 are compatible */
/* flag = 3 means the matrices are compatible and r-strongly compatible */
/* flag = 2 means the matrices are compatible but not r-strongly compatible */
/* flag = 1 means the matrices are r-strongly compatible but not compatible */
/* flag = 0 means the matrices are none of the above */

int mats_compat0(int **imat1a, int imatrank, int **imat2a, int n1, int m1, int q){
  int k;
  int r,rmin,s1,s2;
  long r1,r2,cnk,cmk;
  int **xcombs;
  int **ycombs;
  int cmpt=1; // compatible
  int rscmpt=1; // r-strongly compatible
  int tt;
  bool posprod=0;

  r=imatrank;
  rmin=1;

  for(k=r;k>=rmin;k--){
    xcombs=allcombsgen(n1,k);
    ycombs=allcombsgen(m1,k);
    cnk=comb(n1,k);cmk=comb(m1,k);
    for(r1=0;r1<cnk;r1++){
      for(r2=0;r2<cmk;r2++){
	//	fprintf(stderr, "in here\n");
	if((tt=fixedminorcompat(imat1a, imat2a, n1, m1, xcombs[r1], ycombs[r2],k)) < 0){
	  fprintf(stderr, "\npair of submatrices with determinants of opposite sign:\n\n");
	  for(s1=0;s1<k;s1++){
	    for(s2=0;s2<k;s2++)
	      fprintf(stderr, "%2d  ",imat1a[xcombs[r1][s1]][ycombs[r2][s2]]);
	    fprintf(stderr, "\n");
	  }
	  fprintf(stderr, "*** and ***\n");
	  for(s1=0;s1<k;s1++){
	    for(s2=0;s2<k;s2++)
	      fprintf(stderr, "%2d  ",imat2a[xcombs[r1][s1]][ycombs[r2][s2]]);
	    fprintf(stderr, "\n");
	  }

	  for(s1=0;s1<k;s1++){
	    fprintf(stderr, "[ ");
	    for(s2=0;s2<k;s2++)
	      fprintf(stderr, "%d,%d  ",xcombs[r1][s1],ycombs[r2][s2]);
	    fprintf(stderr, "]\n");
	  }

	  fprintf(stderr, "_____________________\n\n");
	  if(k==imatrank){rscmpt=0;}else{cmpt=0;} // could return here for speed
	  if(q){
	    free_imatrix(xcombs, 0, cnk-1, 0, k-1);
	    free_imatrix(ycombs, 0, cmk-1, 0, k-1);

	    if(rscmpt==0) // not r-strongly
	      return 0;
	    else
	      return 1;
	  }
	}
	else if(tt>0 && k==imatrank){
	  posprod=1;

	}
      }
    }

    free_imatrix(xcombs, 0, cnk-1, 0, k-1);
    free_imatrix(ycombs, 0, cmk-1, 0, k-1);
    if(k==imatrank && !posprod)
      rscmpt=0;
  }

  if(rscmpt && cmpt)
    return 3;
  else if(rscmpt && !cmpt)
    return 1;
  else if(!rscmpt && cmpt)
    return 2;
  return 0;

}

// are two matrices r-strongly negatively compatible?

int mats_negr_compat(int **imat1a, int imatrank, int **imat2a, int n1, int m1){
  int k, s1,s2;
  long r1,r2,cnk,cmk;
  int **xcombs;
  int **ycombs;
  int tt;
  bool posprod=0;

  k=imatrank;
  xcombs=allcombsgen(n1,k);
  ycombs=allcombsgen(m1,k);
  cnk=comb(n1,k);cmk=comb(m1,k);
  for(r1=0;r1<cnk;r1++){
    for(r2=0;r2<cmk;r2++){
      //	fprintf(stderr, "in here\n");
      if((tt=fixedminorcompat(imat1a, imat2a, n1, m1, xcombs[r1], ycombs[r2],k)) > 0){ // always fail immediately
	fprintf(stderr, "\npair of submatrices with determinants of same sign:\n\n");
	for(s1=0;s1<k;s1++){
	  for(s2=0;s2<k;s2++)
	    fprintf(stderr, "%2d  ",imat1a[xcombs[r1][s1]][ycombs[r2][s2]]);
	  fprintf(stderr, "\n");
	}
	fprintf(stderr, "*** and ***\n");
	for(s1=0;s1<k;s1++){
	  for(s2=0;s2<k;s2++)
	    fprintf(stderr, "%2d  ",imat2a[xcombs[r1][s1]][ycombs[r2][s2]]);
	  fprintf(stderr, "\n");
	}

	for(s1=0;s1<k;s1++){
	  fprintf(stderr, "[ ");
	  for(s2=0;s2<k;s2++)
	    fprintf(stderr, "%d,%d  ",xcombs[r1][s1],ycombs[r2][s2]);
	  fprintf(stderr, "]\n");
	}

	fprintf(stderr, "_____________________\n\n");

	free_imatrix(xcombs, 0, cnk-1, 0, k-1);
	free_imatrix(ycombs, 0, cmk-1, 0, k-1);
	return 0;
      }
      else if(tt<0)
	posprod=1;
    }
  }

  free_imatrix(xcombs, 0, cnk-1, 0, k-1);
  free_imatrix(ycombs, 0, cmk-1, 0, k-1);

  if(posprod)
    return 1;

  return 0;

}

// Wrapper for mats_compat0 and mats_negr_compat

int mats_compat(int **imat1, int **imat2, int n, int m, int q){
  int flg, imatrank;
  int **imat1a, **imat2a;
  int n1,m1;

  imatrank=matrank(imat1,n,m);
  fprintf(stderr, "Checking if (stoich and exp) matrices are compatible, %d-strongly compatible or %d-strongly negatively compatible.\n\n",imatrank,imatrank);
  simppair(imat1, imat2, n, m, &imat1a, &imat2a, &n1, &m1); // imat2a is sparser
  printmat(imat1a, n1, m1);printmat(imat2a, n1, m1);

  flg=mats_compat0(imat1a, imatrank, imat2a, n1, m1, q);
  if(flg){
    free_imatrix(imat1a, 0, n1-1, 0, m1-1);
    free_imatrix(imat2a, 0, n1-1, 0, m1-1);
    if(flg==3){ // both succeeded
      fprintf(stderr, "Finished checking if (stoich and exp) matrices are compatible: matrices are compatible and %d-strongly compatible.\n---------------------------------------\n", imatrank);
      return 3;
    }
    else if(flg==1){
      fprintf(stderr, "Finished checking if (stoich and exp) matrices are compatible: matrices are %d-strongly compatible, but not compatible.\n---------------------------------------\n", imatrank);
      return 1;
    }
    else if(flg==2){
      fprintf(stderr, "Finished checking if (stoich and exp) matrices are compatible: matrices are compatible but not %d-strongly compatible.\n---------------------------------------\n", imatrank);
      return 2;
    }
    return flg;
  }


  flg=mats_negr_compat(imat1a, imatrank, imat2a, n1, m1);
  free_imatrix(imat1a, 0, n1-1, 0, m1-1);
  free_imatrix(imat2a, 0, n1-1, 0, m1-1);
  if(flg){
    fprintf(stderr, "Finished checking if (stoich and exp) matrices are compatible: matrices are %d-strongly negatively compatible.\n---------------------------------------\n", imatrank);
    return -1;
  }
  fprintf(stderr, "Finished checking if (stoich and exp) matrices are compatible: matrices are not compatible or %d-strongly compatible or %d-strongly negatively compatible.\n---------------------------------------\n", imatrank, imatrank);

  return 0;
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

matrix isubmattoexmat(int **A, int *vec1, int n1, int *vec2, int m1){
  matrix C(n1,m1);
  int i,j;
  for(i=0;i<n1;i++){
    for(j=0;j<m1;j++)
      C(i,j)=A[vec1[i]][vec2[j]];
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
  if(m!=(*m1))
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
  if(n!=(*n1))
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

int addintstr(int k, int *s, int len1, int ***t, int *ind)
     /* routine to check if an integer string already belongs to a list of strings, and if not to add it to the end. Returns new free position in the list. all strings have length len1 */
{
  int i=0,j=0, flag;
  (*ind)=0;
  if(k==0){
    (*t) = (int**) malloc(sizeof(int*) * 1);
    (*t)[k] = (int*) malloc(sizeof(int) * (len1));
    for(i=0;i<len1;i++)
      (*t)[k][i] = s[i];
    return k+1;
  }

  flag=0;
  while(j<k){
    //    printvec(s, len1); printvec(s1, len2);printvec((*t)[j], (*t)[j][0]+(*t)[j][1]+2);
    flag=1;
    for(i=0;i<len1;i++){
      if((*t)[j][i]!=s[i]){
	flag=0;
	break;
      }
    }
    if(flag){ // found
      (*ind)=j;
      break;
    }
    j++;
  }

  if(!flag){//not found
    //    fprintf(stderr, "not found: %d\n", k);
    //    printvec(s, len1);printvec(s1,len2);
    (*t)=(int**) realloc((*t), sizeof(int*) *(k+1));
    (*t)[k] = (int*) malloc(sizeof(int) * (len1));
    for(i=0;i<len1;i++)
      (*t)[k][i] = s[i];
    (*ind)=k;
    return k+1;
  }
  return k;
}

int addoneintstr(int k, int *s, int len1, int ***t)
     /* routine to check if an integer string already belongs to a list of strings, and if not to add it to the end. Returns new free position in the list. The first entry of each member of the list is the length of the string*/
{
  int i=0,j=0, flag;

  if(k==0){
    (*t) = (int**) malloc(sizeof(int*) * 1);
    (*t)[k] = (int*) malloc(sizeof(int) * (len1+1));
    (*t)[k][0] = len1;
    for(i=1;i<1+len1;i++)
      (*t)[k][i] = s[i-1];
    return k+1;
  }

  flag=0;
  while(j<k){
    //    printvec(s, len1); printvec(s1, len2);printvec((*t)[j], (*t)[j][0]+(*t)[j][1]+2);
    flag=1;
    if(len1!=(*t)[j][0]){
      j++;flag=0;continue;
    }
    for(i=1;i<1+len1;i++){
      if((*t)[j][i]!=s[i-1]){
	flag=0;
	break;
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
    (*t)[k] = (int*) malloc(sizeof(int) * (len1+1));
    (*t)[k][0] = len1;
    for(i=1;i<1+len1;i++)
      (*t)[k][i] = s[i-1];
    return k+1;
  }
  return k;
}

// Is the ordered int list A a subset of the ordered int list B?

bool AsubsB(int *A, int szA, int *B, int szB){
  int lastfree=0;
  int i,j;
  bool flg;

  if(szA > szB)
    return 0;

  //  printvec(A,szA);  printvec(B,szB);
  for(i=0;i<szA;i++){
    if(lastfree>=szB)// off the end of superset
      return 0;
    for(j=lastfree;j<szB;j++){
      if(B[j]>A[i])
	return 0;
      if(B[j]==A[i]){
	lastfree=j+1;
	flg=1;
	break;
      }
    }
    if(!flg)
      return 0;
  }
  return 1;
}

// Check if the set subs of size n1 is a siphon
// by examining the irreversible stoichiometric matrix
// and the matrix of powers

bool issiphon(int **S,int **M,int n,int m,int *subs,int n1){
  int i,j,k;
  int **tmp=cpmat(S,n,m);
  for(i=0;i<n1;i++){// each  member of the siphon
    for(j=0;j<m;j++){
      if(M[subs[i]][j]!=0){//substrate affects reac rate
	//set column of S to zero
	for(k=0;k<n;k++)
	  tmp[k][j]=0;
      }
    }
  }
  for(i=0;i<n1;i++){
    for(j=0;j<m;j++){
      if(tmp[subs[i]][j]!=0){// row not wiped out
	free_imatrix(tmp, 0, n, 0, m-1);
	return 0;
      }
    }
  }

  free_imatrix(tmp, 0, n, 0, m-1);
  return 1;
}

int checksiphons(int **mat1, int **mat2, int n, int m, int ***allsiphons, int *totminsiphons, int ***allminsiphons){

  int j,k;
  bool flg=0;
  int **xcombs;
  long i,cnk;
  int totsiphons=0;
  (*totminsiphons)=0;

  for(k=1;k<n;k++){// each dimension
    xcombs=allcombsgen(n,k);
    cnk=comb(n,k);
    for(i=0;i<cnk;i++){//each subset

      if(issiphon(mat1,mat2,n,m,xcombs[i],k)){
	// add to minimal list if poss.
	flg=1;
	for(j=0;j<totsiphons;j++){//each existing siphon
	  if(AsubsB((*allsiphons)[j]+1,(*allsiphons)[j][0],xcombs[i],k)){
	    flg=0;break;//not minimal
	  }
	}
	if(flg)
	  (*totminsiphons)=addoneintstr((*totminsiphons),xcombs[i],k,allminsiphons);
	// add to siphon list
	totsiphons=addoneintstr(totsiphons,xcombs[i],k,allsiphons);
      }
    }
    free_imatrix(xcombs, 0, cnk-1, 0, k-1);
  }


  /* fprintf(stderr, "totsiphons=%d\n", totsiphons); */
  /* for(j=0;j<totsiphons;j++){//each siphon */
  /*   for(k=1;k<=(*allsiphons)[j][0];k++) */
  /*     fprintf(stderr, "%d ", (*allsiphons)[j][k]); */
  /*   fprintf(stderr, "\n"); */
  /* } */

  /* fprintf(stderr, "totminsiphons=%d\n", *totminsiphons); */
  /* for(j=0;j<(*totminsiphons);j++){//each siphon */
  /*   for(k=1;k<=(*allminsiphons)[j][0];k++) */
  /*     fprintf(stderr, "%d ", (*allminsiphons)[j][k]); */
  /*   fprintf(stderr, "\n"); */
  /* } */
  return totsiphons;


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
    return sgn(imat[i1[0]][i2[0]]);
  for (i=0;i<dim;i++){
    if(imat[i1[0]][i2[i]]!=0){
      r=0;
      for(j=0;j<dim;j++){
	if(j!=i)
	  i2a[r++]=i2[j];
      }
      tot=qualp(tot,qualt(par(i),qualt(sgn(imat[i1[0]][i2[i]]),qualdetsubmat(imat, n, m, i1+1, i2a, dim-1))));
    }

  }
  //fprintf(stderr, "qualdet=%d\n", tot);

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

void free_imat(int **A, int numA){
  int i;
  for(i=0;i<numA;i++)
    free ((char *)(A[i]));
  if(A)
    free((char *) A);
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
    exit(0);
  }

  while((c=getc(fd)) != EOF)
    i++;
 
  i++;
  str = (char*) malloc(sizeof(char) * (i));
  // second parse

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
  if(k>n)
    return 0;
  for(i=0;i<n-k;i++)
    inds[i]=1;
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

// imat1 is the stoichiometric matrix
// imat2 is the pattern matrix of V
// imat3 is the irreversible stoichiometric matrix
// imat4 is the corresponding matrix of powers (assuming mass action)
// stoichl is the left stoichiometric matrix
// stoichr is the right stoichiometric matrix

int getallreacs(char *str, int ***imat1, int ***imat2, int ***imat3, int ***imat4, int ***stoichl, int ***stoichr, char ***chems, bool *haszerocomplex, int *n, int *m, int *cols3, int *allrev, int *allgood){

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

      for(i=0;i<numleft;i++)
	totchems= addv1(totchems, leftchems[i], chems);
      for(i=0;i<numright;i++)
	totchems= addv1(totchems, rightchems[i], chems);

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
  (*stoichl)=imatrix(0, (*n)-1, 0, irrcols-1);
  (*stoichr)=imatrix(0, (*n)-1, 0, irrcols-1);


  //pass two

  pos=0;numlines=0;
  line = getlinefromstr(&pos, str);
  while(strcmp(line, "")!=0){
    if(!iscomline(line) && isreac(line)){
      getreac(line, &leftchems, &leftstoics, &numleft, &rightchems, &rightstoics, &numright, &rev);
      if(numleft==0||numright==0)
	(*haszerocomplex)=1;

      for(i=0;i<(*n);i++){//each substrate
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

	(*stoichl)[i][numcols]=lstoic;
	(*stoichr)[i][numcols]=rstoic;
	if(rev){
	  (*stoichl)[i][numcols+1]=rstoic;
	  (*stoichr)[i][numcols+1]=lstoic;
	}

	if(lstoic>0 && rstoic > 0){ // on both sides of reac
	  (*allgood)=0;
	  (*imat1)[i][numlines]=rstoic-lstoic;
	  if(rev==1){
	    (*imat2)[i][numlines]=2;
	    (*imat3)[i][numcols]=rstoic-lstoic;
	    (*imat3)[i][numcols+1]=-rstoic+lstoic;
	    (*imat4)[i][numcols]=-lstoic;
	    (*imat4)[i][numcols+1]=-rstoic;

	    //	    (*imat4)[i][numcols+1]=-rightstoics[ind];
	  }
	  else{
	    (*imat2)[i][numlines]=-1;
	    (*imat3)[i][numcols]=rstoic-lstoic;
	    (*imat4)[i][numcols]=-lstoic;
	  }

	}
	else if(lstoic>0){ // only on left
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
	else if(rstoic>0){ // only on right
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
  int j,sp;
  char **v1;
  char **tmp;
  int num1;
  int flag=1;
  if((p=strstr(str, "<-->")) || (p=strstr(str, "<==>")) || (p=strstr(str, "<->"))|| (p=strstr(str, "<=>")) || (p=strstr(str, "<>"))){
    if(strstr(p,"<-->") || strstr(p,"<==>"))
      sp=4;
    else if(strstr(p,"<->") || strstr(p,"<=>"))
      sp=3;
    else
      sp=2;
    j=strlen(str)-strlen(p);
    right=strdup(str+j+sp);
    left=strdup(str);
    left[j]=0;
    (*rev)=1;
  }
  else if((p=strstr(str, "-->")) || (p=strstr(str, "==>")) || (p=strstr(str, "->")) || (p=strstr(str, "=>")) || (p=strstr(str, ">"))){
    if(strstr(p,"-->") || strstr(p,"==>"))
      sp=3;
    else if(strstr(p,"->") || strstr(p,"=>"))
      sp=2;
    else
      sp=1;

    j=strlen(str)-strlen(p);
    right=strdup(str+j+sp);
    left=strdup(str);
    left[j]=0;
    (*rev)=0;
  }
  else if((p=strstr(str, "<--")) || (p=strstr(str, "<==")) || (p=strstr(str, "<-")) || (p=strstr(str, "<=")) || (p=strstr(str, "<"))){
    if(strstr(p,"<--") || strstr(p,"<=="))
      sp=3;
    else if(strstr(p,"<-") || strstr(p,"<="))
      sp=2;
    else
      sp=1;

    j=strlen(str)-strlen(p);
    left=strdup(str+j+sp);
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


  (*numleft)=chemgts2a(left, &v1, '+');
  (*leftstoics)=(int*) malloc(sizeof(int) * (*numleft));
  (*leftchems)=(char**) malloc(sizeof(char*) * (*numleft));
  for(j=0;j<(*numleft);j++){
    num1=chemgts2a(v1[j], &tmp, ' ');

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

  (*numright)=chemgts2a(right, &v1, '+');
  (*rightstoics)=(int*) malloc(sizeof(int) * (*numright));
  (*rightchems)=(char**) malloc(sizeof(char*) * (*numright));
  for(j=0;j<(*numright);j++){
    num1=chemgts2a(v1[j], &tmp, ' ');  
  
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

int chemgts2a(char *s, char ***v, char sep){
  // sep is the separator
  int i, j, k;
  int numgets=0;
  char *tmp;
  i=0, k=0;

  (*v)=NULL;
  while(s[k]){
    j=0;
    while((s[k] == sep) || isspace((int) s[k])){k++;} // skip white space
    while((s[k] != sep) && !isend(s[k])){j++;k++;}
    if(j>0){
      tmp=lrtrim(strchop2(s, k-j, j));
      if(!ispureint(tmp) || atoi(tmp)!=0)
	numgets++;
      free(tmp);
    }
  }
  if(numgets>0){
    (*v)=(char**) malloc(sizeof(char*) * numgets);
    i=0;k=0;
    while(s[k]){
      j=0;
      while((s[k] == sep) || isspace((int) s[k])){k++;} // skip white space
      while((s[k] != sep) && !isend(s[k])){j++;k++;}
      if(j>0){
	tmp=lrtrim(strchop2(s, k-j, j));
	if(!ispureint(tmp) || atoi(tmp)!=0)
	  (*v)[i++] = lrtrim(strchop2(s, k-j, j));
	free(tmp);
      }
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

// For each siphon:
// does there exist a nonnegative vector orthogonal to the siphon face 
// and in ker(Gamma^t)?
// In the language of Angeli et. al., Petri nets paper, we check
// whether each siphon contains the support of a P-semiflow.
// 

bool structpersist(int **imat3, int nlen, int cols3, int **allminsiphons, int totminsiphons){
  int k,j,l,m,siphlen;
  int *ia, *ja;
  double *ar,z;
  int r;
  int **siphmat;
  bool goodflag=1;
  glp_prob *lp;
  glp_smcp parm;

  glp_init_smcp(&parm);
  parm.msg_lev = GLP_MSG_OFF;

  ia=(int *)malloc((size_t) ((1+(nlen+1)*(nlen+cols3+1))*sizeof(int)));
  ja=(int *)malloc((size_t) ((1+(nlen+1)*(nlen+cols3+1))*sizeof(int)));
  ar=(double *)malloc((size_t) ((1+(nlen+1)*(nlen+cols3+1))*sizeof(double)));
  for(k=0;k<totminsiphons;k++){

    siphlen=nlen-allminsiphons[k][0];
    //fprintf(stderr, "siphlen=%d\n", siphlen);
    //make a siphon matrix
    siphmat=imatrix(0,nlen-1,0,siphlen-1);
    for(j=0;j<nlen;j++)//initialise
      for(l=0;l<siphlen;l++)
	siphmat[j][l]=0;
    r=0;

    for(j=0;j<nlen;j++){
      if(!isinlist(j,allminsiphons[k]+1,allminsiphons[k][0])){
	siphmat[j][r++]=1;
      }
    }
    /* printmat(siphmat, nlen, siphlen); */

    lp= glp_create_prob();
    glp_set_obj_dir(lp, GLP_MAX);
    //number of constraints (rows) = col dimension of Gamma + boundedness constraint + siphlen
    glp_add_rows(lp, 1+cols3+siphlen);
    glp_set_row_bnds(lp, 1, GLP_UP, 0.0, 10.0);
    for(j=2;j<cols3+2;j++)
      glp_set_row_bnds(lp, j, GLP_FX, 0.0, 0.0);
    for(j=cols3+2;j<cols3+2+siphlen;j++)
      glp_set_row_bnds(lp, j, GLP_FX, 0.0, 0.0);

    //number of variables (cols) = row dimension of Gamma
    glp_add_cols(lp, nlen);
    for(j=1;j<nlen+1;j++){
      glp_set_col_bnds(lp, j, GLP_LO, 0.0, 0.0);
      glp_set_obj_coef(lp, j, 1.0); // objective function
    }

    m=1;
    for(j=1;j<nlen+1;j++){//col indices
      ia[m]=1;ja[m]=j;ar[m]=1.0;//boundedness constraint
      m++;
      for(l=2;l<2+cols3;l++){//row indices
	ia[m]=l;ja[m]=j;ar[m]=imat3[j-1][l-2];
	m++;
      }
      for(l=2+cols3;l<2+cols3+siphlen;l++){//row indices
	ia[m]=l;ja[m]=j;ar[m]=siphmat[j-1][l-2-cols3];// face constraints
	m++;
      }
    }
    //fprintf(stderr, "m=%d\n", m);
    glp_load_matrix(lp, m-1, ia, ja, ar);
    glp_simplex(lp, &parm);//can be NULL
    z = glp_get_obj_val(lp);
    /* fprintf(stderr, "z = %.2f\n", z); */
 
    glp_delete_prob(lp);
    free_imatrix(siphmat, 0,nlen-1,0,siphlen-1);
    if(z<0.01){// assume fail
      goodflag=0;
      break;
    }
  }
  free ((char *)(ia));free ((char *)(ja));free ((char *)(ar));

  return goodflag;

}

bool has_nonz_signed_col(int **imat, int n, int m){
  int i,j,p;
  for(j=0;j<m;j++){// each column
    p=0;// nothing found yet
    for(i=0;i<n;i++){
      if(imat[i][j]){
	if(!p && imat[i][j]<0)
	  p=-1;
	else if(!p && imat[i][j]>0)
	  p=1;
	else if(p*imat[i][j]<0){
	  p=0;break;//mixed column
	}
      }
    }
    if(p)//not a mixed column
      return 1;
  }
  return 0;
}


// Does Gamma^t have a strictly positive vector in its kernel? 
// Checked only up to some tolerance

int hasposkervec(int **imat1, int nlen, int mlen, bool strict){
  int j,l,m,lpstat;
  int *ia, *ja;
  double *ar,z;
  int goodflag=1;
  glp_prob *lp;
  glp_smcp parm;
  double tol=0.0;
  if(strict){
    if(has_nonz_signed_col(imat1, nlen, mlen))
      return 0;
    tol=0.01;
  }
  glp_init_smcp(&parm);
  parm.msg_lev = GLP_MSG_OFF;

  ia=(int *)malloc((size_t) ((1+(nlen+1)*(mlen+1))*sizeof(int)));
  ja=(int *)malloc((size_t) ((1+(nlen+1)*(mlen+1))*sizeof(int)));
  ar=(double *)malloc((size_t) ((1+(nlen)*(mlen+1))*sizeof(double)));

  lp= glp_create_prob();
  glp_set_obj_dir(lp, GLP_MAX);
  //number of constraints (rows) = col dimension of Gamma + boundedness constraint
  glp_add_rows(lp, 1+mlen);
  glp_set_row_bnds(lp, 1, GLP_UP, 0.0, 10.0);
  for(j=2;j<mlen+2;j++)
    glp_set_row_bnds(lp, j, GLP_FX, 0.0, 0.0);

  //number of variables (cols) = row dimension of Gamma
  glp_add_cols(lp, nlen);m=1;
  for(j=1;j<nlen+1;j++){//col indices
    glp_set_col_bnds(lp, j, GLP_LO, tol, 0.0);// to ensure positivity: could introduce false negatives. 
    glp_set_obj_coef(lp, j, 1.0); // objective function

    ia[m]=1;ja[m]=j;ar[m]=1.0;//boundedness constraint
    m++;
    for(l=2;l<2+mlen;l++){//row indices
      ia[m]=l;ja[m]=j;ar[m]=imat1[j-1][l-2];
      m++;
    }
  }

  glp_load_matrix(lp, m-1, ia, ja, ar);

  //glp_write_lp(lp, NULL,"tmpfile.glp");

  glp_simplex(lp, &parm);//can be NULL
  /* fprintf(stderr, "**%d, %d\n", glp_get_status(lp), GLP_OPT); */
  lpstat=glp_get_status(lp);
  if(strict){// probably no positive vector
    if(lpstat!=GLP_OPT && lpstat!=GLP_FEAS)// infeasible
      goodflag=-1;//probably (unless a very marginally positive vector)!
  }
  else{// no nonnegative vector (must have z value 10 if it exists)
    if(lpstat==GLP_OPT && (z=glp_get_obj_val(lp))<0.01) // solution is zero
      goodflag=0;
  }
  /* z = glp_get_obj_val(lp); */
  /* fprintf(stderr, "z = %.2f\n", z); */
 
  glp_delete_prob(lp);
  free ((char *)(ia));free ((char *)(ja));free ((char *)(ar));

  return goodflag;
}

int colsum(int **imat1, int nlen, int mlen, int j){
  int i, tot=0;
  if(j>=mlen){
    fprintf(stderr, "ERROR in colsum: j out of range. Exiting.\n");
    exit(0);
  }
  for(i=0;i<nlen;i++)
    tot+=imat1[i][j];
  return tot;
}

// Does imat1 have a nonnegative vector in its image?
// if yes, stoichiometry classes unbounded...

int hasposimvec(int **imat1, int nlen, int mlen){
  int j,l,m,lpstat,csum;
  int *ia, *ja;
  double *ar;
  int goodflag=1;
  glp_prob *lp;
  glp_smcp parm;

  glp_init_smcp(&parm);
  parm.msg_lev = GLP_MSG_OFF;

  ia=(int *)malloc((size_t) ((1+(mlen+1)*(nlen+1))*sizeof(int)));
  ja=(int *)malloc((size_t) ((1+(mlen+1)*(nlen+1))*sizeof(int)));
  ar=(double *)malloc((size_t) ((1+(mlen)*(nlen+1))*sizeof(double)));

  lp= glp_create_prob();
  glp_set_obj_dir(lp, GLP_MAX);
  //number of constraints (rows) = row dimension of Gamma + boundedness constraint
  glp_add_rows(lp, 1+nlen);
  glp_set_row_bnds(lp, 1, GLP_DB, 1.0, 10.0);
  for(j=2;j<nlen+2;j++)
    glp_set_row_bnds(lp, j, GLP_LO, 0.0, 0.0);

  //number of variables (cols) = col dimension of Gamma
  glp_add_cols(lp, mlen);
  m=1;
  for(j=1;j<mlen+1;j++){
    glp_set_col_bnds(lp, j, GLP_FR, 0.0, 0.0);// free
    csum=colsum(imat1, nlen, mlen, j-1);
    glp_set_obj_coef(lp, j, csum); // objective function
    ia[m]=1;ja[m]=j;ar[m]=csum;//boundedness constraint
    m++;
    for(l=2;l<2+nlen;l++){//row indices
      ia[m]=l;ja[m]=j;ar[m]=imat1[l-2][j-1];
      m++;
    }
  }

  glp_load_matrix(lp, m-1, ia, ja, ar);

  //glp_write_lp(lp, NULL,"tmpfile1.glp");

  glp_simplex(lp, &parm);//can be NULL
  /* fprintf(stderr, "**%d, %d\n", glp_get_status(lp), GLP_OPT); */
  lpstat=glp_get_status(lp);

  if(lpstat==GLP_NOFEAS) // no feasible solution
    goodflag=0;

  /* z = glp_get_obj_val(lp); */
  /* fprintf(stderr, "z = %.2f\n", z); */
 
  glp_delete_prob(lp);
  free ((char *)(ia));free ((char *)(ja));free ((char *)(ar));

  return goodflag;
}

int **transposemat(int **imat, int n, int m){
  int i,j;
  int **tmp=imatrix(0, m-1, 0, n-1);
  for(i=0;i<n;i++)
    for(j=0;j<m;j++)
      tmp[j][i]=imat[i][j];
  return tmp;
}

int hasposrkervec(int **imat1, int nlen, int mlen, bool strict){
  int flg;
  int **tmp=transposemat(imat1,nlen,mlen);
  flg=hasposkervec(tmp,mlen,nlen,strict);
  free_imatrix(tmp,0,mlen-1,0,nlen-1);
  return flg;
}

int hasposlimvec(int **imat1, int nlen, int mlen){
  int flg;
  int **tmp=transposemat(imat1,nlen,mlen);
  flg=hasposimvec(tmp,mlen,nlen);
  free_imatrix(tmp,0,mlen-1,0,nlen-1);
  return flg;
}


void strongconnect(int i, int **mat, int n, int *inds, int *lowlink, int *ind, int ***S, int *totS, int *Sc, int *Scpos){
  int j,k;

  /* fprintf(stderr, "current: "); */
  /* printvec1(Sc, *Scpos); */
  // Set the depth index for v to the smallest unused index
  inds[i]=*ind;lowlink[i]=*ind;(*ind)++;
  Sc[(*Scpos)++]=i;//add to stack

  // Consider successors of v
  for (j=0;j<n;j++){
    if(mat[i][j]){// each out-edge
      if(inds[j]<0){
	strongconnect(j, mat, n, inds, lowlink, ind, S, totS, Sc, Scpos);
	lowlink[i]=min(lowlink[i],lowlink[j]);
      }
      else if(isinlist(j, Sc, *Scpos))
	lowlink[i]=min(lowlink[i],inds[j]);
    }
  }
  // If node i is a root node, pop the stack and generate an SCC
  if (lowlink[i]==inds[i]){
    //Create SCC
    if((*totS)==0)
      (*S) = (int**) malloc(sizeof(int*) * 1);
    else
      (*S) =(int**) realloc((*S), sizeof(int*) *((*totS)+1));

    k=(*Scpos)-1;
    while(Sc[k]!=i)
      k--;

    (*S)[(*totS)] = (int*) malloc(sizeof(int) * ((*Scpos)-k+1));

    (*S)[(*totS)][0]=(*Scpos)-k;// first element is size
    for(j=k;j<(*Scpos);j++)
      (*S)[(*totS)][j-k+1]=Sc[j];
    (*totS)++;(*Scpos)=k;//reset stack
  }

}

// My implementation of Tarjan's algorithm
// Based on pseudocode on Wikipedia
// mat is the n X n incidence matrix of the digraph
// output is the list of SCCs and totS, the number of SCCs


int **Tarjan(int **mat, int n, int *totS){
  //  output: set of strongly connected components (set of vertex-sets)
  int i, ind=0; //ind increments as we encounter vertices
  int *inds=(int *)malloc((size_t) ((n)*sizeof(int)));
  int *lowlink=(int *)malloc((size_t) ((n)*sizeof(int)));
  int **S=NULL;// to store the SCCs
  int *Sc= (int *)malloc((size_t) ((n)*sizeof(int)));//stack
  int Scpos=0;//stack-counter

  for(i=0;i<n;i++)
    inds[i]=-1;

  for(i=0;i<n;i++){
    if(inds[i]<0)
      strongconnect(i,mat,n,inds,lowlink,&ind, &S,totS, Sc, &Scpos);
  }

  free((char *)(Sc));

  /* fprintf(stderr, "Total number of SCCs in complex graph = %d:\n", *totS); */
  /* for(i=0;i<(*totS);i++){ */
  /*   fprintf(stderr, "i=%d, sz=%d\n", i,S[i][0]); */
  /*   for(ind=0;ind<S[i][0];ind++) */
  /*     fprintf(stderr, "%d ", S[i][1+ind]); */
  /*   fprintf(stderr, "\n"); */
  /* } */
  return S;
}


int isinarray4(int **v, int numv, int *s){
  // is the integer vector s a member of the array v? If so return its index, if not return -1
  int i;
  // assume first entry is length of string
  for(i=0;i<numv;i++){
    if(unordlistareeq(v[i]+1,v[i][0],s+1,s[0]))
      return i;
  }
  return -1;
}

// Confirm that each member of A is also in B
// E.g. the CCs are also SCCs

bool is_sublist(int **A, int numA, int **B, int numB){
  int i;
  for(i=0;i<numA;i++){
    if(isinarray4(B,numB,A[i])<0)
      return 0;
  }
  return 1;
}

// Which member of the partition A is integer i a member of?
// Return -1 if it isn't in the partition
// Assume that the first entry of each integer string in A
// is the length of the remainder

int partmember(int i, int **A, int totA){
  int j,k;
  for(j=0;j<totA;j++){
    for(k=1;k<A[j][0]+1;k++){
      if(i==A[j][k])
	return j;
    }
  }
  return -1;
}

int zerocmplx(int **cmplxs, int totcmplx, int n){
  int i,j;
  bool flg;
  for(i=0;i<totcmplx;i++){
    flg=1;
    for(j=0;j<n;j++){
      if(cmplxs[i][j]){// not zero complex
	flg=0;break;
      }
    }
    if(flg)
      return i;
  }
  return -1;
}


//is b a sublist of a

bool unordsublist(int *a, int *b){
  int i;
  if(b[0]>a[0])
    return 0;
  for(i=1;i<b[0]+1;i++){
    if(!isinlist(b[i],a+1,a[0]))
       return 0;
  }
  return 1;
}

// is a given SCC b in a given CC a? - just check first member as
// each SCC is in a CC iff any element of it is

bool SCCinCC(int *SCCj,int *CCj){
  if(SCCj[0]>CCj[0] || CCj[0]<1 || SCCj[0]<1)
    return 0;
  if(!isinlist(SCCj[1],CCj+1,CCj[0]))
    return 0;
  return 1;
}

//Is an SCC terminal in an CC
//Assume we've checked the SCC is in the CC already

int isterm(int *SCCj,int *CCj,int **cmpmat,int totcmplx){
  int i,j;
  if(SCCj[0]==CCj[0]) // same size: sole terminal SCC in CC
    return 2;
  for(i=1;i<SCCj[0]+1;i++){//each complex in SCC
    for(j=0;j<totcmplx;j++){//outedges from SCCj[i]
      if(cmpmat[SCCj[i]][j] && !isinlist(j,SCCj+1,SCCj[0]))//outedge elsewhere
	return 0;//...so not terminal
    }
  }
  return 1;
}

// Check if the system is weakly reversible
// We do this by computing the complex graph, 
// computing all its SCCs and all its CCs
// and checking that each CC is an SCC
// It could be more efficient to check each
// CC is an SCC as it is computed. 

bool weak_rev(int **imatir, int Srank, int **stoichl, int **stoichr, int n, int m, int *numcomp, int *numlink, bool haszero, bool *zeronotterm, bool *zeroinitial, bool *def1flg, int *deficiency, char **chems, bool q, bool statswitch, bool htmlswitch){
  //stoichl is the left stoichiometric matrix; stoichr is the right stoichiometric matrix
  int i,j,jj,k,m1,totcmplx=0,zeroint;
  int **stoichlt, **stoichrt; //for transposed matrices
  int **cmplxs=NULL;
  int ind1,ind2;
  int **cmpmat; //for SCCs (the digraph)
  int **cmpmat1; //for CCs (the graph)
  int edg[m][2];
  int xc[n];//row list
  int yc[m];//col list
  int **SCC=NULL, **CC=NULL;
  int totSCC=0,totCC=0;
  int tcnt,tflg,CCrnk,totdef,CCdef;
  bool flag=0,onetflg;
  (*zeronotterm)=0;
  (*zeroinitial)=0;

  stoichlt=transposemat(stoichl, n, m);
  stoichrt=transposemat(stoichr, n, m);
  for(i=0;i<m;i++){// get the complexes and the edges of the incidence graph
    totcmplx=addintstr(totcmplx, stoichlt[i], n, &cmplxs, &ind1);
    totcmplx=addintstr(totcmplx, stoichrt[i], n, &cmplxs, &ind2);
    edg[i][0]=ind1;edg[i][1]=ind2;
  }
  (*numcomp)=totcmplx;

  free_imatrix(stoichlt,0,m-1,0,n-1);
  free_imatrix(stoichrt,0,m-1,0,n-1);

  if(!q){
    fprintf(stderr, "Complexes:\n");
    for(i=0;i<totcmplx;i++){
      fprintf(stderr, "{ ");
      for(j=0;j<n;j++){
	if(cmplxs[i][j]==1)
	  fprintf(stderr, "%s ", chems[j]);
	else if(cmplxs[i][j])
	  fprintf(stderr, "%d%s ", cmplxs[i][j],chems[j]);
      }
      fprintf(stderr, "}\n");
    }
    fprintf(stderr, "\n");
  }

  //Create the adjacency matrix of the complex digraph (cmpmat) 
  //and its symmetrised version (cmpmat1)
  cmpmat=imatrix(0,totcmplx-1,0,totcmplx-1);
  cmpmat1=imatrix(0,totcmplx-1,0,totcmplx-1);
  for(i=0;i<totcmplx;i++){
    for(j=0;j<totcmplx;j++){
      cmpmat[i][j]=0;
      cmpmat1[i][j]=0;
    }
  }
  for(i=0;i<m;i++){
    cmpmat[edg[i][0]][edg[i][1]]=1;
    cmpmat1[edg[i][0]][edg[i][1]]=1;
    cmpmat1[edg[i][1]][edg[i][0]]=1;
  }

  // Extract the SCCs and CCs using Tarjan's algorithm
  // First element of a CC or SCC is its size
  // Then the list of complexes
  SCC=Tarjan(cmpmat, totcmplx, &totSCC);
  CC=Tarjan(cmpmat1, totcmplx, &totCC);
  (*numlink)=totCC;

  if(haszero){//if the zero complex is there...
    zeroint=zerocmplx(cmplxs, totcmplx, n); //...find its index
    k=partmember(zeroint,SCC,totSCC); //it is in the kth SCC
    for(jj=0;jj<totcmplx;jj++){// for each complex...
      if(partmember(jj,SCC,totSCC)==k){// ...in the same SCC as zero
	for(j=0;j<totcmplx;j++){ // for each complex...
	  if(cmpmat[jj][j] && (partmember(j,SCC,totSCC)!=k)){// outedge kth SCC
	    (*zeronotterm)=1;break;// so not terminal
	  }
	}
      }
      if((*zeronotterm))
	break;
    }
    //check if zero is initial (so we know whether {0} is an equilibrium)
    for(j=0;j<totcmplx;j++){
      if(j!=zeroint && cmpmat[zeroint][j]){// outedge from zero
	(*zeroinitial)=1;break;//some source reactions
      }
    }

  }

  /* //Which SCCs are terminal? */
  /* term=(int *)malloc((size_t) ((totSCC)*sizeof(int))); */
  /* for(i=0;i<totSCC;i++) */
  /*   term[i]=1; */

  /* for(i=0;i<totcmplx;i++){ */
  /*   k=partmember(i,SCC,totSCC); */
  /*   for(j=0;j<totcmplx;j++){ */
  /*     if(cmpmat[i][j] && (partmember(j,SCC,totSCC)!=k)){// outedge */
  /* 	term[i]=0;break; */
  /*     } */
  /*   } */
  /* } */

  if(!q){
    fprintf(stderr, "Complex incidence matrix:\n");
    printmat(cmpmat,totcmplx,totcmplx);
  }

  /* fprintf(stderr, "terminal?\n"); */
  /* printvec(term,totSCC); */

  if(!q){
    fprintf(stderr, "CCs of complex graph:\n");
    for(i=0;i<totCC;i++){
      fprintf(stderr, "    { ");
      for(j=1;j<CC[i][0]+1;j++)
	fprintf(stderr, "%d ", CC[i][j]);
      fprintf(stderr, "}\n");
    }
    fprintf(stderr, "\n");

    fprintf(stderr, "SCCs of complex graph:\n");
    for(i=0;i<totSCC;i++){
      fprintf(stderr, "    { ");
      for(j=1;j<SCC[i][0]+1;j++)
	fprintf(stderr, "%d ", SCC[i][j]);
      fprintf(stderr, "}\n");
    }
    fprintf(stderr, "\n");
  }

  if(is_sublist(CC,totCC,SCC,totSCC)){
    if(statswitch)
      fprintf(stdout, "WRflag\n");//for batch checking
    flag=1; // weakly reversible
  }

  (*deficiency)=totcmplx-totCC-Srank;//rank of CRN
  //fprintf(stderr, "CRNdef=%d\n", (*deficiency));


  if(statswitch &&!(*deficiency)){
    fprintf(stdout, "DEF0\n");//for batch checking
  }

  if((*deficiency)){//nonzero then does CRN satisfy the def. one thm?
    //Does each CC contain exactly one terminal SCC?
    onetflg=1;//assume all good initially
    for(i=0;i<totCC;i++){//each CC
      tcnt=0; //no of terminal SCCs in CC
      for(j=0;j<totSCC;j++){//each SCC
	if(SCCinCC(SCC[j],CC[i]) && (tflg=isterm(SCC[j],CC[i],cmpmat,totcmplx))){//is in there and is terminal
	  //fprintf(stderr, "%d is terminal in %d, flg=%d\n", j, i,tflg);
	  if(tflg==2){//sole terminal CC
	    tcnt=1;
	    break;
	  }
	  else
	    tcnt++;
	  if(tcnt>=2){
	    onetflg=0;
	    break;
	  }
	}
      }
      if(onetflg==0)
	break;
    }
    //does each subnetwork have deficiency <=1?
    (*def1flg)=0;totdef=0;
    if(onetflg){//Each CC contains exactly one terminal SCC
      (*def1flg)=1;
      for(i=0;i<n;i++){xc[i]=i;}//all rows
      //compute the deficiency of each CC = numcomp-1-Srank;
      for(i=0;i<totCC;i++){//each CC
	m1=0;
	//make the reduced matrix from imatir
	for(jj=0;jj<m;jj++){//each reaction
	  if(isinlist(edg[jj][0],CC[i]+1,CC[i][0]))//only check left complex
	    yc[m1++]=jj;
	}
	CCrnk=submatrank(imatir,xc,n,yc,m1);//rank of subnetwork
	CCdef=CC[i][0]-1-CCrnk;//deficiency of subnetwork
	if(CCdef>=2){
	  (*def1flg)=0;
	  break;//
	}
	totdef+=CCdef;//total deficiency
	//fprintf(stderr, "%d subndef=%d\n", i, CCdef);
      }
    }
    //Is the sum of subnetwork deficiencies equal to the network deficiency?
    if((*def1flg) && totdef!=(*deficiency))
      (*def1flg)=0;
  }
  if(statswitch && (*def1flg)){
    fprintf(stdout, "DEF1\n");//for batch checking
  }

  if((*deficiency) && !(*def1flg)){
    if(htmlswitch)
      fprintf(stdout, "The network is not deficiency zero and fails the conditions of the deficiency one theorem (see <a href=\"http://reaction-networks.net/wiki/Deficiency_theory.\" target=\"_blank\">the wiki</a> for more details.)\n\n");
    else
      fprintf(stdout, "The network is not deficiency zero and fails the conditions of the deficiency one theorem.\n\n");
  }

  //http://reaction-networks.net/wiki/Deficiency_theory#Deficiency_one_theorem

  free_imat(SCC,totSCC);
  free_imat(CC,totCC);
  free_imat(cmplxs,totcmplx);

  free_imatrix(cmpmat,0,totcmplx-1,0,totcmplx-1);
  free_imatrix(cmpmat1,0,totcmplx-1,0,totcmplx-1);

  return flag;
}

// Is the DSR graph of imat1 (nXm) and imat2^t (mXn) strongly connected?

bool DSRSC(int **imat1, int **imat2, int n, int m){
  int i,j,**SCC,**imat=imatrix(0,n+m-1,0,n+m-1);
  int totSCC=0;

  for(i=0;i<n;i++){
    for(j=0;j<n;j++)
      imat[i][j]=0;
    for(j=n;j<n+m;j++)
      imat[i][j]=imat1[i][j-n];
  }

  for(i=n;i<n+m;i++){
    for(j=0;j<n;j++)
      imat[i][j]=imat2[j][i-n];
    for(j=n;j<n+m;j++)
      imat[i][j]=0;
  }

  /* printmat(imat,n+m,n+m); */

  SCC=Tarjan(imat, n+m, &totSCC);
  free_imatrix(imat,0,n+m-1,0,n+m-1);
  free_imat(SCC,totSCC);
  if(totSCC==1)
    return 1;
  return 0;
}

// given a matrix with no more than two nonzero entries in each column
// check if its DSR graph has o-cycles: rows can be re-signed to 
// ensure no column has more than one positive entry or more than one 
// negative entry
// The basic idea follows M. Banaji, Cycle structure in SR and DSR graphs: implications for multiple equilibria and stable oscillation in chemical reaction networks, in K. Jensen, S. Donatelli and J. Kleijn (eds): Transactions On Petri Nets and Other models of Concurrency (ToPNoC), volume V, series: LNCS, volume 6900 (2012) 
// An SR-graph with S-degree <= 2 can be R-sorted if and only if it contains no o-cycles. (If it contains an o-cycle it can't be R-sorted is trivial since this would change the o-cycle to an e-cycle; the other direction for a connected graph is Lemma 2) in the ref above

int hasoloops(int **mat, int *resgn, int n, int m){

  int i, j, k, p, q,flg=1;
  resgn[0]=1;
  for(i=1;i<n;i++)
    resgn[i]=0;

  while(flg){
    flg=0;
    for(j=0;j<m;j++){// each column
      p=0;// no nonzero entries found yet
      for(i=0;i<n;i++){ // each entry
	if(!p && mat[i][j]){// first nonzero entry in col j
	  p=mat[i][j];q=i;
	}
	else if(mat[i][j]*p > 0){ // second entry of same sign in col j
	  if(!(resgn[q]) && !(resgn[i])){
	    resgn[i]=-1;
	    resgn[q]=1;
	    for(k=0;k<m;k++) // resign row i
	      mat[i][k]=-mat[i][k];
	    flg=1;// something resigned
	    break;
	  }
	  else if(!(resgn[i])){
	    resgn[i]=-1;
	    for(k=0;k<m;k++) // resign row i
	      mat[i][k]=-mat[i][k];
	    flg=1; // something resigned
	    break;
	  }
	  else if(!(resgn[q])){
	    resgn[q]=-1;
	    for(k=0;k<m;k++) // resign row q
	      mat[q][k]=-mat[q][k];
	    flg=1;// something resigned
	    break;
	  }
	  else // can't resign: failed
	    return 1;
	}
      }
      if(flg) // some row got resigned, start again. 
	break;
    }
  }// exit loop if no changes

  for(i=0;i<n;i++)//non-resigned rows
    if(!(resgn[i]))
      resgn[i]=1;
  return 0;
}

int hasoloopst(int **mat, int *resgn, int n, int m){
  int **tmp=transposemat(mat,n,m);
  int i,j,flg;
  flg=hasoloops(tmp,resgn,m,n);
  for(i=0;i<n;i++)
    for(j=0;j<m;j++)
      mat[i][j]=tmp[j][i];

  free_imatrix(tmp,0,m-1,0,n-1);
  return flg;
}

// row degree of a matrix

int rowdegree(int **imat, int n, int m){
  int i,j,tmp,deg=0;
  for(i=0;i<n;i++){//each row
    tmp=0;
    for(j=0;j<m;j++){
      if(imat[i][j])
	tmp++;
    }
    deg=max(deg, tmp);
  }
  return deg;
}

// column degree of a matrix

int coldegree(int **imat, int n, int m){
  int i,j,tmp,deg=0;
  for(j=0;j<m;j++){//each column
    tmp=0;
    for(i=0;i<n;i++){
      if(imat[i][j])
	tmp++;
    }
    deg=max(deg, tmp);
  }
  return deg;
}

// multiply two integer matrices A and B of dimensions nXr and rXm
// to get the integer matrix C of dimension nXm

void multindAB(int **A, int **B, int **C, int n, int r, int m){

  int i,j,k;
  for(i=0;i<n;i++){
    for(j=0;j<m;j++){
      C[i][j]=0;
      for(k=0;k<r;k++)
	C[i][j]+=A[i][k]*B[k][j];
    }
  }
  return;
}

int **multindAB1(int **A, int **B, int n, int m){
  //output AB^t where A,B are nXm matrices
  int i,j,k;
  int **C=imatrix(0,n-1,0,n-1);
  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      C[i][j]=0;
      for(k=0;k<m;k++)
	C[i][j]+=A[i][k]*B[j][k];
    }
  }
  return C;
}

// symmetrise a pattern matrix

int **symmetrise(int **A, int n){
  int i,j;
  int **B=imatrix(0, n-1, 0, n-1);
  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
      B[i][j]=0;

  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      if(A[i][j]){
	B[i][j]=1;B[j][i]=1;
      }
    }
  }
  return B;
}

bool connected_prod(int **A, int **B, int n, int m){
  // Is the i-graph of the product of AB^t
  // connected (not necessarily strongly)?
  // Assumes that the I-graph is well-defined
  int **C,**C1,**CC=NULL;
  int totCC=0;
  bool flg=0;
  C=multindAB1(A,B,n,m);
  C1=symmetrise(C,n);
  CC=Tarjan(C1,n,&totCC);
  if(totCC==1)
    flg=1;
  free_imatrix(C,0,n-1,0,n-1);
  free_imatrix(C1,0,n-1,0,n-1);
  free_imat(CC,totCC);
  return flg;
}

// From Graph-theoretic characterizations of monotonicity
// of chemical networks in reaction coordinates
// David Angeli, Patrick De Leenheer and Eduardo Sontag
// JMB (2009)
// Assumes it has been checked that no chemicals appear on both 
// sides of a reaction. 
// Assumes we have checked that stoichiometry classes are bounded
// Assumes structural persistence


int ALS_JMB_2009(int **im1, int **im2, int n, int m, int allgood, int bdclass, int persistflag){
  int *resgn=(int *)malloc((size_t) ((m)*sizeof(int)));
  //local copies
  int **imat1=cpmat(im1, n, m);
  int **imat2=cpmat(im2, n, m);
  int i,j;
  int flg=0;

  if(!allgood || bdclass!=1 || !persistflag)
    return 0;

  //  fprintf(stderr, "rowd = %d\n", rowdegree(imat1, n, m));
  if(rowdegree(imat1, n, m)<=2 && !hasoloopst(imat1, resgn, n, m)){//orthant monotone (gets resigned to preserve the nonnegative orthant)
    //fprintf(stderr, "here\n");
    // resign imat2
    for(i=0;i<n;i++)
      for(j=0;j<m;j++)
	imat2[i][j]=resgn[i]*imat2[i][j];

    // is the I-graph of the reverse product connected?
    if(connected_prod(imat2,imat1,m,n)){// strong monotonicity
      //if(hasposrkervec(imat1,n,m,1)==1)//Condition (9) in ref.
      if(!hasposlimvec(imat1,n,m))//Condition (9) in ref.
	flg=1;
      else if(!hasposrkervec(imat1,n,m,0))//Condition (10) in ref.
	flg=2;
    }
  }
  free_imatrix(imat1,0,n-1,0,m-1);
  free_imatrix(imat2,0,n-1,0,m-1);
  free((char*)resgn);
  return flg;
}


int simple_CRN(int **im, int n, int m){
  int rd,cd,flg=0;
  int **tmp,*resgn;
  int **imat=cpmat(im,n,m);//local copy
  // first for the matrix
  if((rd=rowdegree(imat, n, m))<=2 || (cd=coldegree(imat, n, m))<=2){
    if(rd){
      resgn=(int *)malloc((size_t) ((n)*sizeof(int)));
      if(!hasoloops(imat, resgn, n, m))
	flg=1;
      free((char*)resgn);
    }
    else if(cd){
      tmp=transposemat(imat, n, m);
      resgn=(int *)malloc((size_t) ((m)*sizeof(int)));
      if(!hasoloops(tmp, resgn, m, n))
	flg=2;
      free((char*)resgn);
      free_imatrix(tmp,0,m-1,0,n-1);
    }
    
  }
  free_imatrix(imat,0,n-1,0,m-1);
  return flg;
}

void printsiphons(int **allsiphons, int totsiphons, char **chems){
  int i,j;
  //  fprintf(stderr, "tot=%d\n",totsiphons);
  for(i=0;i<totsiphons;i++){
    fprintf(stderr, "    {");
    for(j=1;j<allsiphons[i][0];j++)
      fprintf(stderr, "%s ", chems[allsiphons[i][j]]);
    fprintf(stderr, "%s}\n", chems[allsiphons[i][j]]);
  }
  fprintf(stderr, "\n");
}

//
//Euclid's algorithm
//

int gcd(int a, int b)
{
  int temp;
  if(a<=0 || b <=0){
    fprintf(stderr, "the routine gcd requires positive integers.\n");
    exit(0);
  }
  while(b){
    temp=a%b;a=b;b=temp;
  }

  return a;
}

// Returns -2 if unresignable
// or there are signed columns of different signs

int hasnegpos(int **imat, int n, int m){
  int i,j,p,q=0;
  for(j=0;j<m;j++){// each column
    p=0;// nothing found yet
    for(i=0;i<n;i++){
      if(imat[i][j]){
	if(!p && imat[i][j]<0)
	  p=-1;
	else if(!p && imat[i][j]>0)
	  p=1;
	else if(p*imat[i][j]<0){
	  p=0;break;//mixed column
	}
      }
    }
    if(p){//not a mixed column
      if(!q && p<0)
	q=-1;
      else if(!q && p>0)
	q=1;
      else if(q*p<0) // has signed columns of different sign
	return -2;
    }
  }
  return q;

}

//gcd of entries in a vector

int gcd_vec(int *v, int n){
  int i,j=0,g;
  if(n==1){
    if(v[0])
      return v[0];
    else
      return 1;
  }
  while(!(v[j]))
    j++;
  if(j>=n)
    return 1;
  g=abs(v[j]);
  for(i=j+1;i<n;i++)
    if(v[i])
      g=gcd(g,abs(v[i]));
  return g;
}

int reduce_vec(int *v, int n){
  int i;
  int g=gcd_vec(v, n);
  for(i=0;i<n;i++)
    v[i]/=g;
  return g;
}

// is a row v2 a nonzero multiple of row v1

int get_mult(int *v1, int *v2, int n){
  int i,t=0;
  for(i=0;i<n;i++){
    if(!(v1[i]*v2[i]) && (v1[i] || v2[i]))
      return 0;
    else if(!t && v1[i] && v2[i])
      t=v2[i]/v1[i];
    if(v2[i]!=t*v1[i])
      return 0;
  }
  return t;
}

void addrow(int ***mat, int n, int m){
  if(!n || !(*mat))
    (*mat) = (int**) malloc(sizeof(int*) * 1);
  else
    (*mat) =(int**) realloc((*mat), sizeof(int*)*(n+1));
  (*mat)[n]=(int*) malloc(sizeof(int) * m);
}

// Carry out the factorisation
// partition rows into pairwise linearly dependent sets
// for each such set create a single row in the second factor
// and a column in the first factor
// output the new "middle" dimension and the two matrices

int special_fact(int **mat, int n, int m, int ***F1, int ***F2o){
  int **F2=NULL,**F1t=NULL;
  int i,j,g,t,r=0;
  int dlt[n];
  /* cout << "entering special_fact.\n"; */
  for(i=0;i<n;i++)//initially all rows unused
    dlt[i]=0;

  for(i=0;i<n;i++){//each row
  if(!(dlt[i])){//not yet used: new row and column
    addrow(&F2,r,m);addrow(&F1t,r,n);r++;dlt[i]=1;
      for(j=0;j<m;j++)
	F2[r-1][j]=mat[i][j];
      g=reduce_vec(F2[r-1], m);
      for(j=0;j<i;j++)// set initial entries to zero (definitely all dealt with)
	F1t[r-1][j]=0;
      F1t[r-1][i]=g; // set entry to g
      for(j=i+1;j<n;j++){
	if(dlt[j])//used
	  F1t[r-1][j]=0;
	else if((t=get_mult(F2[r-1],mat[j],m))){//multiple of reduced?
	  F1t[r-1][j]=t;dlt[j]=1;//row now dealt with
	}
	else
	  F1t[r-1][j]=0;
      }
    }
  }
  /* printmat(mat,n,m); */
  /* printmat(F2,r,m); */
  /* printmat(F1t,r,n); */

  (*F2o)=cpmat(F2,r,m);
  (*F1)=transposemat(F1t,r,n);
  free_imat(F2,r);free_imat(F1t,r);
  /* cout << "exiting special_fact.\n"; */
  return r;
}

bool matseq(int **A, int **B, int n, int m){
  int i,j;
  for(i=0;i<n;i++)
    for(j=0;j<m;j++)
      if(A[i][j]!=B[i][j])
	return 0;
  return 1;
}


// Routine to check if a matrix has a PM1 (Pete-Murad-1) 
// factorisation or not. The first factor must be a matrix 
// with exactly one nonzero entry in each row. The second 
// must have left nullspace of dimension 1 and including 
// a strictly positive vector

int hasPM1factors(int **Gamma, int **Vpat, int n, int m, int allgood, int q){
  int i,j,h;
  int **Lambda=NULL, **Theta=NULL;
  int cd,rdim=0;
  int *resgn;
  int **Gamma_new;
  bool genflg=0,scflg=0,Kflg=0;
  if(!allgood)
    return 0;

  rdim=special_fact(Gamma,n,m,&Lambda,&Theta);
  resgn =(int *)malloc((size_t) ((rdim)*sizeof(int)));

  /* printmat(Lambda,n,rdim); */
  /* printmat(Theta,rdim,m); */


  if((cd=coldegree(Theta, rdim, m))<=2 && !hasoloops(Theta, resgn, rdim, m) && matrank(Theta, rdim, m)==rdim-1){//also resigns the second factor
    genflg=1;
    if(!q){
      fprintf(stderr, "This seems to be a monotone system.\n\n");
      //printvec(resgn,rdim);
    }
    //resign the first factor
    for(j=0;j<rdim;j++){//each column
      for(i=0;i<n;i++)//each entry
	Lambda[i][j]=resgn[j]*Lambda[i][j];
    }
    if(!q){
      printmat(Lambda,n,rdim);
      printmat(Theta,rdim,m);
    }
    // if it helps, reverse signs of Lambda and Theta

    if((h=hasnegpos(Lambda,n,rdim))==-1){
      for(j=0;j<rdim;j++){
	for(i=0;i<n;i++)
	  Lambda[i][j]=-Lambda[i][j];
	for(i=0;i<m;i++)
	  Theta[j][i]=-Theta[j][i];
      }
    }
    if(h!=-2)
      Kflg=1;

    //sanity check
    Gamma_new=imatrix(0,n-1,0,m-1);
    multindAB(Lambda,Theta,Gamma_new,n,rdim,m);
    if(!matseq(Gamma,Gamma_new,n,m)){
      fprintf(stderr, "something went wrong: matrix product doesn't work out.\n");
      printmat(Gamma, n, m);
      fprintf(stderr, "doesn't factor as:\n\n");
      printmat(Lambda, n, rdim);printmat(Theta, rdim, m);
      exit(0);
    }
    free_imatrix(Gamma_new,0,n-1,0,m-1);

    if(!q){
      fprintf(stderr, "The stoichiometric matrix:\n\n");
      printmat(Gamma, n, m);
      fprintf(stderr, "factors into:\n\n");
      printmat(Lambda, n, rdim);
      printmat(Theta, rdim, m);
    }

    if(!q){
      if(h!=-2)
	fprintf(stderr, "K(Lambda) does not include negative vectors.\n\n");
      else
	fprintf(stderr, "However K(Lambda) includes negative vectors.\n\n");
    }
    scflg=DSRSC(Gamma, Vpat, n, m);
    if(!q){
      if(scflg)      
	fprintf(stderr, "The DSR graph is strongly connected.\n\n");
      else
	fprintf(stderr, "The DSR graph appears not to be strongly connected.\n\n");
    }
  }


  free((char *)resgn);
  free_imatrix(Lambda, 0,n-1,0,rdim-1);
  free_imatrix(Theta,0,rdim-1,0,m-1);

  if(genflg && scflg && Kflg)
    return 1;
  if(genflg && scflg)//not necessarily bounded SCs
    return 2;

  return 0;
}


int analysereacs(const char fname[], int q, bool htmlswitch, bool statswitch){

  int k,**imat1, **imat2, **imat3, **imat4, **imat1a, **stoichl, **stoichr;
  char *str;
  int mlen=0, nlen=0, mlena=0;
  char **chems;
  long pos=0;
  char *line;
  int type=-1;
  int allrev, allgood;
  int cols3;
  int csdflag=0;
  int ssdflag=0;
  int wsdflag=0;
  int compatflag=0;
  int MAcompatflag=0;
  char mpnestr[100];
  char IC1str[300];
  char IC1ppstr[300];
  char IC1pstr[300];
  char IC3str[300];
  char IC13str[500];
  char IC1pp3str[500];
  char IC1p3str[500];
  char MAIC3[500];
  char MAIC2[500];
  char MAIC2p[500];
  char MAIC2pp[500];
  char MAIC2IC3[500];
  char MAIC2pIC3[500];
  char MAIC2ppIC3[500];
  char notSSD[500];
  char notWSD[500];
  char notrWSD[500];
  char notrcmpt[500];
  char notrcmpt1[500];
  char feinbergdef0[500];
  char panteadef0[500];
  char feinbergdef1[500];
  char ALSstr[500];
  char PMglob[500];
  char PMloc[500];
  char BanajiPanteastr[200];
  char deftheor[200];
  char weakrstr[200];
  char defstr[200];
  char siphonstr[200];
  char structpers[200];
  char bddclass[300];
  char unbddclass[300];
  int totsiphons=1,totminsiphons=0,**allsiphons=NULL,**allminsiphons=NULL;
  bool persistflag=1,def1flg=0;
  int bdclass, notallbd;
  int Srank,numlink,numcomp,deficiency,simp_flg,ALS_flg,PM1_flg;
  bool weakr, haszerocomplex=0,zeronotterm,zeroinitial,SSPO=1,posSSPO=1;
  char *tmpstr;

  //SauroSplit((char *)("datfiles/d_d2s3r3.txt"), (char *)("datfiles/s3r3/s3r3"));
  //return 0;


  //Version x=year 2012+x, .y=month number, .z = revision number
  fprintf(stdout, "Analysereacs version 1.12.3. (Please note that this is work in progress.)\n\n");

  str=readfileintostr(fname);
  if(isonlyspace(str)){
   fprintf(stderr, "\n        **ERROR**\nNothing found in file \"%s\". EXITING.\n          *****\n", fname);
    free(str);
    return -1;
  }



  if(htmlswitch){
    strcpy(BanajiPanteastr, "<a href=\"http://arxiv.org/abs/1309.6771\" target=\"_blank\">Banaji and Pantea, arxiv:1309.6771</a>");
    strcpy(bddclass, "<a href=\"http://reaction-networks.net/wiki/CoNtRol#Stoichiometric_subspace_and_stoichiometry_classes\" target=\"_blank\">Stoichiometry classes are bounded</a>");
    strcpy(unbddclass, "<a href=\"http://reaction-networks.net/wiki/CoNtRol#Stoichiometric_subspace_and_stoichiometry_classes\" target=\"_blank\">Stoichiometry classes are unbounded</a> (the stoichiometric subspace includes a nonnegative vector)");
    strcpy(structpers, "<a href=\"http://reaction-networks.net/wiki/CoNtRol#Persistence_condition_2_.28PC2.29\" target=\"_blank\">structurally persistent</a>");
    strcpy(siphonstr, "<a href=\"http://reaction-networks.net/wiki/CoNtRol#Siphon\" target=\"_blank\">siphons</a>");
    strcpy(defstr, "<a href=\"http://reaction-networks.net/wiki/Deficiency#Deficiency\" target=\"_blank\">deficiency</a>");
    strcpy(weakrstr, "<a href=\"http://reaction-networks.net/wiki/Weakly_reversible#Weak_reversibility\" target=\"_blank\">weakly reversible</a>");
    strcpy(deftheor, "<a href=\"http://reaction-networks.net/wiki/Deficiency_theory\" target=\"_blank\">deficiency theory</a>");
    strcpy(feinbergdef0, "Theorem 6.1.1 in Feinberg (<a href=\"http://www.sciencedirect.com/science/article/pii/0009250987800994\" target=\"_blank\">Chem. Eng. Sci. 42(10), 1987</a>)");
    strcpy(panteadef0, "Theorem 6.3 in Pantea (<a href=\"http://epubs.siam.org/doi/abs/10.1137/110840509\" target=\"_blank\">SIAM J. Math. Anal. 44(3), 2012</a>)");
    strcpy(feinbergdef1, "Theorem 6.2.1 in Feinberg (<a href=\"http://www.sciencedirect.com/science/article/pii/0009250987800994\" target=\"_blank\">Chem. Eng. Sci. 42(10), 1987</a>");
  strcpy(ALSstr, "Theorem 2 in Angeli, De Leenheer and Sontag (<a href=\"http://link.springer.com/article/10.1007/s00285-009-0309-0\" target=\"_blank\">J. Math. Biol. 61(4), 2010</a>)");
  strcpy(PMglob, "Theorem 2.2 in Donnell and Banaji (<a href=\"http://epubs.siam.org/doi/abs/10.1137/120898486\" target=\"_blank\">SIADS, 12(2), 2013</a>)");
  strcpy(PMloc, "Theorem 2.1 in Donnell and Banaji (<a href=\"http://epubs.siam.org/doi/abs/10.1137/120898486\" target=\"_blank\">SIADS, 12(2), 2013</a>)");
  strcpy(mpnestr, "<a title=\"MPNE\" href=\"http://reaction-networks.net/wiki/CoNtRol#MPNE\">MPNE</a>");
  strcpy(IC1ppstr, "General kinetics: no stoichiometry class includes more than one equilibrium. The system satisfies condition <a title=\"IC1++\" href=\"http://reaction-networks.net/wiki/CoNtRol#Injectivity_condition_1.2B.2B_.28IC1.2B.2B.29\">IC1++</a>");
    strcpy(IC1pstr, "General kinetics: no nontrivial stoichiometry class includes more than one equilibrium. The system satisfies condition <a title=\"IC1+\" href=\"http://reaction-networks.net/wiki/CoNtRol#Injectivity_condition_1.2B_.28IC1.2B.29\">IC1+</a>");
    strcpy(IC1str, "General kinetics: no stoichiometry class includes more than one positive equilibrium. The system satisfies condition <a title=\"IC1\" href=\"http://reaction-networks.net/wiki/CoNtRol#Injectivity_condition_1_.28IC1.29\">IC1</a>");
    strcpy(IC3str, "General kinetics: the fully open system is injective. The system satisfies condition <a title=\"IC3\" href=\"http://reaction-networks.net/wiki/CoNtRol#Injectivity_condition_3_.28IC3.29\">IC3</a>");
    strcpy(IC13str, "General kinetics: no stoichiometry class includes more than one positive equilibrium, and the fully open system is injective. The system satisfies conditions <a title=\"IC1\" href=\"http://reaction-networks.net/wiki/CoNtRol#Injectivity_condition_1_.28IC1.29\">IC1</a> and <a title=\"IC3\" href=\"http://reaction-networks.net/wiki/CoNtRol#Injectivity_condition_3_.28IC3.29\">IC3</a>");
    strcpy(IC1p3str, "General kinetics: no nontrivial stoichiometry class includes more than one equilibrium, and the fully open system is injective. The system satisfies conditions <a title=\"IC1+\" href=\"http://reaction-networks.net/wiki/CoNtRol#Injectivity_condition_1.2B_.28IC1.2B.29\">IC1+</a> and <a title=\"IC3\" href=\"http://reaction-networks.net/wiki/CoNtRol#Injectivity_condition_3_.28IC3.29\">IC3</a>");
    strcpy(IC1pp3str, "General kinetics: no stoichiometry class includes more than one equilibrium, and the fully open system is injective. The system satisfies conditions <a title=\"IC1++\" href=\"http://reaction-networks.net/wiki/CoNtRol#Injectivity_condition_1.2B.2B_.28IC1.2B.2B.29\">IC1++</a> and <a title=\"IC3\" href=\"http://reaction-networks.net/wiki/CoNtRol#Injectivity_condition_3_.28IC3.29\">IC3</a>");
    strcpy(MAIC2, "Mass action kinetics: the system satisfies condition <a title=\"IC2\" href=\"http://reaction-networks.net/wiki/CoNtRol#Injectivity_condition_2_.28IC2.29\">IC2</a>");
    strcpy(MAIC2p, "Mass action kinetics: no nontrivial stoichiometry class includes boundary equilibria or has more than one equilibrium. The system satisfies condition <a title=\"IC2+\" href=\"http://reaction-networks.net/wiki/CoNtRol#Injectivity_condition_2.2B_.28IC2.2B.29\">IC2+</a>");
    strcpy(MAIC2pp, "Mass action kinetics: the system has no nonzero boundary equilibria and no stoichiometry class includes more than one equilibrium. The system satisfies condition <a title=\"IC2++\" href=\"http://reaction-networks.net/wiki/CoNtRol#Injectivity_condition_2.2B.2B_.28IC2.2B.2B.29\">IC2++</a>");
    strcpy(MAIC3, "Mass action kinetics: the fully open system is injective. The system satisfies condition <a title=\"IC3\" href=\"http://reaction-networks.net/wiki/CoNtRol#Injectivity_condition_3_.28IC3.29\">IC3</a>");
    strcpy(MAIC2IC3, "Mass action kinetics: no stoichiometry class includes more than one positive equilibrium, and the fully open system is injective. The system satisfies conditions <a title=\"IC2\" href=\"http://reaction-networks.net/wiki/CoNtRol#Injectivity_condition_2_.28IC2.29\">IC2</a> and <a title=\"IC3\" href=\"http://reaction-networks.net/wiki/CoNtRol#Injectivity_condition_3_.28IC3.29\">IC3</a>");
    strcpy(MAIC2pIC3, "Mass action kinetics: no nontrivial stoichiometry class has boundary equilibria or includes more than one equilibrium, and the fully open system is injective. The system satisfies conditions <a title=\"IC2+\" href=\"http://reaction-networks.net/wiki/CoNtRol#Injectivity_condition_2.2B_.28IC2.2B.29\">IC2+</a> and <a title=\"IC3\" href=\"http://reaction-networks.net/wiki/CoNtRol#Injectivity_condition_3_.28IC3.29\">IC3</a>");
    strcpy(MAIC2ppIC3, "Mass action kinetics: the system has no nonzero boundary equilibria, no stoichiometry class includes more than one equilibrium, and the fully open system is injective. The system satisfies conditions <a title=\"IC2++\" href=\"http://reaction-networks.net/wiki/CoNtRol#Injectivity_condition_2.2B.2B_.28IC2.2B.2B.29\">IC2++</a> and <a title=\"IC3\" href=\"http://reaction-networks.net/wiki/CoNtRol#Injectivity_condition_3_.28IC3.29\">IC3</a>");
    strcpy(notSSD, "There exists a choice of power-law kinetics and inflows and outflows such that the fully open system has multiple positive equilibria");
    strcpy(notrWSD, "Mass action kinetics: the system fails condition <a title=\"IC2\" href=\"http://reaction-networks.net/wiki/CoNtRol#Injectivity_condition_2_.28IC2.29\">IC2</a> for some choice of rate constants");
    strcpy(notWSD, "Mass action kinetics: there exists a choice of rate constants and inflows and outflows such that the system fails condition <a title=\"IC3\" href=\"http://reaction-networks.net/wiki/CoNtRol#Injectivity_condition_3_.28IC3.29\">IC3</a> (the fully open system is noninjective)");
    strcpy(notrcmpt, "There exists a choice of power-law kinetics such that the system fails condition <a title=\"IC2\" href=\"http://reaction-networks.net/wiki/CoNtRol#Injectivity_condition_2_.28IC2.29\">IC2</a>");
    strcpy(notrcmpt1, "There exists a choice of power-law kinetics such that the system has multiple positive equilibria on some stoichiometry class");
}
  else{
    strcpy(BanajiPanteastr, "Banaji and Pantea, arxiv:1309.6771");
    strcpy(bddclass, "Stoichiometry classes are bounded");
    strcpy(unbddclass, "Stoichiometry classes are unbounded (the stoichiometric subspace includes a nonnegative vector)");
    strcpy(structpers, "structurally persistent");
    strcpy(siphonstr, "siphons");
    strcpy(defstr, "deficiency");
    strcpy(weakrstr, "weakly reversible");
    strcpy(deftheor, "deficiency theory");
    strcpy(feinbergdef0, "Theorem 6.1.1 in Feinberg (Chem. Eng. Sci. 42(10), 1987)");
    strcpy(panteadef0, "Theorem 6.3 in Pantea (SIAM J. Math. Anal. 44(3), 2012)");
    strcpy(feinbergdef1, "Theorem 6.2.1 in Feinberg (Chem. Eng. Sci. 42(10), 1987)");
    strcpy(ALSstr, "Theorem 2 in Angeli, De Leenheer and Sontag (J. Math. Biol. 61(4), 2010)");
    strcpy(PMglob, "Theorem 2.2 in Donnell and Banaji (SIADS, 12(2), 2013)");
    strcpy(PMloc, "Theorem 2.1 in Donnell and Banaji (SIADS, 12(2), 2013)");

    strcpy(mpnestr, "MPNE");
    strcpy(IC1ppstr, "General kinetics: no stoichiometry class includes more than one equilibrium. The system satisfies condition IC1++");
    strcpy(IC1pstr, "General kinetics: no nontrivial stoichiometry class includes more than one equilibrium. The system satisfies condition IC1+");
    strcpy(IC1str, "General kinetics: no stoichiometry class includes more than one positive equilibrium. The system satisfies condition IC1");
    strcpy(IC3str, "General kinetics: the fully open system is injective. The system satisfies condition IC3");
    strcpy(IC13str, "General kinetics: no stoichiometry class includes more than one positive equilibrium, and the fully open system is injective. The system satisfies conditions IC1 and IC3");
    strcpy(IC1pp3str, "General kinetics: no stoichiometry class includes more than one equilibrium, and the fully open system is injective. The system satisfies conditions IC1++ and IC3");
    strcpy(IC1p3str, "General kinetics: no nontrivial stoichiometry class includes more than one equilibrium, and the fully open system is injective. The system satisfies conditions IC1+ and IC3");
    strcpy(MAIC2, "Mass action kinetics: the system satisfies condition IC2");
    strcpy(MAIC2p, "Mass action kinetics: no nontrivial stoichiometry class includes boundary equilibria or has more than one equilibrium. The system satisfies condition IC2+");
    strcpy(MAIC2pp, "Mass action kinetics: the system has no nonzero boundary equilibria and no stoichiometry class includes more than one equilibrium. The system satisfies condition IC2++");
    strcpy(MAIC3, "Mass action kinetics: the fully open system is injective. The system satisfies condition IC3");
    strcpy(MAIC2IC3, "Mass action kinetics: no stoichiometry class includes more than one positive equilibrium, and the fully open system is injective. The system satisfies conditions IC2 and IC3");
   strcpy(MAIC2pIC3, "Mass action kinetics: no nontrivial stoichiometry class has boundary equilibria or includes more than one equilibrium, and the fully open system is injective. The system satisfies conditions IC2+ and IC3");
   strcpy(MAIC2ppIC3, "Mass action kinetics: the system has no nonzero boundary equilibria, no stoichiometry class includes more than one equilibrium, and the fully open system is injective. The system satisfies conditions IC2++ and IC3");
   strcpy(notSSD, "There exists a choice of power-law kinetics and inflows and outflows such that the fully open system has multiple positive equilibria");
    strcpy(notrWSD, "Mass action kinetics: the system fails condition IC2 for some choice of rate constants");
    strcpy(notWSD, "Mass action kinetics: there exists a choice of rate constants and inflows and outflows such that the system fails condition IC3 (the fully open system is noninjective)");
    strcpy(notrcmpt, "There exists a choice of power-law kinetics such that the system fails condition IC2");
    strcpy(notrcmpt1, "There exists a choice of power-law kinetics such that the system has multiple positive equilibria on some stoichiometry class");

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
    if(!getallreacs(str, &imat1, &imat2, &imat3, &imat4, &stoichl, &stoichr, &chems, &haszerocomplex, &nlen, &mlen, &cols3, &allrev, &allgood)){
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
      free_imatrix(stoichl, 0, nlen-1, 0, cols3-1);
      free_imatrix(stoichr, 0, nlen-1, 0, cols3-1);
      freearraydat(chems, nlen);
      return -1;
    }
    free(str);

    Srank=matrank(imat1, nlen, mlen);
    fprintf(stderr, "The stoichiometric matrix (rank = %d):\n\n", Srank);
    printmat1(chems, imat1, nlen, mlen);
    fprintf(stderr, "The pattern matrix for -V^T:\n\n");
    printmatpat(chems, imat2, nlen, mlen);
    fprintf(stderr, "The irreversible stoichiometric matrix:\n\n");
    printmat1(chems, imat3, nlen, cols3);
    fprintf(stderr, "The matrix of powers (mass action):\n\n");
    printmat1(chems, imat4, nlen, cols3);
    fprintf(stderr, "The left stoichiometric matrix:\n\n");
    printmat1(chems, stoichl, nlen, cols3);
    fprintf(stderr, "The right stoichiometric matrix:\n\n");
    printmat1(chems, stoichr, nlen, cols3);

    fprintf(stderr, "\n_________________________________\n\n");

    // Is the system weakly reversible?
    weakr=weak_rev(imat3, Srank, stoichl, stoichr, nlen, cols3, &numcomp, &numlink, haszerocomplex, &zeronotterm, &zeroinitial, &def1flg, &deficiency, chems, q, statswitch, htmlswitch);
    if(deficiency==0){
      if(weakr){
	if(Srank>3){
	  fprintf(stdout, "This is a %s %s zero network. According to %s, with mass-action kinetics: each nontrivial stoichiometry class admits exactly one positive equilibrium, and this equilibrium is locally asymptotically stable relative to its stoichiometry class. There are no positive, nontrivial periodic orbits.\n\n", weakrstr, defstr, feinbergdef0);
	}
	else{
	  fprintf(stdout, "This is a %s %s zero network with stoichiometric subspace of dimension %d. With mass-action kinetics: (i) According to %s, each nontrivial stoichiometry class admits exactly one positive equilibrium which attracts the relative interior of its stoichiometry class; (ii) By %s, this equilibrium is also locally asymptotically stable relative to its stoichiometry class. \n\n", weakrstr, defstr, Srank, panteadef0, feinbergdef0);
	}
	strcat(notrWSD, ". By ");strcat(notrWSD, deftheor);
	strcat(notrWSD, " however, with mass action kinetics, the system has precisely one positive equilibrium on each nontrivial stoichiometry class");
      }
      else if(zeronotterm){
	SSPO=0;//no steady states of POs
	fprintf(stdout, "This is a %s zero network and the zero complex does not lie in a terminal strong linkage class. According to %s, for general kinetics there are no equilibria at all and no nontrivial periodic orbits, including on the boundary.\n\n", defstr, feinbergdef0);
	strcat(notrcmpt, ". By deficiency theory however, there are no equilibria at all and no nontrivial periodic orbits, including on the boundary (for general kinetics)");
	strcat(notrWSD, ". By deficiency theory however, there are no equilibria at all and no nontrivial periodic orbits, including on the boundary (for general kinetics)");
      }
      //https://reaction-networks.net/wiki/Reaction_graph#Deficiency
      else{
	posSSPO=0;// no positive steady states or POs
	fprintf(stdout, "The network has %s zero, but is not %s. According to %s, for general kinetics there are no positive equilibria or positive nontrivial periodic orbits.\n\n", defstr, weakrstr, feinbergdef0);
	strcat(notrcmpt, ". By deficiency theory however, there are no positive equilibria at all (for general kinetics) - see ");
	strcat(notrcmpt, feinbergdef0);
	strcat(notrWSD, ". By deficiency theory however, there are no positive equilibria at all (for general kinetics) - see ");
	strcat(notrWSD, feinbergdef0);
      }
    }
    else if(def1flg){
      if(weakr){
	fprintf(stdout, "This network is %s and satisfies the conditions of the deficiency one theorem. According to %s, for mass action kinetics, the system has precisely one positive equilibrium on each nontrivial stoichiometry class.\n\n", weakrstr, feinbergdef1);
	strcat(notrWSD, ". By ");strcat(notrWSD, deftheor);
	strcat(notrWSD, " however, with mass action kinetics, the system has precisely one positive equilibrium on each nontrivial stoichiometry class");

      }
      else{
	fprintf(stdout, "This network satisfies the conditions of the deficiency one theorem (but is not %s). According to %s, for mass action kinetics, the system has no more than one positive equilibrium on each nontrivial stoichiometry class.\n\n", weakrstr, feinbergdef1);
	strcat(notrWSD, ". By deficiency theory however, with mass action kinetics, the system has no more than one positive equilibrium on each nontrivial stoichiometry class - see ");
	strcat(notrWSD, feinbergdef1);
      }
    }
    else{
      if(weakr)
	fprintf(stdout, "The network is %s.\n\n", weakrstr);
      else
	fprintf(stdout, "The network is not %s.\n\n", weakrstr);

      fprintf(stdout, "The network has %s %d.\n\n", defstr, deficiency);
    }

    if((simp_flg=simple_CRN(imat1,nlen,mlen))){
      if(simp_flg==1)
	fprintf(stderr, "The system has no o-cycles and S-degree <=2\n\n");
      else
	fprintf(stderr, "The system has no o-cycles and R-degree <=2\n\n");
    }

    /* fprintf(stderr, "numcomp = %d\n", numcomp); */
    /* fprintf(stderr, "numlink = %d\n", numlink); */
    /* fprintf(stderr, "Srank = %d\n", Srank); */

    // check if the stoichiometric matrix has positive vectors in its left-kernel
    bdclass=hasposkervec(imat1, nlen, mlen, 1);//strict
    if(bdclass==-1 && hasposimvec(imat1, nlen, mlen))
      bdclass=0;//definitely unbounded

    if(!bdclass)
      fprintf(stdout, "%s.\n\n",unbddclass);
    else if(bdclass==-1)//probably unbounded (shouldn't occur)
      fprintf(stdout, "Stoichiometry classes may be unbounded.\n\n");
    else
      fprintf(stdout, "%s (each stoichiometry class contains at least one equilibrium).\n\n", bddclass);

    // test if the irreversible stoichiometric matrix has no positive vectors in its kernel
    notallbd=hasposrkervec(imat3, nlen, cols3, 1);//strict
    if(notallbd==-1 && hasposlimvec(imat3, nlen, cols3))
      notallbd=0;

    //count the siphons
    totsiphons=checksiphons(imat3,imat4,nlen,cols3,&allsiphons,&totminsiphons,&allminsiphons);
    if(statswitch && !totsiphons)
      fprintf(stdout, "SIPH0\n");
    if(SSPO && posSSPO){ //deficiency zero tests didn't already answer these questions
      if(!notallbd && !totsiphons){
	if(!zeroinitial)//zero is an equilibrium
	  fprintf(stdout, "The only equilibrium of the system is a trivial one (the system has no siphons, and no positive equilibria).\n\n");
	else
	  fprintf(stdout, "The system has no equilibria (the system has no siphons, no positive equilibria, and zero is not an equilbrium).\n\n");
      }
      /* else if(notallbd==-1 && !totsiphons){//shouldn't execute */
      /* 	if(!zeroinitial)//zero is an equilibrium */
      /* 	  fprintf(stdout, "The only equilibrium of the system appears to be the trivial one.\n\n"); */
      /* 	else */
      /* 	  fprintf(stdout, "The system appears to have no equilibria.\n\n"); */
      /* } */
      else if(!notallbd)
	fprintf(stdout, "The system has no positive equilibria.\n\n");
      /* else if(notallbd==-1)//shouldn't execute */
      /* 	fprintf(stdout, "All equilibria appear to be boundary equilibria.\n\n"); */
      else if(!totsiphons)
	fprintf(stdout, "This system has no %s (the boundary includes no nontrivial equilbria).\n\n", siphonstr);
    }

    if(totsiphons){
      fprintf(stderr, "%s of the system:\n", siphonstr);
	printsiphons(allsiphons, totsiphons, chems);
	fprintf(stderr, "Minimal %s of the system:\n", siphonstr);
      printsiphons(allminsiphons, totminsiphons, chems);

      persistflag=structpersist(imat3, nlen, cols3, allminsiphons, totminsiphons);
      if(persistflag)
	fprintf(stdout, "The system has %s, but is %s: no nontrivial stoichiometry class includes boundary equilibria.\n\n", siphonstr, structpers);
    }

    if((ALS_flg=ALS_JMB_2009(imat1,imat2,nlen,mlen,allgood,bdclass,persistflag))){
      if(ALS_flg==1){
	if(statswitch)
	  fprintf(stdout, "ALS1 (global convergence)\n");
	fprintf(stdout, "According to %s, the system (with general kinetics) is globally convergent in the following sense: all positive initial conditions converge to an equilibrium which is the unique equilibrium on its stoichiometry class.\n\n", ALSstr);
      }
      else if(ALS_flg==2){
	if(statswitch)
	  fprintf(stdout, "ALS2 (generic QC)\n");
	fprintf(stdout, "According to %s, the system (with general kinetics) is generically quasiconvergent: almost all positive initial conditions converge to the set of equilibria (the measure of the set of possibly non-convergent initial conditions is zero).\n\n", ALSstr);
      }
    }

    if((PM1_flg=hasPM1factors(imat1, imat2, nlen, mlen, allgood, 1))){
      if(PM1_flg==1){
	if(!totsiphons){
	  if(statswitch)
	    fprintf(stdout, "PMflag1 (global no siph)\n");
	  fprintf(stdout, "According to %s, the system (with general kinetics) is globally convergent in the following sense: all initial conditions converge to an equilibrium which is the unique equilibrium on its stoichiometry class, and is (for stoichiometry classes other than {0}) positive.\n\n", PMglob);
	}
	else if(persistflag){
	  if(statswitch)
	    fprintf(stdout, "PMflag2 (global with siph)\n");
	  fprintf(stdout, "According to %s, the system (with general kinetics) is globally convergent in the following sense: all initial conditions on any nontrivial stoichiometry class converge to an equilibrium which is positive and is the unique equilibrium on its stoichiometry class.\n\n", PMglob);
	}
      }
      else if (PM1_flg==2){
	if(statswitch)
	  fprintf(stdout, "PMflag3 (local)\n");
	fprintf(stdout, "According to %s, the system (with general kinetics) is locally convergent in the following sense: every positive initial condition is locally asymptotically stable.\n\n", PMloc);
      }
    }


    fprintf(stdout, "Claims based on theory in %s:\n\n", BanajiPanteastr);
    //free siphons
    for(k=0;k<totsiphons;k++)
      free ((char *)(allsiphons[k]));
    if(allsiphons)
      free((char *) allsiphons);
    for(k=0;k<totminsiphons;k++)
      free ((char *)(allminsiphons[k]));
    if(allminsiphons)
      free((char *) allminsiphons);

    if(allgood && allrev){ // no reactants on both sides and all reversible
      imat1a=redmat(imat1, nlen, mlen, &mlena);//remove redundant cols
      csdflag=isCSD(imat1a, nlen, mlena, q);
      if(csdflag!=2){
	ssdflag=isSSD(imat1a, nlen, mlena, q);
	if(ssdflag!=2)
	  wsdflag=doubleisWSD(imat1a, nlen, mlena, q);
      }
      if(csdflag==2){//CSD
	if(totsiphons==0)
	  fprintf(stdout, "%s. (In fact, this reaction structure satisfies conditions IC1++ and IC3 with *any stoichiometries*.)\n", IC1pp3str);
	else if(persistflag)
	  fprintf(stdout, "%s. (In fact, this reaction structure satisfies conditions IC1+ and IC3 with *any stoichiometries*.)\n", IC1p3str);
	else
	  fprintf(stdout, "%s. (In fact, this reaction structure satisfies conditions IC1 and IC3 with *any stoichiometries*.)\n", IC13str);
      }
      else if(ssdflag){
	if(ssdflag==2 && csdflag==1){//r-CSD but not CSD
	  if(totsiphons==0)
	    fprintf(stdout, "%s. (In fact, this reaction structure satisfies conditions IC1++ with *any stoichiometries*.)\n", IC1pp3str);
	  else if(persistflag)
	    fprintf(stdout, "%s. (In fact, this reaction structure satisfies conditions IC1+ with *any stoichiometries*.)\n", IC1p3str);
	  else
	    fprintf(stdout, "%s. (In fact, this reaction structure satisfies conditions IC1 with *any stoichiometries*.)\n", IC13str);
	}
	else if(ssdflag==2){//SSD but not CSD
	  if(totsiphons==0)
	    fprintf(stdout, "%s.\n", IC1pp3str);
	  else if(persistflag)
	    fprintf(stdout, "%s.\n", IC1p3str);
	  else
	    fprintf(stdout, "%s.\n", IC13str);
	}
	else if(csdflag==1){//r-CSD (and r-SSD) but not CSD or SSD
	  if(totsiphons==0)
	    fprintf(stdout, "%s. (In fact, this reaction structure satisfies condition IC1++ with *any stoichiometries*.) %s.\n", IC1ppstr, notSSD);
	  else if(persistflag)
	    fprintf(stdout, "%s. (In fact, this reaction structure satisfies condition IC1+ with *any stoichiometries*.) %s.\n", IC1pstr, notSSD);
	  else
	    fprintf(stdout, "%s. (In fact, this reaction structure satisfies condition IC1 with *any stoichiometries*.) %s.\n", IC1str, notSSD);
	}
	else{//r-SSD (but not SSD)
	  if(totsiphons==0)
	    fprintf(stdout, "%s. %s.\n", IC1ppstr, notSSD);
	  else if(persistflag)
	    fprintf(stdout, "%s. %s.\n", IC1pstr, notSSD);
	  else
	    fprintf(stdout, "%s. %s.\n", IC1str, notSSD);
	}
      }
      else{// neither SSD nor r-SSD (as reversible, nonzero vector in ker(Gamma_ir) is automatic)
	fprintf(stdout, "General kinetics: %s. %s.\n", notrcmpt1,notSSD);
      }

      if(csdflag!=2 && ssdflag!=2){ //Not strongest conclusions: check mass action
	if(wsdflag==3){ //WSD and r-strongly WSD
	  if(totsiphons==0)
	    fprintf(stdout, "%s.\n", MAIC2ppIC3);
	  else if(persistflag)
	    fprintf(stdout, "%s.\n", MAIC2pIC3);
	  else
	    fprintf(stdout, "%s.\n", MAIC2IC3);
	}
	else if(wsdflag==2) //WSD but not r-strongly WSD
	  fprintf(stdout, "%s.\n%s.\n", MAIC3, notrWSD);
	else if(wsdflag==1){ //r-strongly WSD but not WSD
	  if(totsiphons==0)
	    fprintf(stdout, "%s.\n%s.\n", MAIC2pp, notWSD);
	  else if(persistflag)
	    fprintf(stdout, "%s.\n%s.\n", MAIC2p, notWSD);
	  else
	    fprintf(stdout, "%s.\n%s.\n", MAIC2, notWSD);
	}
      }
      free_imatrix(imat1a, 0, nlen-1, 0, mlena-1);
      if(csdflag==0 && ssdflag==0 && wsdflag==0)
	fprintf(stdout, "%s.\n%s.\n", notrWSD,notWSD);
    }//finished with the simply reversible case

    if(!allgood || (!allrev && csdflag==0 && ssdflag==0 && wsdflag==0)){
  
      /* flag = 3 means the matrices are compatible and r-strongly compatible */
      /* flag = 2 means the matrices are compatible but not r-strongly compatible */
      /* flag = 1 means the matrices are r-strongly compatible but not compatible */
      /* flag = 0 means the matrices are none of the above */

      compatflag=arecompat(imat1, imat2, nlen, mlen, q);
      if(compatflag==3){//matrix and sign-pattern are compatible and r-strongly compatible
	if(statswitch)
	  fprintf(stdout, "NR***\nGKIC1\nGKIC3\nMAIC2\nMAIC3\n");//for batch checking
	if(totsiphons==0)
	  fprintf(stdout, "%s.\n", IC1pp3str);
	else if(persistflag)
	  fprintf(stdout, "%s.\n", IC1p3str);
	else
	  fprintf(stdout, "%s.\n", IC13str);
      }
      else if(compatflag==2){ //matrix and sign-pattern are compatible but not r-strongly compatible
	if(statswitch)
	  fprintf(stdout, "NR***\nGKIC3\nMAIC3\n");//for batch checking
	if(SSPO && posSSPO && notallbd)//possibly interior equilibria
	  tmpstr=notrcmpt1;
	else
	  tmpstr=notrcmpt;

	MAcompatflag=mats_compat(imat3, imat4, nlen, cols3, q);
	if(MAcompatflag==3 || MAcompatflag==1 || MAcompatflag==-1){// (stoich and exp) matrices are r-strongly compatible or r-strongly negatively compatible
	  if(statswitch)
	    fprintf(stdout, "MAIC2\n");//for batch checking
	  if(totsiphons==0)
	    fprintf(stdout, "%s.\n%s.\n%s.\n", IC3str, tmpstr, MAIC2pp);
	  else if(persistflag)
	    fprintf(stdout, "%s.\n%s.\n%s.\n", IC3str, tmpstr, MAIC2p);
	  else
	    fprintf(stdout, "%s.\n%s.\n%s.\n", IC3str, tmpstr, MAIC2);
	}
	else{
	  if(SSPO && posSSPO && notallbd)//possibly interior equilibria
	    fprintf(stdout, "%s.\n%s.\n%s.\n", IC3str,notrWSD,notrcmpt1);
	  else
	    fprintf(stdout, "%s.\n%s.\n", IC3str,notrWSD);
	}
      }
      else if(compatflag==1 || compatflag==-1){//  matrix and sign-pattern are r-strongly compatible or r-strongly negatively compatible but are not compatible
	if(statswitch)
	  fprintf(stdout, "NR***\nGKIC1\nMAIC2\n");//for batch checking
	MAcompatflag=mats_compat(imat3, imat4, nlen, cols3, q);
	if(MAcompatflag==3 || MAcompatflag==2){// (stoich and exp) matrices are compatible
	  if(statswitch)
	    fprintf(stdout, "MAIC3\n");//for batch checking
	  if(totsiphons==0)
	    fprintf(stdout, "%s.\n%s.\n%s.\n", IC1ppstr, notSSD, MAIC3);
	  else if(persistflag)
	    fprintf(stdout, "%s.\n%s.\n%s.\n", IC1pstr, notSSD, MAIC3);
	  else
	    fprintf(stdout, "%s.\n%s.\n%s.\n", IC1str, notSSD, MAIC3);
	}
	else//not compatible
	  if(totsiphons==0)
	    fprintf(stdout, "%s.\n%s.\n", IC1ppstr, notWSD);
	  else if(persistflag)
	    fprintf(stdout, "%s.\n%s.\n", IC1pstr, notWSD);
	  else
	    fprintf(stdout, "%s.\n%s.\n", IC1str, notWSD);
      }
      else{//check MA

	if(SSPO && posSSPO && notallbd)//possibly interior equilibria
	  tmpstr=notrcmpt1;
	else
	  tmpstr=notrcmpt;

	MAcompatflag=mats_compat(imat3, imat4, nlen, cols3, q);
	if(MAcompatflag==3){// (stoich and exp) matrices are compatible and r-strongly compatible
	  if(statswitch)
	    fprintf(stdout, "NR***\nMAIC2\nMAIC3\n");//for batch checking
	  if(totsiphons==0)
	    fprintf(stdout, "%s.\n%s.\n%s.\n", notSSD, tmpstr, MAIC2ppIC3);
	  else if(persistflag)
	    fprintf(stdout, "%s.\n%s.\n%s.\n", notSSD, tmpstr, MAIC2pIC3);
	  else
	    fprintf(stdout, "%s.\n%s.\n%s.\n", notSSD, tmpstr, MAIC2IC3);
	}
	else if(MAcompatflag==2){ //(stoich and exp) matrices are compatible but not r-strongly compatible or r-strongly negatively compatible
	  if(statswitch)
	    fprintf(stdout, "NR***\nMAIC3\n");//for batch checking
	  fprintf(stdout, "%s.\n%s.\n%s.\n%s.\n", notSSD, tmpstr, MAIC3, notrWSD);
	}
	else if(MAcompatflag==1 || MAcompatflag==-1){ //(stoich and exp) matrices are r-strongly compatible or r-strongly negatively compatible, but not compatible
	  if(statswitch)
	    fprintf(stdout, "NR***\nMAIC2\n");//for batch checking
	  if(totsiphons==0)
	    fprintf(stdout, "%s.\n%s.\n%s.\n", tmpstr, MAIC2pp, notWSD);
	  else if(persistflag)
	    fprintf(stdout, "%s.\n%s.\n%s.\n", tmpstr, MAIC2p, notWSD);
	  else
	    fprintf(stdout, "%s.\n%s.\n%s.\n", tmpstr, MAIC2, notWSD);
	}
	else{
	  if(statswitch)
	    fprintf(stdout, "NR***\n");//for batch checking
	  //	allminorsigns(imat3, imat4, nlen, cols3, q);
	  if(SSPO && posSSPO && notallbd)//possibly interior equilibria
	    fprintf(stdout, "%s.\n%s.\n%s.\n", notrcmpt1, notrWSD, notWSD);
	  else
	    fprintf(stdout, "%s.\n%s.\n", notrWSD, notWSD);
	}
      }
    }
 

    free_imatrix(imat1, 0, nlen-1, 0, mlen-1);
    free_imatrix(imat2, 0, nlen-1, 0, mlen-1);
    free_imatrix(imat3, 0, nlen-1, 0, cols3-1);
    free_imatrix(imat4, 0, nlen-1, 0, cols3-1);
    free_imatrix(stoichl, 0, nlen-1, 0, cols3-1);
    free_imatrix(stoichr, 0, nlen-1, 0, cols3-1);
    freearraydat(chems, nlen);
  }
  //
  // The file contains two matrices
  //
  else if (type==1){
    fprintf(stderr, "The current version of this routine doesn't support this format. Please enter the reactions in reaction format.\n");
    /* if(!readmatpairfromstr(str, &nlen, &mlen, &imat1, &imat2)){  */
    /*   fprintf(stderr, "ERROR: Expecting two matrices in file \"%s\". Couldn't find these. EXITING. \n", fname); */
    /*   free(str); */
    /*   return -1; */
    /* } */
    /* fprintf(stderr, "Assuming that these are the stoichiometric matrix and the pattern matrix for -V^T:\n\n"); */
    /* printmat(imat1, nlen, mlen); */
    /* printmat(imat2, nlen, mlen); */

    /* compatflag=arecompat(imat1, imat2, nlen, mlen, q); */
    /* if(compatflag==3) */
    /*   fprintf(stdout, "%s.\n", IC13str); */
    /* else if(compatflag==2) */
    /*   fprintf(stdout, "%s.\n", IC3str); */
    /* else if(compatflag==1) */
    /*   fprintf(stdout, "%s\n", IC1str); */
    /* else */
    /*   fprintf(stdout, "No injectivity claims seem possible for mass-action kinetics or for general kinetics.\n"); */
    /* free_imatrix(imat1, 0, nlen-1, 0, mlen-1); */
    /* free_imatrix(imat2, 0, nlen-1, 0, mlen-1); */
  }
  //
  // The file contains a single matrix
  //
  else{ // single matrix, assumed to be a stoichiometric matrix
    fprintf(stderr, "The current version of this routine doesn't support this format. Please enter the reactions in reaction format.\n");
    /* imat1=readmatrixfromstr(str, &nlen, &mlen); */
    /* if(!(*imat1)){ */
    /*   fprintf(stderr, "ERROR: Couldn't find reactions, a pair of matrices or a matrix in file \"%s\". EXITING. \n", fname); */
    /*   free(str); */
    /*   free_imatrix(imat1, 0, nlen-1, 0, mlen-1); */
    /*   return -1; */
    /* } */
    /* imat1a=redmat(imat1, nlen, mlen, &mlena);//remove redundant cols */

    /* if(mlena!=mlen){ */
    /*   fprintf(stderr, "WARNING: This matrix appears to have redundant columns. Assuming that this is a stoichiometric matrix (rank = %d), that no reactants appear on both sides of a reaction, and that all reactions are reversible. If you aren't happy with these assumptions, use the reaction format rather than the matrix format. Will use the following matrix for analysis:\n",matrank(imat1, nlen, mlen)); */
    /*   printmat(imat1a, nlen, mlena); */
    /* } */
    /* else{ */
    /*   fprintf(stderr, "WARNING: Assuming that this is a stoichiometric matrix (rank = %d), that no reactants appear on both sides of a reaction, and that all reactions are reversible. If you aren't happy with these assumptions, use the reaction format rather than the matrix format.\n",matrank(imat1, nlen, mlen)); */
    /*   printmat(imat1a, nlen, mlena); */
    /* } */

    /* csdflag=isCSD(imat1a, nlen, mlena, q); */
    /* if(csdflag!=2){ */
    /*   ssdflag=isSSD(imat1a, nlen, mlena, q); */
    /*   if(ssdflag!=2) */
    /* 	wsdflag=doubleisWSD(imat1a, nlen, mlena, q); */
    /* } */
    /* if(csdflag==2) */
    /*   fprintf(stdout, "%s. (In fact, this reaction structure satisfies conditions IC1 and IC3 with *any stoichiometries*.)\n", IC13str); */
    /* else if(ssdflag){ */
    /*   if(ssdflag==2 && csdflag==1) */
    /* 	fprintf(stdout, "%s. (In fact, this reaction structure satisfies condition IC1, but not necessarily IC3, with *any stoichiometries*.)\n", IC13str); */
    /*   else if(ssdflag==2) */
    /* 	fprintf(stdout, "%s.\n", IC13str); */
    /*   else if(csdflag==1) */
    /* 	fprintf(stdout, "%s. (In fact, this reaction structure satisfies condition IC1 with *any stoichiometries*.)\n", IC1str); */
    /*   else */
    /* 	fprintf(stdout, "%s\n", IC1str); */
    /* } */

    /* if(csdflag!=2 && ssdflag!=2){ // check mass action */
    /*   if(wsdflag==3) */
    /* 	fprintf(stdout, "%s.\n", MAIC2IC3); */
    /*   else if(wsdflag==2) */
    /* 	fprintf(stdout, "%s.\n", MAIC3); */
    /*   else if(wsdflag==1) */
    /* 	fprintf(stdout, "%s.\n", MAIC2); */
    /* } */

    /* if(csdflag==0 && ssdflag==0 && wsdflag==0) */
    /*   fprintf(stderr, "Without further information on reversibility, etc. no claims about multistationarity seem possible.\n"); */

    /* free_imatrix(imat1, 0, nlen-1, 0, mlen-1); */
    /* free_imatrix(imat1a, 0, nlen-1, 0, mlena-1); */
  }

  fprintf(stdout, "\n");
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

//
// Get the rank of an integer matrix
//

int matrank(int **A, int n, int m){
  matrix J = imattoexmat(A, n, m);
  return J.rank();
}

int submatrank(int **A, int *vec1, int n, int *vec2, int m){
  matrix J = isubmattoexmat(A, vec1, n, vec2, m);
  return J.rank();
}



// Get a reaction as a single line in Sauro format
// Print it as a CRN to file fout. 

int getintline(char *s, char *fout){
  int i, j, jj, k;
  char *tmp;
  int nR,nS;
  int dg1, dg2;//assume reac then substrate
  int **matl, **matr;
  FILE *fd;
  bool flg;
  fd = fopen(fout, "w");
  if(!fd){
    fprintf(stderr, "ERROR in getintline: \"%s\" could not be opened for writing.\n", fout);
    exit(0);
  }
  //create a left array and a right array
  i=0;k=0;j=0;
  while(s[k] && !isdigit((int) s[k])){k++;}
  while(s[k] && isdigit((int) s[k])){j++;k++;}
  tmp=strchop2(s, k-j, j);nR=atoi(tmp);free(tmp);//number of reactions
  j=0;
  while(s[k] && !isdigit((int) s[k])){k++;}
  while(s[k] && isdigit((int) s[k])){j++;k++;}
  tmp=strchop2(s, k-j, j);nS=atoi(tmp);free(tmp);//number of species
  j=0;

  if(!nR || !nS){
    fprintf(stderr, "ERROR in getintline: \"%s\" is not a valid CRN as either no reactions or no species.\n", s);
    exit(0);
  }

  //initialise matrices
  matl=imatrix(0, nS-1, 0, nR-1);
  matr=imatrix(0, nS-1, 0, nR-1);
  for(i=0;i<nS;i++){
    for(jj=0;jj<nR;jj++){
      matl[i][jj]=0;matr[i][jj]=0;
    }
  }

  //create left and right stoichiometric matrices
  while(s[k]){
    while(s[k] && !isdigit((int) s[k])){k++;}
    while(s[k] && isdigit((int) s[k])){j++;k++;}
    tmp=strchop2(s, k-j, j);dg1=atoi(tmp);free(tmp);
    j=0;

    while(s[k] && !isdigit((int) s[k])){k++;}
    while(s[k] && isdigit((int) s[k])){j++;k++;}
    tmp=strchop2(s, k-j, j);dg2=atoi(tmp);free(tmp);
    j=0;

    if(dg1<dg2)//R-to-S
      (matr[dg2-nR][dg1])++;
    else//S-to-R
      (matl[dg1-nR][dg2])++;
  }

  //output to file
  for(i=0;i<nR;i++){
    flg=0;
    for(j=0;j<nS;j++){
      if(matl[j][i]){
	if(!flg){
	  if(matl[j][i]!=1)
	    fprintf(fd, "%dA%d",matl[j][i],j);
	  else
	    fprintf(fd, "A%d",j);
	  flg=1;
	}
	else{
	  if(matl[j][i]!=1)
	    fprintf(fd, "+%dA%d",matl[j][i],j);
	  else
	    fprintf(fd, "+A%d",j);
	}
      }
    }
    fprintf(fd, " --> ");
    flg=0;
    for(j=0;j<nS;j++){
      if(matr[j][i]){
	if(!flg){
	  if(matr[j][i]!=1)
	    fprintf(fd, "%dA%d",matr[j][i],j);
	  else
	    fprintf(fd, "A%d",j);
	  flg=1;
	}
	else{
	  if(matr[j][i]!=1)
	    fprintf(fd, "+%dA%d",matr[j][i],j);
	  else
	    fprintf(fd, "+A%d",j);
	}
      }
    }
    fprintf(fd, "\n");
  }

  free_imatrix(matl, 0, nS-1, 0, nR-1);
  free_imatrix(matr, 0, nS-1, 0, nR-1);

  fclose(fd);
  return 0;
}

// Read a list of reactions in Sauro Format 
// (http://128.208.17.26/NetworkEnumeration)
// and output these as individual files

int SauroSplit(char *fname, char *prefix){
  FILE *fd;
  int i;
  char oneline[1000];
  char fout[30];
  fd = fopen(fname, "r");
  if(!fd){
    fprintf(stderr, "ERROR in getintline: \"%s\" could not be opened for reading.\n", fname);
    exit(0);
  }
  while(getline0(fd, oneline, 1000) > 0){//get each line
    if(!iscomline(oneline)){
      i++;sprintf(fout, "%s%4.4d",prefix,i);
      fprintf(stderr, "fout = %s\n", fout);
      getintline(oneline, fout);
    }
  }
  fclose(fd);
  return 0;
}

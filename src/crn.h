#include <stdio.h>
#include <math.h>
#include <stdlib.h> /* for srand() */
#include <string.h>
#include <time.h> /* for random seeding */
#include <ctype.h>
#include <ginac/ginac.h>
#include <iostream>
#include <fstream>
using namespace std;
using namespace GiNaC;


#define NR_END 1
#define FREE_ARG char*

#define max(A, B) ((A) > (B) ? (A) : (B))
#define min(A, B) ((A) < (B) ? (A) : (B))

int **imatrix(long nrl, long nrh, long ncl, long nch);
void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch);
int **submat1(int **imat, int n, int i1, int i2);
int det(int **imat, int n);
int detsubmat(int **imat, int n, int m, int *i1, int *i2, int dim);
char *readfileintostr(const char fname[]);
char *getlinefromstr(long *pos, char *block);
int getnumints(char *s);
char *getnthint(char *s, int n);
unsigned long factorial(int x);
int **allperms(int n);
int par(int k);
int isinlist(int i, int lst[], int tot);
int **submatgen(int **imat, int n, int m, int *i1, int *i2, int dim);
int **allcombsgen(int n, int n1);
int isSSD2(int **imat, int n, int m);
int isSSD1(int **imat, int n, int m, bool allm, int q);
int isSSD(int **mat, int n, int m, int q);
int isonlyspace(char *s);

int isWSD2(int **imat, int n, int m);
int isWSD1(int **imat, int n, int m, int q);
int isCSD(int **imat, int n, int m, int q);
int isWS(int **mat, int n);
void simppair(int **imat, int **imat1, int n, int m, int ***imat2, int ***imat3, int *n1, int *m1);
void WSpair(int **imat, int n, int m, int ***imat2, int ***imat3, int *n1, int *m1);
int nextperm(int **vec, int **veclr, int *par, int n);
int strmatisSNSSS(char *fname);
int strmatisSSD(char *fname, int q);
int strmatcheckS2minors(char *fname, char *fnameout);
int strmatisCSD(char *fname, int q);
int minorisWS(int **imat1, int **imat2, int n, int m, int *vec1, int *vec2, int k, unsigned long fk);
int minorhas0rc(int **imat, int n, int m, int *vec1, int *vec2, int k);
int **submatminus(int **mat, int n, int m);
int **doublemat(int **imat, int n, int m);
int doubleisWSD2(int **imat, int n, int m);
int doubleisWSD1(int **imat, int n, int m);
int doubleisWSD(int **imat, int n, int m, int q);
int isWSD(int **mat, int n, int m, int q);
int mats_compat(int **imat1, int **imat2, int n, int m, int q);
int **detsk(int **imat, int n, int m, int k);
int iscomline(char s[]);
long getpos(int *vec, int n, int r);
long getpos1(int *vec, int omit, int n, int r);
int minor1(int **imat, int *xcombs, int *ycombs, int n, int m, int k, int **dets);
int **detsk1(int **imat, int n, int m, int k, int **dets);
int ***allminors(int **imat, int n, int m);
void free_allminors(int ***t, int n, int m);
int **allperms1(int *vec, int n);
long comb(int n, int k);
void printmat(int **imat, int n, int m);
void printbmat(bool **imat, int n, int m);
void printsubmat(int **imat, int *vec1, int *vec2, int k1, int k2);
void printindsubmat(int **imat, int *vec1, int *vec2, int k1, int k2);
void printmaximaindsubmat(FILE *fd, int **imat, int *vec1, int *vec2, int k1, int k2);
void printmat1(char **cmat, int **imat, int n, int m);
void printexmat(matrix imat, int n, int m);
int **simpmat(int **mat, int n, int m, int *n1, int *m1);
int **readmatrixfromstr(char *str, int *nlen, int *mlen);
int readmatpairfromstr(char *str, int *nlen, int *mlen, int ***mat1, int ***mat2);
int mat_signpat_compat(int **imat1, int **imat2, int n, int m, int q);
int getreac(char *str, char ***leftchems, int **leftstoics, int *numleft, char ***rightchems, int **rightstoics, int *numright, int *rev);
int chemgts2(char *s, char ***v, char sep);
int freearraydat(char **array, int lim);
int addv1(int k, char *s, char ***t);
int isinarray(char *v[], int numv, char *s);
int getallreacs(char *str, int ***imat1, int ***imat2, int ***imat3, int ***imat4, char ***chems, int *n, int *m, int *cols3, int *allrev, int *allgood);
int analysereacs(const char fname[], int q, bool htmlswitch);
int fixedminorcompat(int **imat1, int **imat2, int n, int m, int *vec1, int *vec2, int k);
int allminorsigns(int **imat1, int **imat2, int n, int m, int q);
int genMAXMAreacs(char *fname, int **imat1, int **imat2, int n, int m);
int S2(int **S, int n, int m, int ***imat1, int ***imat2, int *n1, int *m1);
int arecompat(int **imat1, int **imat2, int n, int m, int q);
int S2a(int **S, int **V, int n, int m, int ***imat1, int ***imat2, int *n1, int *m1);
int S2b(int **S, int **V, int n, int m, int ***imat1, int **valsvec, int **indsvec, int *base, int *basetot, int *n1, int *m1);
int areequal(int *vec1, int *vec2, int n);
void nextnum(int *vec, int n, int base);
void firstcomb(int *vec, int n, int n1);
int nextcomb(int *vec, int n, int n1);
void order(int *vec, int n);
void printvec(int *ivec, int n);
void printvec1(int *ivec, int n);
int **getallinds(int *longv, int *longv1, int *base, int *basetot, long numtot, int *comb, int n, long *k);
int S2arecompat(int **S, int **V, int n, int m, int q);
void nextcombk(int *vec, int n, int n1, long k);
void nextnumk(int *vec, int n, int base, long k);
int S2isWSD(char *fname, int q, int maxonly);
int minorisSNSsing2(int **imat, int n, int m, int *vec1, int *vec2, int k);
int qualdetsubmat(int **imat, int n, int m, int *i1, int *i2, int dim);
void printdotindsubmat(FILE *fd, int **imat, int *vec1, int *vec2, int k1, int k2, long lab);
int cuttails(int **mat, int n, int m, int *rw, int *col, int **n1, int **m1);
int addtwointstr(int k, int *s, int *s1, int len1, int len2, int ***t);
int S2detstocheck(char *fname, int q, int maxonly);
ex **exmatrix(long nrl, long nrh, long ncl, long nch);
bool **bmatrix(long nrl, long nrh, long ncl, long nch);
int findbad(int ** lists, int numlists, int *rowset, int latestcol, int level, int dim, bool **bmat, int **imat, matrix mmat, int n, int m, long *totbad, long *totall, int **hist);
ex intpairtopoly(int *l1, int cf1, int n1, int *l2, int cf2, int n2);
void ifree_l(int **imat, long n);
int **cpmat(int **mat, int n, int m);
matrix imattoexmat(int **A, int n, int m);
int matrank(int **A, int n, int m);

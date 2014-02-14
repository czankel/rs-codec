/*
 * fec.c -- forward error correction based on Vandermonde matrices
 * 980624
 * (C) 1997-98 Luigi Rizzo (luigi@iet.unipi.it)
 *
 * Portions derived from code by Phil Karn (karn@ka9q.ampr.org),
 * Robert Morelos-Zaragoza (robert@spectra.eng.hawaii.edu) and Hari
 * Thirumoorthy (harit@spectra.eng.hawaii.edu), Aug 1995
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above
 *    copyright notice, this list of conditions and the following
 *    disclaimer in the documentation and/or other materials
 *    provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHORS ``AS IS'' AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHORS
 * BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
 * OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,
 * OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
 * TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
 * OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY
 * OF SUCH DAMAGE.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "rs-codec.h"

/*
 * The following parameter defines how many bits are used for
 * field elements. The code supports any value from 2 to 16
 * but fastest operation is achieved with 8 bit elements
 * This is the only parameter you may want to change.
 */

#define GF_BITS  8	/* code over GF(2**GF_BITS) - change to suit */
#define	GF_SIZE ((1 << GF_BITS) - 1)	/* powers of \alpha */

/*
 * Primitive polynomials - see Lin & Costello, Appendix A,
 * and  Lee & Messerschmitt, p. 453.
 */

static const char *allPp[] = {  /* GF_BITS	polynomial    */
  NULL,                   /*  0 no code             */
  NULL,                   /*  1 no code             */
  "111",                  /*  2 1+x+x^2             */
  "1101",                 /*  3 1+x+x^3             */
  "11001",                /*  4 1+x+x^4             */
  "101001",               /*  5 1+x^2+x^5           */
  "1100001",              /*  6 1+x+x^6             */
  "10010001",             /*  7 1 + x^3 + x^7       */
  "101110001",            /*  8 1+x^2+x^3+x^4+x^8   */
  "1000100001",           /*  9 1+x^4+x^9           */
  "10010000001",          /* 10 1+x^3+x^10          */
  "101000000001",         /* 11 1+x^2+x^11          */
  "1100101000001",        /* 12 1+x+x^4+x^6+x^12    */
  "11011000000001",       /* 13 1+x+x^3+x^4+x^13    */
  "110000100010001",      /* 14 1+x+x^6+x^10+x^14   */
  "1100000000000001",     /* 15 1+x+x^15            */
  "11010000000010001"     /* 16 1+x+x^3+x^12+x^16   */
};

/*
 * To speed up computations, we have tables for logarithm, exponent
 * and inverse of a number. If GF_BITS <= 8, we use a table for
 * multiplication as well (it takes 64K, no big deal even on a PDA,
 * especially because it can be pre-initialized an put into a ROM!),
 * otherwhise we use a table of logarithms.
 * In any case the macro gf_mul(x,y) takes care of multiplications.
 */

static unsigned char gf_exp[2*GF_SIZE]; /* index->poly form conversion table	*/
static int gf_log[GF_SIZE + 1];         /* Poly->index form conversion table	*/
static unsigned char inverse[GF_SIZE+1];/* inverse of field elem.		*/
/* inv[\alpha**i]=\alpha**(GF_SIZE-i-1)	*/

/*
 * modnn(x) computes x % GF_SIZE, where GF_SIZE is 2**GF_BITS - 1,
 * without a slow divide.
 */
static inline unsigned char modnn(int x)
{
  while (x >= GF_SIZE) {
    x -= GF_SIZE;
    x = (x >> GF_BITS) + (x & GF_SIZE);
  }
  return x;
}

#define SWAP(a,b,t) {t tmp; tmp=a; a=b; b=tmp;}

/*
 * gf_mul(x,y) multiplies two numbers. If GF_BITS<=8, it is much
 * faster to use a multiplication table.
 *
 * USE_GF_MULC, GF_MULC0(c) and GF_ADDMULC(x) can be used when multiplying
 * many numbers by the same constant. In this case the first
 * call sets the constant, and others perform the multiplications.
 * A value related to the multiplication is held in a local variable
 * declared with USE_GF_MULC . See usage in addmul1().
 */
static unsigned char gf_mul_table[GF_SIZE + 1][GF_SIZE + 1];

#define gf_mul(x,y) gf_mul_table[x][y]

#define USE_GF_MULC register unsigned char * __gf_mulc_
#define GF_MULC0(c) __gf_mulc_ = gf_mul_table[c]
#define GF_ADDMULC(dst, x) dst ^= __gf_mulc_[x]

static void init_mul_table()
{
  int i, j;
  for (i=0; i< GF_SIZE+1; i++)
    for (j=0; j< GF_SIZE+1; j++)
      gf_mul_table[i][j] = gf_exp[modnn(gf_log[i] + gf_log[j]) ] ;

  for (j=0; j< GF_SIZE+1; j++)
    gf_mul_table[0][j] = gf_mul_table[j][0] = 0;
}

/*
 * Generate GF(2**m) from the irreducible polynomial p(X) in p[0]..p[m]
 * Lookup tables:
 *     index->polynomial form		gf_exp[] contains j= \alpha^i;
 *     polynomial form -> index form	gf_log[ j = \alpha^i ] = i
 * \alpha=x is the primitive element of GF(2^m)
 *
 * For efficiency, gf_exp[] has size 2*GF_SIZE, so that a simple
 * multiplication of two numbers can be resolved without calling modnn
 */

  static void
generate_gf(void)
{
  int i;
  unsigned char mask;
  const char *Pp =  allPp[GF_BITS] ;

  mask = 1;	/* x ** 0 = 1 */
  gf_exp[GF_BITS] = 0; /* will be updated at the end of the 1st loop */
  /*
   * first, generate the (polynomial representation of) powers of \alpha,
   * which are stored in gf_exp[i] = \alpha ** i .
   * At the same time build gf_log[gf_exp[i]] = i .
   * The first GF_BITS powers are simply bits shifted to the left.
   */
  for (i = 0; i < GF_BITS; i++, mask <<= 1 ) {
    gf_exp[i] = mask;
    gf_log[gf_exp[i]] = i;
    /*
     * If Pp[i] == 1 then \alpha ** i occurs in poly-repr
     * gf_exp[GF_BITS] = \alpha ** GF_BITS
     */
    if ( Pp[i] == '1' )
      gf_exp[GF_BITS] ^= mask;
  }
  /*
   * now gf_exp[GF_BITS] = \alpha ** GF_BITS is complete, so can als
   * compute its inverse.
   */
  gf_log[gf_exp[GF_BITS]] = GF_BITS;
  /*
   * Poly-repr of \alpha ** (i+1) is given by poly-repr of
   * \alpha ** i shifted left one-bit and accounting for any
   * \alpha ** GF_BITS term that may occur when poly-repr of
   * \alpha ** i is shifted.
   */
  mask = 1 << (GF_BITS - 1 ) ;
  for (i = GF_BITS + 1; i < GF_SIZE; i++) {
    if (gf_exp[i - 1] >= mask)
      gf_exp[i] = gf_exp[GF_BITS] ^ ((gf_exp[i - 1] ^ mask) << 1);
    else
      gf_exp[i] = gf_exp[i - 1] << 1;
    gf_log[gf_exp[i]] = i;
  }
  /*
   * log(0) is not defined, so use a special value
   */
  gf_log[0] =	GF_SIZE ;
  /* set the extended gf_exp values for fast multiply */
  for (i = 0 ; i < GF_SIZE ; i++)
    gf_exp[i + GF_SIZE] = gf_exp[i] ;

  /*
   * again special cases. 0 has no inverse. This used to
   * be initialized to GF_SIZE, but it should make no difference
   * since noone is supposed to read from here.
   */
  inverse[0] = 0 ;
  inverse[1] = 1;
  for (i=2; i<=GF_SIZE; i++)
    inverse[i] = gf_exp[GF_SIZE-gf_log[i]];
}

/*
 * Various linear algebra operations that i use often.
 */

/*
 * addmul() computes dst[] = dst[] + c * src[]
 * This is used often, so better optimize it! Currently the loop is
 * unrolled 16 times, a good value for 486 and pentium-class machines.
 * The case c=0 is also optimized, whereas c=1 is not. These
 * calls are unfrequent in my typical apps so I did not bother.
 * 
 * Note that gcc on
 */

#define addmul(dst, src, c, size) \
  if (c != 0) addmul1(dst, src, c, size)

#define UNROLL 16 /* 1, 4, 8, 16 */
static void addmul1(unsigned char *dst1,
    unsigned char *src1,
    unsigned char c,
    int size)
{
  USE_GF_MULC ;
  register unsigned char *dst = dst1, *src = src1 ;
  unsigned char *lim = &dst[size - UNROLL + 1] ;

  GF_MULC0(c) ;

#if (UNROLL > 1) /* unrolling by 8/16 is quite effective on the pentium */
  for (; dst < lim ; dst += UNROLL, src += UNROLL ) {
    GF_ADDMULC( dst[0] , src[0] );
    GF_ADDMULC( dst[1] , src[1] );
    GF_ADDMULC( dst[2] , src[2] );
    GF_ADDMULC( dst[3] , src[3] );
#if (UNROLL > 4)
    GF_ADDMULC( dst[4] , src[4] );
    GF_ADDMULC( dst[5] , src[5] );
    GF_ADDMULC( dst[6] , src[6] );
    GF_ADDMULC( dst[7] , src[7] );
#endif
#if (UNROLL > 8)
    GF_ADDMULC( dst[8] , src[8] );
    GF_ADDMULC( dst[9] , src[9] );
    GF_ADDMULC( dst[10] , src[10] );
    GF_ADDMULC( dst[11] , src[11] );
    GF_ADDMULC( dst[12] , src[12] );
    GF_ADDMULC( dst[13] , src[13] );
    GF_ADDMULC( dst[14] , src[14] );
    GF_ADDMULC( dst[15] , src[15] );
#endif
  }
#endif
  lim += UNROLL - 1 ;
  for (; dst < lim; dst++, src++ )		/* final components */
    GF_ADDMULC( *dst , *src );
}


/*
 * computes C = AB where A is n*k, B is k*m, C is n*m
 */

static void matmul(unsigned char *a,
    unsigned char *b,
    unsigned char *c,
    int n, int k, int m)
{
  int row, col, i ;

  for (row = 0; row < n ; row++) {
    for (col = 0; col < m ; col++) {
      unsigned char *pa = &a[ row * k ];
      unsigned char *pb = &b[ col ];
      unsigned char acc = 0 ;
      for (i = 0; i < k ; i++, pa++, pb += m )
        acc ^= gf_mul( *pa, *pb ) ;
      c[ row * m + col ] = acc ;
    }
  }
}


/*
 * invert_mat() takes a matrix and produces its inverse
 * k is the size of the matrix.
 * (Gauss-Jordan, adapted from Numerical Recipes in C)
 * Return non-zero if singular.
 */

static int invert_mat(unsigned char *src, int k)
{
  unsigned char c, *p ;
  int irow, icol, row, col, i, ix ;

  int error = 1 ;

  // FIXME statically allocate theses??
  int *indxc = (int*)alloca(k*sizeof(int));
  int *indxr = (int*)alloca(k*sizeof(int));
  int *ipiv = (int*)alloca(k*sizeof(int));
  unsigned char *id_row = (unsigned char*)alloca(k);

  bzero(id_row, k*sizeof(unsigned char));
  /*
   * ipiv marks elements already used as pivots.
   */
  for (i = 0; i < k ; i++)
    ipiv[i] = 0 ;

  for (col = 0; col < k ; col++) {
    unsigned char *pivot_row ;
    /*
     * Zeroing column 'col', look for a non-zero element.
     * First try on the diagonal, if it fails, look elsewhere.
     */
    irow = icol = -1 ;
    if (ipiv[col] != 1 && src[col*k + col] != 0) {
      irow = col ;
      icol = col ;
      goto found_piv ;
    }

    for (row = 0 ; row < k ; row++) {
      if (ipiv[row] != 1) {
        for (ix = 0 ; ix < k ; ix++) {
          if (ipiv[ix] == 0) {
            if (src[row*k + ix] != 0) {
              irow = row ;
              icol = ix ;
              goto found_piv ;
            }
          } else if (ipiv[ix] > 1) {
            fprintf(stderr, "singular matrix\n");
            goto fail ; 
          }
        }
      }
    }
    if (icol == -1) {
      fprintf(stderr, "XXX pivot not found!\n");
      goto fail ;
    }

found_piv:
    ++(ipiv[icol]) ;
    /*
     * swap rows irow and icol, so afterwards the diagonal
     * element will be correct. Rarely done, not worth
     * optimizing.
     */
    if (irow != icol) {
      for (ix = 0 ; ix < k ; ix++ ) {
        SWAP( src[irow*k + ix], src[icol*k + ix], unsigned char) ;
      }
    }
    indxr[col] = irow ;
    indxc[col] = icol ;
    pivot_row = &src[icol*k] ;
    c = pivot_row[icol] ;
    if (c == 0) {
      fprintf(stderr, "singular matrix 2\n");
      goto fail ;
    }
    if (c != 1 ) { /* otherwhise this is a NOP */
      /*
       * this is done often , but optimizing is not so
       * fruitful, at least in the obvious ways (unrolling)
       */
      c = inverse[ c ] ;
      pivot_row[icol] = 1 ;
      for (ix = 0 ; ix < k ; ix++ )
        pivot_row[ix] = gf_mul(c, pivot_row[ix] );
    }
    /*
     * from all rows, remove multiples of the selected row
     * to zero the relevant entry (in fact, the entry is not zero
     * because we know it must be zero).
     * (Here, if we know that the pivot_row is the identity,
     * we can optimize the addmul).
     */
    id_row[icol] = 1;
    if (bcmp(pivot_row, id_row, k*sizeof(unsigned char)) != 0) {
      for (p = src, ix = 0 ; ix < k ; ix++, p += k ) {
        if (ix != icol) {
          c = p[icol] ;
          p[icol] = 0 ;
          addmul(p, pivot_row, c, k );
        }
      }
    }
    id_row[icol] = 0;
  } /* done all columns */
  for (col = k-1 ; col >= 0 ; col-- ) {
    if (indxr[col] <0 || indxr[col] >= k)
      fprintf(stderr, "AARGH, indxr[col] %d\n", indxr[col]);
    else if (indxc[col] <0 || indxc[col] >= k)
      fprintf(stderr, "AARGH, indxc[col] %d\n", indxc[col]);
    else
      if (indxr[col] != indxc[col] ) {
        for (row = 0 ; row < k ; row++ ) {
          SWAP(src[row*k + indxr[col]], src[row*k + indxc[col]], unsigned char) ;
        }
      }
  }
  error = 0 ;
fail:
  return error ;
}

/*
 * fast code for inverting a vandermonde matrix.
 * XXX NOTE: It assumes that the matrix
 * is not singular and _IS_ a vandermonde matrix. Only uses
 * the second column of the matrix, containing the p_i's.
 *
 * Algorithm borrowed from "Numerical recipes in C" -- sec.2.8, but
 * largely revised for my purposes.
 * p = coefficients of the matrix (p_i)
 * q = values of the polynomial (known)
 */

static int invert_vdm(unsigned char *src, int k)
{
  int i, j, row, col ;
  unsigned char *b, *c, *p;
  unsigned char t, xx ;

  if (k == 1) 	/* degenerate case, matrix must be p^0 = 1 */
    return 0 ;

  // FIXME: statically allocate these?
  /*
   * c holds the coefficient of P(x) = Prod (x - p_i), i=0..k-1
   * b holds the coefficient for the matrix inversion
   */
  c = (unsigned char*)alloca(k);
  b = (unsigned char*)alloca(k);
  p = (unsigned char*)alloca(k);

  for ( j=1, i = 0 ; i < k ; i++, j+=k ) {
    c[i] = 0 ;
    p[i] = src[j] ;    /* p[i] */
  }
  /*
   * construct coeffs. recursively. We know c[k] = 1 (implicit)
   * and start P_0 = x - p_0, then at each stage multiply by
   * x - p_i generating P_i = x P_{i-1} - p_i P_{i-1}
   * After k steps we are done.
   */
  c[k-1] = p[0] ;	/* really -p(0), but x = -x in GF(2^m) */
  for (i = 1 ; i < k ; i++ ) {
    unsigned char p_i = p[i] ; /* see above comment */
    for (j = k-1  - ( i - 1 ) ; j < k-1 ; j++ )
      c[j] ^= gf_mul( p_i, c[j+1] ) ;
    c[k-1] ^= p_i ;
  }

  for (row = 0 ; row < k ; row++ ) {
    /*
     * synthetic division etc.
     */
    xx = p[row] ;
    t = 1 ;
    b[k-1] = 1 ; /* this is in fact c[k] */
    for (i = k-2 ; i >= 0 ; i-- ) {
      b[i] = c[i+1] ^ gf_mul(xx, b[i+1]) ;
      t = gf_mul(xx, t) ^ b[i] ;
    }
    for (col = 0 ; col < k ; col++ )
      src[col*k + row] = gf_mul(inverse[t], b[col] );
  }
  return 0 ;
}

/*
 * shuffle move src packets in their position
 */

static int shuffle(unsigned char *pkt[], int index[], int k)
{
  int i;

  for ( i = 0 ; i < k ; ) {

    if (index[i] >= k || index[i] == i)
      i++;
    else {
      /*
       * put pkt in the right position (first check for conflicts).
       */
      int c = index[i] ;

      if (index[c] == c)
        return 1;

      SWAP(index[i], index[c], int) ;
      SWAP(pkt[i], pkt[c], unsigned char *) ;
    }
  }
  return 0 ;
}

/*
 * build_decode_matrix constructs the encoding matrix given the
 * indexes. The matrix must be already allocated as
 * a vector of k*k elements, in row-major order
 */
static unsigned char * build_decode_matrix(struct rs_codec_parms *code,
    unsigned char *pkt[],
    int index[])
{
  int i , k = code->k ;
  unsigned char *matrix = code->dec_matrix;
  unsigned char *p;

  for (i = 0, p = matrix ; i < k ; i++, p += k ) {
    if (index[i] < k) {
      memset(p, 0, k*sizeof(unsigned char));
      p[i] = 1 ;
    } else {
      memcpy(p, &(code->enc_matrix[index[i]*k]), k*sizeof(unsigned char) ); 
    }
  }
  if (invert_mat(matrix, k)) {
  }
  return matrix ;
}



// ---------------------------------------------------------------------------

static int rs_codec_initialized = 0 ;

static void rs_codec_init()
{
  generate_gf();
  init_mul_table();
  rs_codec_initialized = 1;
}

/*
 * create a new encoder, returning a descriptor. This contains k,n and
 * the encoding matrix.
 */
struct rs_codec_parms * rs_codec_create(int k, int n)
{
  int row, col ;
  unsigned char *p, *tmp_m ;

  struct rs_codec_parms *codec ;

  if (rs_codec_initialized == 0)
    rs_codec_init();

  if (k > GF_SIZE + 1 || n > GF_SIZE + 1 || k > n ) {
    fprintf(stderr, "Invalid parameters k %d n %d GF_SIZE %d\n",
        k, n, GF_SIZE );
    return NULL ;
  }

  /* Allocate encoder and decoder matrix */
  codec = new rs_codec_parms();
  codec->k = k ;
  codec->n = n ;
  codec->enc_matrix = new unsigned char[n * k];
  codec->dec_matrix = new unsigned char[k * k];

  /*
   * fill the matrix with powers of field elements, starting from 0.
   * The first row is special, cannot be computed with exp. table.
   */

  tmp_m = new unsigned char[256 * 256];

  if (codec->enc_matrix == 0 || codec->dec_matrix == 0 || tmp_m == 0)
    goto errout;

  tmp_m[0] = 1 ;
  for (col = 1; col < k ; col++)
    tmp_m[col] = 0 ;
  for (p = tmp_m + k, row = 0; row < n-1 ; row++, p += k) {
    for ( col = 0 ; col < k ; col ++ )
      p[col] = gf_exp[modnn(row*col)];
  }

  /*
   * quick code to build systematic matrix: invert the top
   * k*k vandermonde matrix, multiply right the bottom n-k rows
   * by the inverse, and construct the identity matrix at the top.
   */
  invert_vdm(tmp_m, k); /* much faster than invert_mat */
  matmul(tmp_m + k*k, tmp_m, codec->enc_matrix + k*k, n - k, k, k);
  /*
   * the upper matrix is I so do not bother with a slow multiply
   */
  memset(codec->enc_matrix, 0, k*k*sizeof(unsigned char) );
  for (p = codec->enc_matrix, col = 0 ; col < k ; col++, p += k+1 )
    *p = 1 ;

  if (tmp_m)
    delete[] tmp_m;

  return codec;

errout:
  if (tmp_m != 0)
    delete[] tmp_m;
  if (codec != 0)
  {
    if (codec->enc_matrix != 0)
      delete[] codec->enc_matrix;
    if (codec->dec_matrix != 0)
      delete[] codec->dec_matrix;
  }
  return 0;
}


/*
 * This section contains the proper FEC encoding/decoding routines.
 * The encoding matrix is computed starting with a Vandermonde matrix,
 * and then transforming it into a systematic matrix.
 */

void rs_codec_free(struct rs_codec_parms *codec)
{
  if (codec != 0)
  {
    if (codec->enc_matrix != 0)
      delete[] codec->enc_matrix;
    if (codec->dec_matrix != 0)
      delete[] codec->dec_matrix;
    delete codec;
  }
}


/*
 * fec_encode accepts as input pointers to n data packets of size,
 * and produces as output a packet pointed to by fec, computed
 * with index "index".
 */
void rs_codec_encode(struct rs_codec_parms *codec,
    unsigned char *src[],
    unsigned char *fec,
    int index,
    int size)
{
  int i, k = codec->k ;
  unsigned char *p ;

  if (GF_BITS > 8)
    size /= 2 ;

  if (index < k || index >= codec->n) {
    fprintf(stderr, "Invalid index %d (max %d)\n",
        index, codec->n - 1 );
    return;
  }

  p = &(codec->enc_matrix[index*k] );
  memset(fec, 0, size*sizeof(unsigned char));
  for (i = 0; i < k ; i++)
    addmul(fec, src[i], p[i], size ) ;
}

/*
 * rs_decode receives as input a vector of packets, the indexes of
 * packets, and produces the correct vector as output.
 *
 * Input:
 *	code: pointer to code descriptor
 *	pkt:  pointers to received packets and . They are modified
 *	      to store the output packets (in place)
 *	index: pointer to packet indices ()
 *	size:    size of each packet
 */
int rs_codec_decode(struct rs_codec_parms *code,
    unsigned char *pkt[],
    int index[],
    int size)
{
  unsigned char *m_dec ; 
  int row, col, k = code->k ;

  if (GF_BITS > 8)
    size /= 2 ;

  if (shuffle(pkt, index, k))	/* error if true */
    return 1 ;

  m_dec = build_decode_matrix(code, pkt, index);

  if (m_dec == NULL)
    return 1 ; /* error */

  /*
   * do the actual decoding
   * use packets pkt[k...] for reconstructed packets
   */
  int p = k;
  for (row = 0 ; row < k ; row++ ) {
    if (index[row] >= k) {
      unsigned char *new_pkt = pkt[p++];
      memset(new_pkt, 0, size * sizeof(unsigned char) ) ;
      for (col = 0 ; col < k ; col++ )
        addmul(new_pkt, pkt[col], m_dec[row*k + col], size) ;
    }
  }

  /*
   * move pkts to their final destination
   */
  p = k;
  for (row = 0 ; row < k ; row++ ) {
    if (index[row] >= k) {
      pkt[row] = pkt[p++];
    }
  }

  return 0;
}

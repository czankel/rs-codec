/*
 * fec.c -- forward error correction based on Vandermonde matrices
 * 980614
 * (C) 1997-98 Luigi Rizzo (luigi@iet.unipi.it)
 *
 * Portions derived from code by Phil Karn (karn@ka9q.ampr.org),
 * Robert Morelos-Zaragoza (robert@spectra.eng.hawaii.edu) and Hari
 * Thirumoorthy (harit@spectra.eng.hawaii.edu), Aug 1995
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:

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


#ifndef RS_CODEC_H
#define RS_CODEC_H

struct rs_codec_parms {
    int k;
    int n;
    unsigned char *enc_matrix ;
    unsigned char *dec_matrix ;
};

/* Create an 'rs-codec' */
extern struct rs_codec_parms* rs_codec_create(int n, int k);

/* Free the specified 'rs-codec' */
extern void rs_codec_free(struct rs_codec_parms* codec);

/*
 * rs_codec_encode accepts as input pointers to n data packets of size,
 * and produces as output a packet pointed to by fec, computed
 * with index "index".
 */
extern void rs_codec_encode(struct rs_codec_parms *code,
                            unsigned char *src[],
                            unsigned char *fec,
                            int index,
                            int size);

extern void rs_codec_encode_one(struct rs_codec_parms* codec,
    unsigned char* dst[],
    unsigned char* src,
    unsigned int index,
    unsigned int offset,
    unsigned int size);

/*
 * rs_decode receives as input a vector of packets, the indexes of
 * packets, and produces the correct vector as output.
 *
 * Input:
 *	code: pointer to code descriptor
 *	pkt:  pointers to received packets.
 *	index: pointer to packet indexes
 *	size:  pointer to sizes of each packet
 *	max_size: max size of all packets
 */
extern int rs_codec_decode(struct rs_codec_parms* code,
                           unsigned char* pkt[],
                           int index[],
                           int size[],
                           int max_size);

#endif /* RS_CODEC_H */

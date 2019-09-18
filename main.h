#ifndef HEADER_MAIN_H_
#define HEADER_MAIN_H_

#include <stdlib.h>
#include <stdio.h>
//#include <string.h>
#include <stdint.h>
#include <math.h>
//#include <unistd.h>
//#include <stddef.h>
#include <time.h>
//#include <openssl/evp.h>
/////////////////////////////////////////
#define bool int
#define true 1
#define false 0

#define q 4294967295//2^32-1
#define p 65535//2^16-1
//#define qm 32
//#define pm 15
#define qp 65537
#define nn 1024
#define B 1

/*
#define q 1073741824 
#define p 536870912
#define qm 30
#define pm 29
#define qp 2
#define nn 1024
#define B 1
*/
/////////////////////////////////////////
#define modadd(c,a,b) \
do { \
  uint32_t _t = a+b; \
  c = _t + (_t < a); \
} while (0)

#define modsub(c,a,b) c = (a-b) - (b > a)

#define modmul(c,a,b) \
do { \
  uint64_t _T = (uint64_t) a * (uint64_t) b; \
  modadd (c, ((uint32_t) _T), ((uint32_t) ((uint64_t) _T >> (uint64_t) 32))); \
} while (0)


#define modmuladd(c,a,b) \
do { \
  uint64_t _T = (uint64_t) a * (uint64_t) b + c; \
  modadd (c, ((uint32_t) _T), ((uint32_t) ((uint64_t) _T >> (uint64_t) 32))); \
} while (0)

#define div2(c,a) c= (uint32_t) (((uint64_t) (a) + (uint64_t) ((uint32_t)(0-((a)&1))&0xFFFFFFFF))>>1)
#define normalize(c,a) c = (a) + ((a) == 0xFFFFFFFF)
#define DATATYPE uint32_t

#define SET_ZERO(x) (x)=0
#define add(c,a,b) modadd(c,a,b)
#define sub(c,a,b) modsub(c,a,b)
#define mul(c,a,b) modmul(c,a,b)
#define moddiv2(c,a)  normalize(c,a); div2(c,c)
#define neg(c,a)   (c)=0xFFFFFFFF-(a); normalize(c,c)
#define squ(c,a)   mul(c,a,a)
#define set(c,a)   (c)=(a)
#define NTESTS 1000

struct fft_ctx {
	uint32_t **x1;
	uint32_t **y1;
	uint32_t **z1;
	uint32_t *t1;
};
typedef struct fft_ctx FFT_CTX;

int n=1024,m=1024;
long long cpucycles(void);
//static void print_results(const char *s, unsigned long long *t, size_t tlen);
void print_results(const char *s, unsigned long long *t, size_t tlen);
void read(uint32_t a[nn], uint32_t b[nn]);
void readsharedb(uint32_t sp[nn],uint32_t R[nn]);
void getqpbR(uint32_t bR[nn],uint32_t r[nn],uint32_t R[nn]);
void reducefun(uint32_t bR[nn],uint64_t wp[nn]);
//void getwpp(uint32_t wpp[nn],uint64_t wp[nn]);
//void getkmp(uint32_t kmp[nn],uint64_t wp[nn]);
void getwppandkmp(uint32_t wpp[nn],uint32_t kmp[nn],uint64_t wp[nn]);
void recc(uint32_t km[nn],uint32_t w[nn],uint32_t wpp[nn]);
int FFT_CTX_init(FFT_CTX *ctx);
void FFT_mul(uint32_t *z, const uint32_t *x, const uint32_t *y, FFT_CTX *ctx);
void SameOrNo(uint32_t km[nn],uint32_t kmp[nn]);
#endif

//n=m=1024
//ÏµÊý2^16-1
#include"main.h"
//record the running time
unsigned long long t[NTESTS];

int main()
{

uint32_t  b[nn],bRsp[nn],asp[nn],bps[nn],sp[nn],rp[nn],a[nn],s[nn],bR[nn];//r[nn],w[nn],bp[nn];
uint32_t   wpp[nn],kmp[nn],km[nn];
uint64_t   wp[nn];

int i;

FFT_CTX ctx;
	if (!FFT_CTX_init(&ctx)) {//
		printf("Memory allocation error.");
		return -1;
	}


/*
int iterations = 0;
time_t starttime = time(NULL);
while (1) {
iterations++;
*/


  
   
        srand(time(NULL));

	//keygen
        read(a,s);//produce a and s at randoM 
         
        FFT_mul(b, a, s, &ctx);////multiply a and s

        for(i=0;i<n;i++){
           b[i]=b[i]/qp*qp;//qp=65537;*q/p-xiaquzhen-*q/p
	}




/*
if ((iterations % 100) == 0) {
				printf("Iterations: %d,  elapsed time: %ld\n",
				       iterations,  time(NULL) - starttime);
				if (iterations > (1 << 20)) break;
			    }
}
*/




/*
int iterations = 0;
time_t starttime = time(NULL);
while (1) {
iterations++;
*/

  
	//sharedb
	readsharedb(sp,rp);//produce sp and rp at random
	getqpbR(bR,b,rp);//bR=r+rp
	FFT_mul(bRsp, bR, sp, &ctx);////multiply bR and sp
	reducefun(bRsp,wp);

	//getwpp(wpp,wp);//get <w'> by w'
	//getkmp(kmp,wp);//get km by w'
	getwppandkmp(wpp,kmp,wp);
    
	FFT_mul(asp, a, sp, &ctx);////b'=asp

        for(i=0;i<n;i++){
           asp[i]=asp[i]/qp*qp;//b'=asp
	}
    

/*
if ((iterations % 100) == 0) {
				printf("Iterations: %d,  elapsed time: %ld\n",
				       iterations,  time(NULL) - starttime);
				if (iterations > (1 << 20)) break;
			    }
}
*/
	
    

        





int iterations = 0;
time_t starttime = time(NULL);
while (1) {
iterations++;


	
        FFT_mul(bps, asp, s, &ctx);////multiply b' and s, get w=bps

	recc(km,bps,wpp);//rec(w,<w'>)=km


if ((iterations % 100) == 0) {
				printf("Iterations: %d,  elapsed time: %ld\n",
				       iterations,  time(NULL) - starttime);
				if (iterations > (1 << 20)) break;
			    }
}



/*
printf("km\n");
for (i=0; i<nn; i++)
	{  		
            printf("%u ",km[i]);		
	} 

printf("\n");
*/

     
        SameOrNo(km,kmp);



return 0;


}

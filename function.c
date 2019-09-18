
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
//#include <unistd.h>
//#include <stddef.h>
#include <time.h>
//#include"main.h"
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


#define setbit(a,x) ((a)[(x)/64] |= (((uint64_t) 1) << (uint64_t) ((x)%64)))
#define getbit(a,x) (((a)[(x)/64] >> (uint64_t) ((x)%64)) & 1)

int np=1024,mp=1024;

void read(uint32_t a[nn], uint32_t b[nn]){
   
    int i,temp;

    for(i=0;i<np;i++){
        temp=rand()+rand();
	a[i]=temp;
	}
	
    for(i=0;i<mp;i++){
        temp=rand()&1;
	b[i]=temp;	
	}
	
}

void readsharedb(uint32_t sp[nn],uint32_t R[nn]){
    int i,temp;
	
    for(i=0;i<np;i++){
        temp=rand()&1;
	sp[i]=temp;       
	}
    for(i=0;i<mp;i++){
        temp=rand()%qp;
	R[i]=temp;
	}
}

void getqpbR(uint32_t bR[nn],uint32_t r[nn],uint32_t R[nn])
{
	int i;
        uint64_t temp,rr,RR;
	for(i=0;i<np;i++){
	rr=r[i];
	RR=R[i];
	temp=rr+RR;
        bR[i]=temp%q;
	}

}

void reducefun(uint32_t bR[nn],uint64_t wp[nn])
{

	int i;


	for(i=0;i<np;i++)
	{
                wp[i]=bR[i];
		
	}
}

void getwppandkmp(uint32_t wpp[nn],uint32_t kmp[nn],uint64_t wp[nn])
{
	int i;
 	//uint64_t wp;
	for(i=0;i<np;i++){
		//wp=brsp[i];
		wpp[i]=((wp[i]*15)/q)%5;
		kmp[i]=((wp[i]*3)/q)%3;
	}

}
/*
void getwpp(uint32_t wpp[nn],uint64_t wp[nn])
{
	int i;
	for(i=0;i<np;i++){
		wpp[i]=((wp[i]*15)/q)%5;
		//wpp[i]=wpp[i]&1;
	}
}
void getkmp(uint32_t kmp[nn],uint64_t wp[nn])
{
	int i;
	for(i=0;i<np;i++){	
		kmp[i]=((wp[i]*3)/q)%3;
	}

}
*/
/*
void recc(uint32_t km[nn],uint32_t w[nn],uint32_t wpp[nn])
{
   int32_t k;
   uint32_t i;

   int64_t ww,ii,jj,j;

   uint32_t temp;
   uint32_t deadline=8388608;

   for(k=0;k<np;k++)
   {
       for(i=0;i<deadline;i++)
       {

      		ww=w[k];
		ii=i;

		jj=ww+ii;

                if(jj>=q){j=jj%q;}
		else{j=jj;}
       	 	
		
                temp=(j*15/q)%5; 
                     
       		if(temp==wpp[k])
      		{         
          		 km[k]=((j*3)/q)%3; 
           		 break;      
       		}

       		jj=ww-ii;
       		if(jj<0){j=jj+q;}
		else{j=jj;}

		temp=(j*15/q)%5;
       		if(temp==wpp[k])
      		{     
           		km[k]=((j*3)/q)%3;
          		break;       
       		}
       } 
       if (i==deadline){km[k]=3;}  
   }
}
*/
/*
void recc(uint32_t km[nn],uint32_t w[nn],uint32_t wpp[nn])
{
   int32_t i,k,t;
   uint64_t ww;
   uint32_t wf,wj;
   uint32_t f[15]={0,0,0,0,0,1,1,1,1,1,2,2,2,2,2};
   uint32_t j[15]={0,1,2,3,4,0,1,2,3,4,0,1,2,3,4};
   for(i=0;i<np;i++)
   {
	ww=w[i];
        wf=(ww*3/q)%3;
	wj=(ww*15/q)%5;
	for(k=0;k<15;k++)
	{
	  if((f[k]==wf)&&(j[k]==wj))
	  {
		if(j[k]==wpp[i]){km[i]=f[k];break;}

                if(k+1>=15){
                t=k+1-15;
		if(j[t]==wpp[i]){
		km[i]=f[t];
		break;}}
		if(j[k+1]==wpp[i]){km[i]=f[k+1];break;}

		if(k-1<0){
                t=k-1+15;
		if(j[t]==wpp[i]){km[i]=f[t];
		break;}}
		if(j[k-1]==wpp[i]){km[i]=f[k-1];break;}

		if(k+2>=15){
                t=k+2-15;
		if(j[t]==wpp[i]){km[i]=f[t];
		break;}}
		if(j[k+2]==wpp[i]){km[i]=f[k+2];break;}

		if(k-2<0){
                t=k-2+15;
		if(j[t]==wpp[i]){km[i]=f[t];
		break;}}
		if(j[k-2]==wpp[i]){km[i]=f[k-2];break;}
                
          }
	}
    }
        

}
*/

void recc(uint32_t km[nn],uint32_t w[nn],uint32_t wpp[nn])
{
   int32_t i,k,t;
   uint64_t ww;
   
   uint32_t f[15]={0,0,0,0,0,1,1,1,1,1,2,2,2,2,2};
   uint32_t j[15]={0,1,2,3,4,0,1,2,3,4,0,1,2,3,4};
   for(i=0;i<np;i++)
   {
	ww=w[i];
        
	
	k=ww*15/q;
      
	  
		if(j[k]==wpp[i]){km[i]=f[k];continue;}

                if(k+1>=15){
                t=k+1-15;
		if(j[t]==wpp[i]){
		km[i]=f[t];
		continue;}}
		if(j[k+1]==wpp[i]){km[i]=f[k+1];continue;}

		if(k-1<0){
                t=k-1+15;
		if(j[t]==wpp[i]){km[i]=f[t];
		continue;}}
		if(j[k-1]==wpp[i]){km[i]=f[k-1];continue;}

		if(k+2>=15){
                t=k+2-15;
		if(j[t]==wpp[i]){km[i]=f[t];
		continue;}}
		if(j[k+2]==wpp[i]){km[i]=f[k+2];continue;}

		if(k-2<0){
                t=k-2+15;
		if(j[t]==wpp[i]){km[i]=f[t];
		continue;}}
		if(j[k-2]==wpp[i]){km[i]=f[k-2];continue;}
                
         
    }
}

void SameOrNo(uint32_t km[nn],uint32_t kmp[nn])
{
        int i;
        bool flag = true;

        for (i=0; i<np; i++)
        {            
                   if(km[i]!=kmp[i])
                   {
			printf("%d\n",i);
                        flag = false;

                        break;
                   }          
        }
    
        if(flag) 
        {
                printf("successful\n");
        }
        else
        {
                printf("fail\n");
        }

}


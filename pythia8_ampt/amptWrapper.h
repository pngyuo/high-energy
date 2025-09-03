/*
c++ wrapper converts ampt data and functions to accessible c++
data types and functions
*/

#ifndef AMPT_WRAPPER_H
#define AMPT_WRAPPER_H

#include <ctype.h>

//all fortran converted data common must be written in lower case!!!
//otherwise, the converted declarations can not be linked to
//their definitions in the fortran code

//-------------------------------------------------------------
//AMPT Common Block Declarations
const long maxstr=150001;
const long maxptn=400001;
const int maxidl=4001;

//index in fortran arrays starts from 1 while that in c starts 
//from 0!!!!!! katt, patt assignment will be affected
extern "C" {
	// COMMON /PARA1/ MUL
	extern struct 
    {
        int mul;    
    } para1_;
}
#define para1 para1_

extern "C" {
      // common /para2/ xmp, xmu, alpha, rscut2, cutof2
	extern struct 
    {
        double xmp, xmu, alpha, rscut2, cutof2; 
    } para2_;
}
#define para2 para2_

extern "C" {
    //  COMMON /INPUT2/ ILAB, MANYB, NTMAX, ICOLL, INSYS, IPOT, MODE, 
    // &   IMOMEN, NFREQ, ICFLOW, ICRHO, ICOU, KPOTEN, KMUL
    extern struct
    {
        int ilab, manyb, ntmax, icoll, insys, ipot, mode, imomen, nfreq, icflow, icrho, icou, kpoten, kmul;
    } input2_;
}
#define input2 input2_

extern "C" {
    // common/anim/nevent,isoft,isflag,izpc
    extern struct
    {
        int nevent, isoft, isflag, izpc;
    } anim_;
}
#define anim anim_

extern "C" {
    //  COMMON/HMAIN1/EATT,JATT,NATT,NT,NP,N0,N01,N10,N11
    extern struct
    {
        float eatt;
        int jatt,natt,nt,np,n0,n01,n10,n11;     
    } hmain1_;  
}
#define hmain1 hmain1_

extern "C" {
    //  COMMON/HMAIN2/KATT(MAXSTR,4),PATT(MAXSTR,4)
    extern struct
    {
        int katt[4][maxstr];
        float patt[4][maxstr];     
    } hmain2_;
}
#define hmain2 hmain2_

extern "C" {
    //  common /precpa/vxp0(MAXPTN),vyp0(MAXPTN),vzp0(MAXPTN),
    // 1       xstrg0(MAXPTN),ystrg0(MAXPTN),
    // 2       xstrg(MAXPTN),ystrg(MAXPTN),istrg0(MAXPTN),istrg(MAXPTN)
    extern struct
    {
        double vxp0[maxptn],vyp0[maxptn],vzp0[maxptn],xstrg0[maxptn],ystrg0[maxptn],
        xstrg[maxptn],ystrg[maxptn],istrg0[maxptn],istrg[maxptn];       
    } precpa_;  
}
#define precpa precpa_

extern "C" {
     //   COMMON /ARPRC/ ITYPAR(MAXSTR),
     // &     GXAR(MAXSTR), GYAR(MAXSTR), GZAR(MAXSTR), FTAR(MAXSTR),
     // &     PXAR(MAXSTR), PYAR(MAXSTR), PZAR(MAXSTR), PEAR(MAXSTR),
     // &     XMAR(MAXSTR)
    extern struct 
    {
        int itypar[maxstr];
        float gxar[maxstr],gyar[maxstr],gzar[maxstr],ftar[maxstr],pxar[maxstr],
        pyar[maxstr],pzar[maxstr],pear[maxstr],xmar[maxstr];        
    } arprc_;

}
#define arprc arprc_

extern "C" {
    //COMMON/ROOTDT/nline,id(MAXSTR),proot(4,MAXSTR), xroot(4,MAXSTR)
    extern struct 
    {
				int nline;
        int id[maxstr];
        float proot[4][maxstr],xroot[4][maxstr];        
    } rootdt_;
}
#define rootdt rootdt_


//-------------------------------------------------------------
//AMPT routines Declarations
#define run_ampt run_ampt_
#define ini_ampt ini_ampt_

extern "C" {
	void run_ampt(int*, int*, double*);
	void ini_ampt(float*, float*, int*);
}

inline void call_ini_ampt(float pzA, float pzB, int nseed) { ini_ampt(&pzA, &pzB, &nseed); }
inline void call_run_ampt(int event, int ncoll, double bimp) { run_ampt(&event, &ncoll, &bimp); }


#endif  // AMPT_WRAPPER_H

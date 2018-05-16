#ifndef PROTO_H_INCLUDED
#define PROTO_H_INCLUDED
#include <stdlib.h>
#include <math.h>
#include <stdio.h>


typedef struct {float t; float a; float f;} vect;
typedef struct {float xe; float ye; float ze; float xl; float yl; float zl;} pts;
typedef struct {float ae; float be; float ce; float ml; float nl; float ll; float xl; float yl; float zl;} data;

float proof (char mess) {
    float val;
    do {
        printf("%c = ", mess);
        if(!scanf("%f", &val)) {printf("Input error! Please try again.\n"); fflush(stdin);}
        else return val;
    } while(1);
}

void calccoord (vect p, data d, pts* res) {

    float cosa = cos(p.a); float sina = sin(p.a);
    float cosf = cos(p.f); float sinf = sin(p.f);

    (*res).xl = d.xl+p.t*d.ml;
    (*res).yl = d.yl+p.t*d.nl;
    (*res).zl = d.zl+p.t*d.ll;

    float rho = sqrt (d.ae*d.ae*d.be*d.be*d.ce*d.ce / (sinf*sinf*cosa*cosa*d.be*d.be*d.ce*d.ce + sinf*sinf*sina*sina*d.ae*d.ae*d.ce*d.ce + cosf*cosf*d.ae*d.ae*d.be*d.be));

    (*res).xe = rho*sinf*cosa;
    (*res).ye = rho*sinf*sina;
    (*res).ze = rho*cosf;

}

void calcgrad (vect p, data d, pts pt, vect* res) {

    float cosa = cos(p.a); float sina = sin(p.a);
    float cosf = cos(p.f); float sinf = sin(p.f);

    (*res).t = (-2)*( (pt.xe - pt.xl)*d.ml + (pt.ye - pt.yl)*d.nl + (pt.ze - pt.zl)*d.ll);

    float ga = 2* ( sinf*sinf*d.be*d.be*d.ce*d.ce*cosa*(-sina) + sinf*sinf*d.ae*d.ae*d.ce*d.ce*sina*cosa );
    float gf = 2* ( cosa*cosa*d.be*d.be*d.be*d.ce*d.ce*sinf + sina*sina*d.ae*d.ae*d.ce*d.ce*sina*cosf + d.ae*d.ae*d.be*d.be*cosf*(-sinf) );
    float abc = sqrt(d.ae*d.ae*d.be*d.be*d.ce*d.ce);
    float rt = sqrt(sinf*sinf*cosa*cosa*d.be*d.be*d.ce*d.ce + sinf*sinf*sina*sina*d.ae*d.ae*d.ce*d.ce + cosf*cosf*d.ae*d.ae*d.be*d.be);

    (*res).f = 2* ( (pt.xe - pt.xl)*abc*((cosf*cosa*rt)-(sinf*cosa*gf/(2*rt)))/(rt*rt) + (pt.ye - pt.yl)*abc*((cosf*sina*rt)-(sinf*sina*gf/(2*rt)))/(rt*rt) + (pt.ze - pt.zl)*abc*((-sinf*rt)-(cosf*gf/(2*rt)))/(rt*rt));
    (*res).a = 2* ( (pt.xe - pt.xl)*abc*((-sinf*sina*rt)-(sinf*cosa*ga/(2*rt)))/(rt*rt) + (pt.ye - pt.yl)*abc*((sinf*cosa*rt)-(sinf*sina*ga/(2*rt)))/(rt*rt) + (pt.ze - pt.zl)*abc*(-(cosf*ga/(2*rt)))/(rt*rt));

}

float calcdist (pts p) {
    float r = sqrt((p.xe-p.xl)*(p.xe-p.xl)+(p.ye-p.yl)*(p.ye-p.yl)+(p.ze-p.zl)*(p.ze-p.zl));
    return r;
}

#endif // PROTO_H_INCLUDED really hope this works

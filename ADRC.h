#pragma once
#ifndef LADRC_H
#define LADRC_H
#include "Fpu87.h"
class ADRC {
	//TD parameters
	float r;
	float h;
	float h0;
	float td_x1;
	float td_x2;
	//ESO parameters
	float delta;
	float b;
	float beta01;
	float beta02;
	float beta03;
	float z1 ;
	float z2 ;
	float z3 ;
	//NLSEF parameters
	float alpha1;
	float alpha2;
	float beta1;
	float beta2;
	float last_u;
	// upper & lower_limit
	float out_upper_limit;
	float out_lower_limit;
public:
	 ADRC();
	void setTDpara(const float &r, const float &h, const float &h0);
	void setESOpara(const float &delta, const float &b, const float &beta01, const float &beta02, const float &beta03);
	void setNLSEApara(const float &alpha1, const float &alpha2, const float &beta1, const float &beta2);
	void setupperlimit(const float &out_upper_limit);
	void setlowerlimit(const float &out_lower_limit);
	float fst(float x1, float x2, float v);
	float fal(float e, float alfa, float delta);
	float sign(float x);
	float refresh(float v, float y);
	};

#endif // LADRC_H

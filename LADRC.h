#pragma once
#ifndef LADRC_H
#define LADRC_H
#include "Fpu87.h"
class LADRC {
	//TD parameters
	float r;
	float h;
	float h0;
	float td_x1;
	float td_x2;
	//ESO parameters
	float omega_o;
	float b;
	float beta01;
	float beta02;
	float beta03;
	float z1 ;
	float z2 ;
	float z3 ;
	//PD parameters
	float omega_c;
	float Kp;
	float Kd;
	float last_u;
	// upper & lower_limit
	float out_upper_limit;
	float out_lower_limit;
public:
	 LADRC();
	void setTDpara(const float &r, const float &h, const float &h0);
	void setESOpara(const float &b, const float &omega_o);
	void LADRC::setPDpara(const float &omega_c);
	void setupperlimit(const float &out_upper_limit);
	void setlowerlimit(const float &out_lower_limit);
	float fst(float x1, float x2, float v);
	float sign(float x);
	float refresh(float v, float y);
	};

#endif // LADRC_H

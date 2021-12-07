#include "LADRC.h"
	LADRC::LADRC() :
	//TD parameters
	td_x1(0.f),
	td_x2(0.f),
	z1(0.f),
	z2(0.f),
	z3(0.f),
	last_u(0.f)
	{}
	
	void LADRC::setTDpara(const float &r, const float &h, const float &h0)
	{
		this->r = r;
		this->h = h;
		this->h0 = h0;
	}
	void LADRC::setESOpara(const float &b, const float &omega_o)
	{
		
		this->b = b;
		this->beta01 = 3 * omega_o;
		this->beta02 = 3 * omega_o*omega_o;
		this->beta03 = omega_o*omega_o*omega_o;
	}
	void LADRC::setPDpara(const float &omega_c)
	{
		this->omega_c = omega_c;
		this->Kp=omega_c*omega;
		this->kd=2*omega-1;
	}
	void LADRC::setupperlimit(const float &out_upper_limit)
	{
		this->out_upper_limit = out_upper_limit;
	}
	void LADRC::setlowerlimit(const float &out_lower_limit)
	{
		this->out_lower_limit = out_lower_limit;
	}
	float LADRC::fst(float x1, float x2, float v)
	{
		float td_y = 0;
		float a0 = 0;
		float a = 0;
		float fhan = 0;
		float d = 0;
		float d0 = 0;
	

		d = r*h;
		d0 = h*d;
		td_y = x1 - v + h*x2;
		a0 = sqrt_(d*d + 8 * r*fabs_(td_y));


		if (fabs_(td_y)>d0)
			a = x2 + 0.5*(a0 - d)*sign(td_y);
		else
			a = x2 + td_y / h;


		if (fabs_(a)>d)
			fhan = -r*sign(a);
		else
			fhan = -r*a / d;
		return(fhan);
	}
	float LADRC::sign(float x)
	{
		if (x>0)
			return(1);
		if (x<0)
			return(-1);
		if (x == 0)
			return(0);
	}

	float LADRC::refresh(float v, float y)
	{
		float u0 = 0;
		float e = 0;
		float e1 = 0;
		float e2 = 0;

		//**********   TD  ************跟踪微分器
		td_x1 = td_x1 + h*td_x2;                        //td_x1=v1;
		td_x2 = td_x2 + h*fst(td_x1, td_x2, v);           //td_x2=v2;
														  //********  leso  *************观测器
		e = z1 - y;
		z1 = z1 + h*(z2 - beta01*e);
		z2 = z2 + h*(z3 - beta02*e + b*last_u);
		z3 = z3 - h*beta03*e;
		//***********  PD *************控制器
		e1 = td_x1 - z1;   				//e1=v1-z1;   
		e2 = td_x2 - z2;				//e2=v2-z2;
										//    e1=v-z1;
										//	e2=-z2;
		last_u = (Kp*e1 - Kd*z2)-z3/b;
		if (last_u<out_lower_limit) last_u = out_lower_limit;
		if (last_u>out_upper_limit)  last_u = out_upper_limit;
		return(last_u);                 //u=u0-z3/b;
	}

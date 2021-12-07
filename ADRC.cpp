#include "ADRC.h"
	ADRC::ADRC() :
	//TD parameters
	td_x1(0.f),
	td_x2(0.f),
	z1(0.f),
	z2(0.f),
	z3(0.f),
	last_u(0.f)
	{}
	
	void ADRC::setTDpara(const float &r, const float &h, const float &h0)
	{
		this->r = r;
		this->h = h;
		this->h0 = h0;
	}
	void ADRC::setESOpara(const float &delta, const float &b, const float &beta01, const float &beta02, const float &beta03)
	{
		this->delta = delta;
		this->b = b;
		this->beta01 = beta01;
		this->beta02 = beta02;
		this->beta03 = beta03;
	}
	void ADRC::setNLSEApara(const float &alpha1, const float &alpha2, const float &beta1, const float &beta2)
	{
		this->alpha1 = alpha1;
		this->alpha2 = alpha2;
		this->beta1 = beta1;
		this->beta2 = beta2;
	}
	void ADRC::setupperlimit(const float &out_upper_limit)
	{
		this->out_upper_limit = out_upper_limit;
	}
	void ADRC::setlowerlimit(const float &out_lower_limit)
	{
		this->out_lower_limit = out_lower_limit;
	}
	float ADRC::fst(float x1, float x2, float v)
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
	float ADRC::fal(float e, float alfa, float delta)
	{
		
		float y = 0.0;
		if (fabs_(e)>delta) y = pow_(fabs_(e), alfa)*sign(e);
		else			  y = e / pow_(delta, 1.0 - alfa);
		return(y);
	}
	float ADRC::sign(float x)
	{
		if (x>0)
			return(1);
		if (x<0)
			return(-1);
		if (x == 0)
			return(0);
	}

	float ADRC::refresh(float v, float y)
	{
		float u0;
		float e = 0;
		float e1 = 0;
		float e2 = 0;

		//**********   TD  ************¸ú×ÙÎ¢·ÖÆ÷
		td_x1 = td_x1 + h*td_x2;                        //td_x1=v1;
		td_x2 = td_x2 + h*fst(td_x1, td_x2, v);           //td_x2=v2;
														  //********  eso  *************¹Û²âÆ÷
		e = z1 - y;
		z1 = z1 + h*(z2 - beta01*e);
		z2 = z2 + h*(z3 - beta02*fal(e, 0.5, delta) + b*last_u);
		z3 = z3 - h*beta03*fal(e, 0.25, delta);
		//***********  NLSEF *************¿ØÖÆÆ÷
		e1 = td_x1 - z1;   				//e1=v1-z1;   
		e2 = td_x2 - z2;				//e2=v2-z2;
										//    e1=v-z1;
										//	e2=-z2;
		last_u = beta1*fal(e1, alpha1, delta) + beta2*fal(e2, alpha2, delta)-z3/b;// u = u0 - z3 / b;
		if (last_u<out_lower_limit) last_u = out_lower_limit;
		if (last_u>out_upper_limit)  last_u = out_upper_limit;
		return(last_u);                 //
	}
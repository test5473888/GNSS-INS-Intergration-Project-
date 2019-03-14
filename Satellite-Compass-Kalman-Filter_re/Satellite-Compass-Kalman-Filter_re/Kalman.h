#pragma once
#include<stdlib.h>

struct Kalman
{
	double	*X,	**P,**Q,				
		**F,*z,**R,**H,
		**K,**G,*u;
}KF;


enum KFmode
{ 
	Model_set,
	Pridict_mode,
	Update_mode
};



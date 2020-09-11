/*
	Author: Colton Bryant 
	Initially created 06/11/20

	For setting up an initial interface of the form h(x) = a|x|^3 + b|x|^2 + c|x| + d
*/

#include "mirroredcubic.h"

namespace levelset {
	void MirroredCubic::SetParams(InputParams *params) 
	{
		if (parameter) delete[] parameter;
		parameter = new double[4];
		parameter[A] = params->GetDoubleParam("x^3 coeff");
		parameter[B] = params->GetDoubleParam("x^2 coeff");
		parameter[C] = params->GetDoubleParam("x coeff");
		parameter[D] = params->GetDoubleParam("constant");
	}

	double MirroredCubic::X(const double t) const
	{
		return t;
	}

	double MirroredCubic::Y(const double t) const
	{
		return parameter[A]*fabs(t)*fabs(t)*fabs(t) + 
				parameter[B]*fabs(t)*fabs(t) + 
				parameter[C]*fabs(t) + 
				parameter[D];
	}

}





#pragma once
#include "stdafx.h"

class CarlsonFunction
{
	public:
		CarlsonFunction();
		~CarlsonFunction();
		double rcFunction(const double& x, const double& y);
		double rdFunction(double x, double y, double z);
		double rfFunction(double x, double y, double z);
		double rjFunction(const double& x, const double& y, const double& z, const double& p);

	private:
		void iterationForRCFunction(double& x, double& y, double& ave, double& sol);
		void iterationForRDFunction(double& x, double& y, double& z, double& dX, double& dY, double& dZ, double& ave, double& sol, double& f);
		void iterationForRFFunction(double& x, double& y, double& z, double& dX, double& dY, double& dZ, double& ave);
		void iterationForRJFunction(double& x, double& y, double& z, double& dX, double& dY, double& dZ, double& dP, double& ave, double& so, double& f, double& p);
};
#pragma once
#include "stdafx.h"

class CarlsonFunction
{
	public:
		CarlsonFunction();
		~CarlsonFunction();
		double_t rcFunction(const double_t& x, const double_t& y);
		double_t rdFunction(double_t x, double_t y, double_t z);
		double_t rfFunction(double_t x, double_t y, double_t z);
		double_t rjFunction(const double_t& x, const double_t& y, const double_t& z, const double_t& p);

	private:
		void iterationForRCFunction(double_t& x, double_t& y, double_t& ave, double_t& sol);
		void iterationForRDFunction(double_t& x, double_t& y, double_t& z, double_t& dX, double_t& dY, double_t& dZ, double_t& ave, double_t& sol, double_t& f);
		void iterationForRFFunction(double_t& x, double_t& y, double_t& z, double_t& dX, double_t& dY, double_t& dZ, double_t& ave);
		void iterationForRJFunction(double_t& x, double_t& y, double_t& z, double_t& dX, double_t& dY, double_t& dZ, double_t& dP, double_t& ave, double_t& so, double_t& f, double_t& p);
};
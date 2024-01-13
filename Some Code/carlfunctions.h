#pragma once
#include "stdafx.h"

class CarlsonFunction
{
	public:
		CarlsonFunction();
		~CarlsonFunction();
		double_t rcFunction(const double_t& x, const double_t& y);
		double_t rdFunction(double_t x, double_t y, double_t z);

	private:
		void iterationForRCFunction(double_t& x, double_t& y, double_t& ave, double_t& sol);
		void iterationForRDFunction(double_t& x, double_t& y, double_t& z, double_t& dX, double_t& dY, double_t& dZ, double_t& ave, double_t& sol, double_t& f);
};
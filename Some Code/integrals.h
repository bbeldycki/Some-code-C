#pragma once
#include "stdafx.h"

class Integrals
{
	public:
		Integrals();
		~Integrals();
		double_t ellIntCubicAllRootsReal(std::vector<int>& pList, std::vector<double_t>& aList, std::vector<double_t>& bList, double_t& ffr, double_t& y, double_t& x);

	private:
		std::unique_ptr<CarlsonFunction> cFunctions;
		std::unique_ptr<IntegralsHelper> intHelper;
};
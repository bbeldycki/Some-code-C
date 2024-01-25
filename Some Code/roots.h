#pragma once
#include "stdafx.h"

class Roots
{
	public:
		Roots();
		~Roots();
		void laguer(std::vector<std::complex<double>>& a, std::complex<double>& x, int& m);

	private:
		//void laguer(std::vector<std::complex<double>>& a, std::complex<double>& x, int& m);
};
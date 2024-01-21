#pragma once
#include "integralshelper.h"
#include "carlfunctions.h"

class Integrals
{
	public:
		Integrals();
		~Integrals();
		double ellIntCubicAllRootsReal(std::vector<int>& pList, std::vector<double>& aList, std::vector<double>& bList, double& ffr, double& y, double& x);
		double ellIntCubicOneRealTwoComplexRoots(std::vector<int>& pList, std::vector<double>& aList, std::vector<double>& bList, std::vector<double>& fghList, double& ffr, double& y, double& x);
		double ellIntQuarticAllRootsReal(std::vector<int>& pList, std::vector<double>& aList, std::vector<double>& bList, double& ffr, double& y, double& x);
		double ellIntQuarticTwoRealTwoComplexRoots(std::vector<int>& pList, std::vector<double>& aList, std::vector<double>& bList, std::vector<double>& fghList, double& ffr, double& y, double& x);

	private:
		std::unique_ptr<CarlsonFunction> cFunctions;
		std::unique_ptr<IntegralsHelper> intHelper;
};
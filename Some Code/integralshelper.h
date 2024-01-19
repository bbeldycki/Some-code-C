#pragma once
#include "stdafx.h"

class IntegralsHelper
{
	public:
		IntegralsHelper();
		~IntegralsHelper();
		double computeAPOnePn(std::vector<int>& pList, std::vector<double>& aList, std::vector<double>& bList, double& y, double& x);
		double computeAPOnePnDoubleQuadratic(std::vector<int>& pList, std::vector<double>& aList, std::vector<double>& bList, std::vector<double>& fghOne, std::vector<double>& fghTwo, double& y, double& x);
		double computeAPOnePnQuadratic(std::vector<int>& pList, std::vector<double>& aList, std::vector<double>& bList, std::vector<double>& fghList, double& y, double& x);
		double computeCoeffAlfai(std::vector<double>& fList, std::vector<double>& gList, double& a, double& b, int k);
		double computeCoeffBetai(std::vector<double>& gList, std::vector<double>& hList, double& a, double& b, int k);
		double computeCijValue(std::vector<double>& aList, std::vector<double>& bList, std::vector<double>& fghList, int k, int l);
		double computeDeltaij(std::vector<double>& fList, std::vector<double>& gList, std::vector<double>& hList, int k, int l);
		double ComputeGammai(std::vector<double>& fList, std::vector<double>& gList, std::vector<double>& hList, double& a, double& b, int k);

};
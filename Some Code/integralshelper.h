#pragma once
#include "stdafx.h"

class IntegralsHelper
{
	public:
		IntegralsHelper();
		~IntegralsHelper();
		double_t computeAPOnePn(std::vector<int>& pList, std::vector<double_t>& aList, std::vector<double_t>& bList, double_t& y, double_t& x);
		double_t computeAPOnePnDoubleQuadratic(std::vector<int>& pList, std::vector<double_t>& aList, std::vector<double_t>& bList, std::vector<double_t>& fghOne, std::vector<double_t>& fghTwo, double_t& y, double_t& x);
		double_t computeAPOnePnQuadratic(std::vector<int> pList, std::vector<double_t> aList, std::vector<double_t> bList, std::vector<double_t> fghList, double_t y, double_t x);
		double_t computeCoeffAlfai(std::vector<double_t> fList, std::vector<double_t> gList, double_t a, double_t b, int k);
		double_t computeCoeffBetai(std::vector<double_t> gList, std::vector<double_t> hList, double_t a, double_t b, int k);
		double_t computeCijValue(std::vector<double_t> aList, std::vector<double_t> bList, std::vector<double_t> fghList, int k, int l);
		double_t computeDeltaij(std::vector<double_t> fList, std::vector<double_t> gList, std::vector<double_t> hList, int k, int l);
		double_t ComputeGammai(std::vector<double_t> fList, std::vector<double_t> gList, std::vector<double_t> hList, double_t a, double_t b, int k);

};
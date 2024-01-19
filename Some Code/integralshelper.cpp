#include "integralshelper.h"

IntegralsHelper::IntegralsHelper()
{

}

IntegralsHelper::~IntegralsHelper()
{

}

double IntegralsHelper::computeAPOnePn(std::vector<int>& pList, std::vector<double>& aList, std::vector<double>& bList, double& y, double& x)
{
	std::vector<double> xi, yi;
	double handlerX, handlerY, ax, ay;
	try
	{
		if (x < y)
		{
			throw std::runtime_error("Error: x is lower than y. computeAPOnePn.");
		}
	}
	catch (const std::exception& e)
	{
		std::cerr << "Exception caught: " << e.what() << std::endl;
		return 0;
	}
	
	for (int i = 0; i < aList.size(); i++)
	{
		handlerX = aList[i] + bList[i] * x;
		handlerY = aList[i] + bList[i] * y;
		try 
		{
			if (handlerX < 0)
			{
				throw std::runtime_error("Error: handlerX is lower than 0. computeAPOnePn.");
			}
			if (handlerY < 0)
			{
				throw std::runtime_error("Error: handlerX is lower than 0. computeAPOnePn.");
			}
		}
		catch (const std::exception& e)
		{
			std::cerr << "Exception caught: " << e.what() << std::endl;
			return 0;
		}
		xi.push_back(sqrt(handlerX));
		yi.push_back(sqrt(handlerY));
	}
	ax = ay = 1;
	for (int i = 0; i < aList.size(); i++)
	{
		ax = ax * pow(xi[i], pList[i]);
		ay = ay * pow(yi[i], pList[i]);
	}
	return ax - ay;
}

double IntegralsHelper::computeAPOnePnDoubleQuadratic(std::vector<int>& pList, std::vector<double>& aList, std::vector<double>& bList, std::vector<double>& fghOne, std::vector<double>& fghTwo, double& y, double& x)
{
	std::vector<double> ksi, eta;
	std::vector<double> fi, gi, hi;
	double handlerKsi, handlerEta, ax, ay;
	double ksiFive = 0;
	double etaFive = 0;
	try
	{
		if (x < y)
		{
			throw std::runtime_error("Error: x is lower than y. computeAPOnePnDoubleQuadratic.");
		}
	}
	catch (const std::exception& e)
	{
		std::cerr << "Exception caught: " << e.what() << std::endl;
		return 0;
	}

	fi.push_back(fghOne[0]);
	fi.push_back(fghTwo[0]);
	gi.push_back(fghOne[1]);
	gi.push_back(fghTwo[1]);
	hi.push_back(fghOne[2]);
	hi.push_back(fghTwo[2]);

	for (int i = 0; i < fi.size(); i++)
	{
		handlerKsi = fi[i] + gi[i] * x + hi[i] * x * x;
		handlerEta = fi[i] + gi[i] * y + hi[i] * y * y;
		try
		{
			if (handlerKsi < 0)
			{
				throw std::runtime_error("Error: handlerKsi is lower than 0. computeAPOnePnDoubleQuadratic.");
			}
			if (handlerEta < 0)
			{
				throw std::runtime_error("Error: handlerEta is lower than 0. computeAPOnePnDoubleQuadratic.");
			}
		}
		catch (const std::exception& e)
		{
			std::cerr << "Exception caught: " << e.what() << std::endl;
			return 0;
		}
		ksi.push_back(sqrt(handlerKsi));
		eta.push_back(sqrt(handlerEta));
	}
	ax = ay = 1;
	if (pList[4] != 0)
	{
		ksiFive = aList[4] + bList[4] * x;
		etaFive = aList[4] + bList[4] * y;
		for (int i = 0; i < pList.size(); i++)
		{
			if (i < 2)
			{
				ax = ax * pow(ksi[i], pList[i]);
				ay = ay * pow(eta[i], pList[i]);
			}
			else if (i == 4)
			{
				ax = ax * sqrt(pow(ksiFive, pList[i]));
				ay = ay * sqrt(pow(etaFive, pList[i]));
			}
		}
	}
	else
	{
		for (int i = 0; i < ksi.size(); i++)
		{
			ax = ax * pow(ksi[i], pList[i]);
			ay = ay * pow(eta[i], pList[i]);
		}
	}
	return ax - ay;
}

double IntegralsHelper::computeAPOnePnQuadratic(std::vector<int>& pList, std::vector<double>& aList, std::vector<double>& bList, std::vector<double>& fghList, double& y, double& x)
{
	std::vector<double> xi, yi;
	std::vector<double> fi, gi, hi;
	double handlerX, handlerY, handlerKsi, handlerEta, ksi, eta, ax, ay;
	double ksiFive = 0;
	double etaFive = 0;
	try
	{
		if (x < y)
		{
			throw std::runtime_error("Error: x is lower than y. computeAPOnePnQuadratic.");
		}
	}
	catch (const std::exception& e)
	{
		std::cerr << "Exception caught: " << e.what() << std::endl;
		return 0;
	}

	for (int i = 0; i < aList.size(); i++)
	{
		handlerX = aList[i] + bList[i] * x;
		handlerY = aList[i] + bList[i] * y;
		try
		{
			if (handlerX < 0)
			{
				throw std::runtime_error("Error: handlerX is lower than 0. computeAPOnePn.");
			}
			if (handlerY < 0)
			{
				throw std::runtime_error("Error: handlerX is lower than 0. computeAPOnePn.");
			}
		}
		catch (const std::exception& e)
		{
			std::cerr << "Exception caught: " << e.what() << std::endl;
			return 0;
		}
		xi.push_back(sqrt(handlerX));
		yi.push_back(sqrt(handlerY));
	}
	handlerKsi = fghList[0] + fghList[1] * x + fghList[2] * x * x;
	handlerEta = fghList[0] + fghList[1] * y + fghList[2] * y * y;
	try
	{
		if (handlerKsi < 0)
		{
			throw std::runtime_error("Error: handlerKsi is lower than 0. computeAPOnePnDoubleQuadratic.");
		}
		if (handlerEta < 0)
		{
			throw std::runtime_error("Error: handlerEta is lower than 0. computeAPOnePnDoubleQuadratic.");
		}
	}
	catch (const std::exception& e)
	{
		std::cerr << "Exception caught: " << e.what() << std::endl;
		return 0;
	}
	ksi = sqrt(handlerKsi);
	eta = sqrt(handlerEta);
	ax = ay = 1;
	for (int i = 0; i < aList.size(); i++)
	{
		if (i != 1)
		{
			if (i != 2)
			{
				ax = ax * pow(xi[i], pList[i]);
				ay = ay * pow(yi[i], pList[i]);
			}
		}
		else
		{
			ax = ax * pow(ksi, pList[i]);
			ay = ay * pow(eta, pList[i]);
		}
	}
	return ax - ay;
}

double IntegralsHelper::computeCoeffAlfai(std::vector<double>& fList, std::vector<double>& gList, double& a, double& b, int k)
{
	return double(2.0 * fList[k] * b - gList[k] * a);
}

double IntegralsHelper::computeCoeffBetai(std::vector<double>& gList, std::vector<double>& hList, double& a, double& b, int k)
{
	return double(gList[k] * b - 2.0 * hList[k] * a);
}

double IntegralsHelper::computeCijValue(std::vector<double>& aList, std::vector<double>& bList, std::vector<double>& fghList, int k, int l)
{
	return double(2.0 * bList[k] * bList[l] - fghList[1] * (aList[k] * bList[l] + aList[l] * bList[k]) + 2.0 * fghList[2] * aList[k] * aList[l]);
}

double IntegralsHelper::computeDeltaij(std::vector<double>& fList, std::vector<double>& gList, std::vector<double>& hList, int k, int l)
{
	return double(2.0 * fList[k] * hList[l] + 2.0 + fList[l] * hList[k] - gList[k] * gList[l]);
}

double IntegralsHelper::ComputeGammai(std::vector<double>& fList, std::vector<double>& gList, std::vector<double>& hList, double& a, double& b, int k)
{
	double gamma = 0;
	gamma = fList[k] * pow(b, 2) - gList[k] * a * b + hList[k] * pow(a, 2);
	try
	{
		if (gamma <= 0)
		{
			throw std::runtime_error("Error: gamma is lower than 0 or equal to 0. ComputeGammai.");
		}
	}
	catch (const std::exception& e)
	{
		std::cerr << "Exception caught: " << e.what() << std::endl;
		return 0;
	}
	return gamma;
}
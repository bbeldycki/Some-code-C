#include "integrals.h"

Integrals::Integrals()
{
	std::unique_ptr<CarlsonFunction> cFunctions = std::make_unique<CarlsonFunction>();
	std::unique_ptr<IntegralsHelper> intHelper = std::make_unique<IntegralsHelper>();
}

Integrals::~Integrals()
{

}

double Integrals::ellIntCubicAllRootsReal(std::vector<int>& pList, std::vector<double>& aList, std::vector<double>& bList, double& ffr, double& y, double& x)
{
	double dOnTw, dOnTr, dOnFo, dTwFo, dTrFo;
	double uOne, uTwo, uThree, wTwoTwo, qTwoTwo, pTwoTwo;
	double handlerX, handlerY;
	double iTwC, iTrC, rOnTw, rOnTr, rTwFo, rTrFo, kTwC;
	std::vector<double> xi, yi;
	try
	{
		if (x < y)
		{
			throw std::runtime_error("Error: x is lower than y. ellIntCubicAllRootsReal.");
		}
	}
	catch (const std::exception& e)
	{
		std::cerr << "Exception caught: " << e.what() << std::endl;
		return 0;
	}
	dOnTw = aList[0] * bList[1] - aList[1] * bList[0];
	dOnTr = aList[0] * bList[2] - aList[2] * bList[0];
	dOnFo = aList[0] * bList[3] - aList[3] * bList[0];
	dTwFo = aList[1] * bList[3] - aList[3] * bList[1];
	dTrFo = aList[2] * bList[4] - aList[4] * bList[2];
	for (int i = 0; i < aList.size(); i++)
	{
		handlerX = aList[i] + bList[i] * x;
		handlerY = aList[i] + bList[i] * y;
		try
		{
			if (handlerX < 0)
			{
				throw std::runtime_error("Error: handlerX is lower than 0. ellIntCubicAllRootsReal.");
			}
			if (handlerY < 0)
			{
				throw std::runtime_error("Error: handlerX is lower than 0. ellIntCubicAllRootsReal.");
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
	uOne = (xi[0] * yi[1] * yi[2] + yi[0] * xi[1] * xi[2]) / (x - y);
	uTwo = (xi[1] * yi[0] * yi[2] + yi[1] * xi[0] * xi[2]) / (x - y);
	uThree = (xi[2] * yi[0] * yi[1] + yi[2] * xi[0] * xi[1]) / (x - y);
	try
	{
		if (bList[0] == 0 || bList[1] == 0 || bList[2] == 0 || bList[3] == 0 || dOnFo == 0 || dTwFo == 0 || dTrFo == 0 || uOne == 0 || xi[0] == 0 || yi[0] == 0)
		{
			throw std::runtime_error("Error: All the numbers in if statement will appear in denominator later on. We must be sure they are not 0. ellIntCubicAllRootsReal.");
		}
	}
	catch (const std::exception& e)
	{
		std::cerr << "Exception caught: " << e.what() << std::endl;
		return 0;
	}
	wTwoTwo = pow(uOne, 2) - bList[3] * dOnTw * dOnTr / dOnFo;
	qTwoTwo = wTwoTwo * pow((xi[3] * yi[3] / xi[0] / yi[0]), 2);
	pTwoTwo = qTwoTwo + bList[3] * dTrFo * dTrFo / dOnFo;

	if (pList[3] == 0)
	{
		return 2.0 * cFunctions->rcFunction(pTwoTwo, qTwoTwo);
	}
	else
	{
		iTrC = -1.0 * bList[0] * dOnTw * dOnTr * cFunctions->rjFunction(pow(uThree, 2), pow(uTwo, 2), pow(uOne, 2), wTwoTwo) / 3.0 / dOnFo;
		iTrC = 2.0 * cFunctions->rcFunction(pTwoTwo, qTwoTwo) + iTrC;
		if (pList[3] == -2)
		{
			return (bList[3] * iTrC - bList[0] * 2.0 * cFunctions->rcFunction(pTwoTwo, qTwoTwo)) / dOnFo;
		}
		rOnTw = dOnTw / bList[0] / bList[1];
		rOnTr = dOnTr / bList[0] / bList[2];
		rTwFo = dTwFo / bList[1] / bList[3];
		rTrFo = dTrFo / bList[2] / bList[3];
		iTwC = 2.0 * dOnTw * dOnTr * cFunctions->rdFunction(pow(uThree, 2), pow(uTwo, 2), pow(uOne, 2)) / 3.0;
		iTwC = 2.0 * xi[0] * yi[0] / uOne + iTwC;
		kTwC = bList[1] * bList[2] * iTwC - 2.0 * bList[3] * intHelper->computeAPOnePn(pList, aList, bList, y, x);
		return (bList[3] * kTwC / (2.0 * dOnFo * dTwFo * dTrFo) + 2.0 * cFunctions->rcFunction(pTwoTwo, qTwoTwo) * (1.0 - rOnTw * rOnTr / 2.0 / rTwFo / rTrFo) * pow(bList[0] / dOnFo, 2));
	}
}

double Integrals::ellIntCubicOneRealTwoComplexRoots(std::vector<int>& pList, std::vector<double>& aList, std::vector<double>& bList, std::vector<double>& fghList, double& ffr, double& y, double& x)
{
	double dOnFo, beta, cOnOnTw, cOnFoTw, cFoFoTw, ksi, eta, mTw;
	double lmTw, lpTw, wpTw, uTw, wTw, qTw, pTw, rho;
	double handlerX, handlerY;
	double iTrC;
	double kTwC, nTwC, amOnOnOnmTw, rTwFo, rOnTw;
	std::vector<double> xi, yi;
	try
	{
		if (x < y)
		{
			throw std::runtime_error("Error: x is lower than y. ellIntCubicOneRealTwoComplexRoots.");
		}
	}
	catch (const std::exception& e)
	{
		std::cerr << "Exception caught: " << e.what() << std::endl;
		return 0;
	}
	dOnFo = aList[0] * bList[3] - aList[3] * bList[0];
	for (int i = 0; i < aList.size(); i++)
	{
		handlerX = aList[i] + bList[i] * x;
		handlerY = aList[i] + bList[i] * y;
		try
		{
			if (handlerX < 0)
			{
				throw std::runtime_error("Error: handlerX is lower than 0. ellIntCubicOneRealTwoComplexRoots.");
			}
			if (handlerY < 0)
			{
				throw std::runtime_error("Error: handlerX is lower than 0. ellIntCubicOneRealTwoComplexRoots.");
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
	try
	{
		if (xi[0] == 0 || yi[0] == 0)
		{
			throw std::runtime_error("Error: Later on both xi[0] and yi[0] will appear in the denominator so we must be sure they are not equal to 0. ellIntCubicOneRealTwoComplexRoots.");
		}
	}
	catch (const std::exception& e)
	{
		std::cerr << "Exception caught: " << e.what() << std::endl;
		return 0;
	}
	beta = fghList[1] * bList[0] - 2.0 * fghList[2] * aList[0];
	cOnOnTw = intHelper->computeCijValue(aList, bList, fghList, 0, 0);
	cOnFoTw = intHelper->computeCijValue(aList, bList, fghList, 0, 3);
	cFoFoTw = intHelper->computeCijValue(aList, bList, fghList, 3, 3);
	ksi = sqrt(fghList[0] + fghList[1] * x + fghList[2] * pow(x, 2));
	eta = sqrt(fghList[0] + fghList[1] * y + fghList[2] * pow(y, 2));
	mTw = pow((xi[0] + yi[0]) * sqrt(pow(ksi + eta, 2) - fghList[2] * pow(x - y, 2)) / (x - y), 2);
	lmTw = mTw - beta - sqrt(2.0 * fghList[2] * cOnOnTw);
	lpTw = mTw - beta + sqrt(2.0 * fghList[2] * cOnOnTw);

	if (pList[3] == 0)
	{
		return 4.0 * cFunctions->rfFunction(mTw, lmTw, lpTw);
	}
	else
	{
		wpTw = mTw - bList[0] * (cOnFoTw + sqrt(cOnOnTw * cFoFoTw)) / dOnFo;
		uTw = pow((xi[0] * eta + yi[0] * ksi) / (x - y), 2);
		wTw = uTw - cOnOnTw * bList[3] / 2.0 / dOnFo;
		qTw = wTw * pow(xi[3] * yi[3] / xi[0] / yi[0], 2);
		pTw = qTw + cFoFoTw * bList[3] / 2.0 / dOnFo;
		rho = sqrt(2.0 * fghList[2] * cOnOnTw) - beta;
		if (pList[3] == -2)
		{
			iTrC = 2.0 * cFunctions->rcFunction(pTw, qTw) + 3.0 * cFunctions->rcFunction(uTw, wTw) - 6.0 * cFunctions->rfFunction(mTw, lmTw, lpTw);
			iTrC = (2.0 * sqrt(cOnOnTw) / 3.0 / sqrt(cFoFoTw)) * ((-4.0 * bList[0] / dOnFo) * (cOnFoTw + sqrt(cOnOnTw * cFoFoTw))) * cFunctions->rjFunction(mTw, lmTw, lpTw, wpTw) + iTrC;
			return (bList[3] * iTrC - bList[0] * ffr) / dOnFo;
		}
		else if (pList[3] == -4)
		{
			amOnOnOnmTw = intHelper->computeAPOnePnQuadratic(pList, aList, bList, fghList, y, x);
			nTwC = sqrt(8.0 * fghList[2] / 9.0 / cOnOnTw) * (4.0 * rho * cFunctions->rdFunction(mTw, lmTw, lpTw) - 6.0 * cFunctions->rfFunction(mTw, lmTw, lpTw) + 3.0 / sqrt(uTw) + 2.0 / (xi[0] * yi[0] * sqrt(uTw)));
			kTwC = cOnOnTw * nTwC / 2.0 - 2.0 * dOnFo * amOnOnOnmTw;
			rTwFo = cFoFoTw / 2.0 / pow(bList[3], 2) / fghList[2];
			rOnTw = cOnOnTw / 2.0 / pow(bList[0], 2) / fghList[2];
			return bList[3] * kTwC / 4.0 / dOnFo / cFoFoTw + pow(bList[0] / dOnFo, 2) * (1.0 - rOnTw / rTwFo) * ffr;
		}
	}
	return 0;
}
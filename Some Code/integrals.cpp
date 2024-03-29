#include "integrals.h"

Integrals::Integrals()
{
	cFunctions = std::make_unique<CarlsonFunction>();
	intHelper = std::make_unique<IntegralsHelper>();
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

double Integrals::ellIntQuarticAllRootsReal(std::vector<int>& pList, std::vector<double>& aList, std::vector<double>& bList, double& ffr, double& y, double& x)
{
	double dOnTw, dOnTr, dOnFo, dOnFv;
	double uOnTw, uOnTr, uOnFo;
	double wOnTwTw, qOnTwTw, pOnTwTw, iTrP;
	double handlerX, handlerY;
	double dTwFv, dTrFv, dFoFv, wTwTw, qTwTw, pTwTw;
	double iTwo, iThree;
	double dTwFo, dTwFv, dTrFo, dTrFv, dFoFv, rOnTw, rOnTr, rTwFv, rTrFv;
	double aOnOnOnmOnmTw;
	double result = 0;
	std::vector<double> xi, yi;
	try
	{
		if (x < y)
		{
			throw std::runtime_error("Error: x is lower than y. ellIntQuarticAllRootsReal.");
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
	dOnFv = aList[0] * bList[4] - aList[4] * bList[0];
	for (int i = 0; i < aList.size(); i++)
	{
		handlerX = aList[i] + bList[i] * x;
		handlerY = aList[i] + bList[i] * y;
		try
		{
			if (handlerX < 0)
			{
				throw std::runtime_error("Error: handlerX is lower than 0. ellIntQuarticAllRootsReal.");
			}
			if (handlerY < 0)
			{
				throw std::runtime_error("Error: handlerX is lower than 0. ellIntQuarticAllRootsReal.");
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
	uOnTw = (xi[0] * yi[1] * yi[2] * yi[3] + yi[0] * yi[1] * xi[2] * xi[3]) / (x - y);
	uOnTr = (xi[1] * yi[2] * yi[1] * yi[2] + yi[0] * yi[2] * xi[1] * xi[3]) / (x - y);
	uOnFo = (xi[1] * yi[3] * yi[1] * yi[1] + yi[0] * yi[3] * xi[1] * xi[2]) / (x - y);
	try
	{
		if (bList[0] == 0 || bList[1] == 0 || bList[2] == 0 || bList[4] == 0 || dOnFv == 0 || uOnFo == 0 || xi[0] == 0 || yi[0] == 0 || xi[3] == 0 || yi[3] == 0)
		{
			throw std::runtime_error("Error: All the numbers in if statement will appear in denominator later on. We must be sure they are not 0. ellIntQuarticAllRootsReal.");
		}
	}
	catch (const std::exception& e)
	{
		std::cerr << "Exception caught: " << e.what() << std::endl;
		return 0;
	}

	if (aList[4] == 1.0 && bList[4] == 0.0)
	{
		wOnTwTw = pow(uOnTw, 2) - bList[1] * dOnTr * dOnFo / bList[0];
		qOnTwTw = wOnTwTw / pow(xi[0] * yi[0], 2);
		pOnTwTw = qOnTwTw + bList[1] * bList[2] * bList[3] / bList[0];
		iTrP = -2.0 * dOnTw * dOnTr * dOnFo * cFunctions->rjFunction(pow(uOnTw, 2), pow(uOnTr, 2), pow(uOnFo, 2), wOnTwTw) / 3.0 / bList[0] + 2.0 * cFunctions->rcFunction(pOnTwTw, qOnTwTw);
		if (pList[4] == 0)
		{
			return 2.0 * cFunctions->rcFunction(pOnTwTw, qOnTwTw);
		}
		else if (pList[4] == -2)
		{
			return (bList[4] * iTrP - bList[0] * 2.0 * ffr) / dOnFv;
		}
	}
	else
	{
		dTwFv = aList[1] * bList[4] - aList[4] * bList[1];
		dTrFv = aList[2] * bList[4] - aList[4] * bList[2];
		dFoFv = aList[3] * bList[4] - aList[4] * bList[3];
		wTwTw = pow(uOnTw, 2) - dOnTr * dOnFo * dTwFv / dOnFv;
		qTwTw = wTwTw * pow(xi[4] * yi[4] / xi[0] / yi[0], 2);
		pTwTw = qTwTw + dTwFv * dTrFv * dFoFv / dOnFv;
		iTwo = 2.0 * dOnTw * dOnTr * cFunctions->rdFunction(pow(uOnTw, 2), pow(uOnTr, 2), pow(uOnFo, 2)) / 2.0 + 2.0 * xi[0] * yi[0] / xi[3] / yi[3] / uOnFo;
		iThree = 2.0 * dOnTw * dOnTr * dOnFo * cFunctions->rjFunction(pow(uOnTw, 2), pow(uOnTr, 2), pow(uOnFo, 2), wTwTw) / 3.0 / dOnFv + 2.0 * cFunctions->rcFunction(pTwTw, qTwTw);
		if (pList[4] == 0)
		{
			return 2.0 * cFunctions->rfFunction(pow(uOnTw, 2), pow(uOnTr, 2), pow(uOnFo, 2));
		}
		else if (pList[4] == -2)
		{
			return (bList[4] * iThree - bList[0] * 2.0 * ffr) / dOnFv;
		}
		else if(pList[4] == -4)
		{
			dTwFo = aList[1] * bList[3] - aList[3] * bList[1];
			dTwFv = aList[1] * bList[4] - aList[4] * bList[1];
			dTrFo = aList[2] * bList[3] - aList[3] * bList[2];
			dTrFv = aList[2] * bList[4] - aList[4] * bList[2];
			dFoFv = aList[3] * bList[4] - aList[4] * bList[3];
			rOnTw = dOnTw / bList[0] / bList[1];
			rOnTr = dOnTr / bList[0] / bList[2];
			rTwFv = dTwFv / bList[1] / bList[4];
			rTrFv = dTrFv / bList[2] / bList[4];
			aOnOnOnmOnmTw = intHelper->computeAPOnePn(pList, aList, bList, y, x);
			result = pow(bList[4], 2) * dTwFo * dTrFo * iTwo / (2.0 * dOnFv * dTwFv * dTrFv * dFoFv);
			result = pow(bList[0], 2) * 2.0 * ffr * (1.0 - rOnTw * rOnTr / (2.0 * rTwFv * rTrFv)) + result;
			result = result - pow(bList[4], 2) * aOnOnOnmOnmTw / dOnFv / dTwFv / dTrFv;
			return result;
		}
	}
}

double Integrals::ellIntQuarticTwoRealTwoComplexRoots(std::vector<int>& pList, std::vector<double>& aList, std::vector<double>& bList, std::vector<double>& fghList, double& ffr, double& y, double& x)
{
	double dOnFo, dOnFv, dFoFv, cOnOn, cOnFo, cOnFv, cFoFo, cFvFv, ksi, eta, mTw, lmTw, lpTw, wpTw, uTw;
	double wOnTw, rho, mTwRho, qOnTw, pOnTw, iTrC;
	double handlerX, handlerY;
	double wTw, qTw, pTw, iTwo, iThree, result;
	std::vector<double> xi, yi;
	try
	{
		if (x < y)
		{
			throw std::runtime_error("Error: x is lower than y. ellIntQuarticTwoRealTwoComplexRoots.");
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
				throw std::runtime_error("Error: handlerX is lower than 0. ellIntQuarticTwoRealTwoComplexRoots.");
			}
			if (handlerY < 0)
			{
				throw std::runtime_error("Error: handlerX is lower than 0. ellIntQuarticTwoRealTwoComplexRoots.");
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
	dOnFo = aList[0] * bList[3] - aList[3] * bList[0];
	dOnFv = aList[0] * bList[4] - aList[4] * bList[0];
	dFoFv = aList[3] * bList[4] - aList[4] * bList[3];
	cOnOn = intHelper->computeCijValue(aList, bList, fghList, 0, 0);
	cOnFo = intHelper->computeCijValue(aList, bList, fghList, 0, 3);
	cOnFv = intHelper->computeCijValue(aList, bList, fghList, 0, 4);
	cFoFo = intHelper->computeCijValue(aList, bList, fghList, 3, 3);
	cFvFv = intHelper->computeCijValue(aList, bList, fghList, 4, 4);
	try
	{
		if (dOnFv == 0.0 || bList[0] == 0.0 || fghList[2] == 0.0 || xi[0] == 0.0 || yi[0] == 0.0 || xi[3] == 0.0 || yi[3] == 0.0 || cFoFo == 0.0 || cFvFv == 0.0)
		{
			throw std::runtime_error("Error: Later on a number of variables will appear in the denominator so we must be sure they are not equal to 0. ellIntCubicOneRealTwoComplexRoots.");
		}
	}
	catch (const std::exception& e)
	{
		std::cerr << "Exception caught: " << e.what() << std::endl;
		return 0;
	}
	ksi = sqrt(fghList[0] + fghList[1] * x + fghList[2] * pow(x, 2));
	eta = sqrt(fghList[0] + fghList[1] * y + fghList[2] * pow(y, 2));
	mTw = pow((xi[0] * yi[3] + yi[0] * xi[3]) * sqrt(pow(ksi + eta, 2) - fghList[2] * pow(x - y, 2)) / (x - y), 2);
	lmTw = mTw + cOnFo - sqrt(cOnOn * cFoFo);
	lpTw = mTw + cOnFo + sqrt(cOnOn * cFoFo);
	wpTw = mTw + dOnFo * (cOnFv * sqrt(cOnOn * cFvFv)) / dOnFv;
	uTw = pow((xi[0] * xi[3] * eta + yi[0] * yi[3] * ksi) / (x - y), 2);
	try
	{
		if (uTw == 0.0)
		{
			throw std::runtime_error("Error: Later on u2 will appear in the denominator so we must be sure they are not equal to 0. ellIntCubicOneRealTwoComplexRoots.");
		}
	}
	catch (const std::exception& e)
	{
		std::cerr << "Exception caught: " << e.what() << std::endl;
		return 0;
	}
	if (aList[4] == 1.0 and bList[4] == 0.0)
	{
		wOnTw = uTw - cOnOn * bList[3] / 2.0 / bList[0];
		rho = dOnFo * (fghList[1] * bList[0] - 2.0 * fghList[2] * aList[0] - sqrt(2.0 * fghList[2] * cOnOn)) / bList[0];
		mTwRho = mTw + rho;
		qOnTw = wOnTw / pow(xi[0], 2) / pow(yi[0], 2);
		pOnTw = qOnTw + fghList[2] * bList[3] / bList[0];
		if (pList[4] == 0)
		{
			return 4.0 * cFunctions->rfFunction(mTw, lmTw, lpTw);
		}
		else if (pList[4] == -2)
		{
			iTrC = sqrt(2.0 * cOnOn / 9.0 / fghList[2]) * (4.0 * rho * cFunctions->rjFunction(mTw, lmTw, lpTw, mTwRho) - 6.0 * cFunctions->rfFunction(mTw, lmTw, lpTw) + 3.0 * cFunctions->rcFunction(uTw, wOnTw));
			iTrC = 2.0 * cFunctions->rcFunction(pOnTw, qOnTw) + iTrC;
			return (bList[4] * iTrC - bList[0] * ffr) / dOnFv;
		}
	}
	else
	{
		wTw = uTw - cOnOn * dFoFv / 2.0 / dOnFv;
		qTw = wTw * pow(xi[4] * yi[4] / xi[0] / yi[0], 2);
		pTw = qTw + cFvFv * dFoFv / 2.0 / dOnFv;
		if (pList[4] == 0)
		{
			return 4.0 * cFunctions->rfFunction(mTw, lmTw, lpTw);
		}
		else if (pList[4] == -2)
		{
			iThree = (2.0 * sqrt(cOnOn) / 3.0 / sqrt(cFvFv));
			iThree = (4.0 * dOnFo / dOnFv * (cOnFv + sqrt(cOnOn * cFvFv)) * cFunctions->rjFunction(mTw, lmTw, lpTw, wpTw) - 6.0 * cFunctions->rfFunction(mTw, lmTw, lpTw) + 3.0 * cFunctions->rcFunction(uTw, wTw)) * iThree;
			iThree = 2.0 * cFunctions->rcFunction(pTw, qTw) + iThree;
			return (bList[4] * iThree - bList[0] * ffr) / dOnFv;
		}
		else if (pList[4] == -4)
		{
			iTwo = (2.0 * sqrt(cOnOn) / 3.0 / sqrt(cFoFo)) * (4.0 * (cOnFo + sqrt(cOnOn * cFoFo)) * cFunctions->rdFunction(mTw, lmTw, lpTw) - 6.0 * cFunctions->rfFunction(mTw, lmTw, lpTw) + 3.0 / sqrt(uTw));
			iTwo = 2.0 * xi[0] * yi[0] / xi[3] / yi[3] / sqrt(uTw) + iTwo;
			result = pow(bList[4], 2) * cFoFo * iTwo / 2.0 / dOnFv / dFoFv / cFvFv + pow(bList[0], 2) * ffr * (1.0 - cOnOn * pow(bList[4], 2) / 2.0 / cFvFv / pow(bList[0], 2));
			result = result - 2.0 * pow(bList[4], 2) * intHelper->computeAPOnePnQuadratic(pList, aList, bList, fghList, y, x) / dOnFv / cFvFv;
			return result;
		}
	}
}

double Integrals::ellIntQuarticAllComplexRoots(std::vector<int>& pList, std::vector<double>& aList, std::vector<double>& bList, std::vector<double>& fghOneList, std::vector<double>& fghTwoList, double& ffr, double& y, double& x)
{;
	double handlerKsi, handlerEta;
	std::vector<double> fi, gi, hi, ksi, eta, theta, dzeta;
	double ksiOnep, etaOnep, u, m, dOnOn, dTwTw, dOnTw, delta, deltaP, deltaAm, lmTw, lpTw;
	double gie, alfaOnFv, alfaTwFv, betaOnFv, betaTwFv, gammaOn, gammaTw;
	double bLambda, bSigTw, psi, bix, bes, mu, bet, vTw, bTw, aTw, beh, omega;
	double ksiFv, etaFv;
	double result;
	try
	{
		if (x < y)
		{
			throw std::runtime_error("Error: x is lower than y. ellIntQuarticAllComplexRoots.");
		}
	}
	catch (const std::exception& e)
	{
		std::cerr << "Exception caught: " << e.what() << std::endl;
		return 0;
	}
	
	fi.push_back(fghOneList[0]);
	fi.push_back(fghTwoList[0]);
	gi.push_back(fghOneList[1]);
	gi.push_back(fghTwoList[1]);
	hi.push_back(fghOneList[2]);
	hi.push_back(fghTwoList[2]);
	for (int i = 0; i < fi.size(); i++)
	{
		handlerKsi = fi[i] + gi[i] * x + hi[i] * pow(x, 2);
		handlerEta = fi[i] + gi[i] * y + hi[i] * pow(y, 2);
		try
		{
			if (handlerKsi < 0)
			{
				throw std::runtime_error("Error: handlerX is lower than 0. ellIntQuarticAllComplexRoots.");
			}
			if (handlerEta < 0)
			{
				throw std::runtime_error("Error: handlerX is lower than 0. ellIntQuarticAllComplexRoots.");
			}
		}
		catch (const std::exception& e)
		{
			std::cerr << "Exception caught: " << e.what() << std::endl;
			return 0;
		}
		ksi.push_back(sqrt(handlerKsi));
		eta.push_back(sqrt(handlerEta));
		theta.push_back(handlerKsi + handlerEta - hi[i] * pow(x - y, 2));
	}
	ksiOnep = (gi[0] + 2.0 * hi[0] * x) / 2.0 / ksi[0];
	etaOnep = (gi[0] + 2.0 * hi[0] * y) / 2.0 / eta[0];
	for (int i = 0; i < fi.size(); i++)
	{
		handlerKsi = pow(ksi[i] + eta[i], 2) - 2.0 * hi[0] * pow(x - y, 2);
		try
		{
			if (handlerKsi < 0)
			{
				throw std::runtime_error("Error: handlerX is lower than 0. ellIntQuarticAllComplexRoots.");
			}
		}
		catch (const std::exception& e)
		{
			std::cerr << "Exception caught: " << e.what() << std::endl;
			return 0;
		}
		dzeta.push_back(sqrt(handlerKsi));
	}
	u = (ksi[0] * eta[1] + eta[0] * ksi[0]) / (x - y);
	m = dzeta[0] * dzeta[1] / (x - y);
	dOnOn = intHelper->computeDeltaij(fi, gi, hi, 0, 0);
	dTwTw = intHelper->computeDeltaij(fi, gi, hi, 1, 1);
	dOnTw = intHelper->computeDeltaij(fi, gi, hi, 0, 1);
	delta = sqrt(pow(dOnTw, 2) - dOnOn * dTwTw);
	deltaP = dOnTw + delta;
	deltaAm = dOnTw - delta;
	lmTw = pow(m, 2) + deltaAm;
	lpTw = pow(m, 2) + deltaP;
	if (pList[4] == 0)
	{
		return 4.0 * cFunctions->rfFunction(pow(m, 2), lmTw, lpTw);
	}
	try
	{
		if (u == 0.0)
		{
			throw std::runtime_error("Error: u is equal to 0 and we will use it in denominator. ellIntQuarticAllComplexRoots.");
		}
	}
	catch (const std::exception& e)
	{
		std::cerr << "Exception caught: " << e.what() << std::endl;
		return 0;
	}
	gie = 2.0 * delta * deltaP * cFunctions->rdFunction(pow(m, 2), lmTw, lpTw) / 3.0 + delta / 2.0 / u + (dOnTw * theta[0] - dOnOn * theta[1]) / 4.0 / ksi[0] / eta[0] / u;
	alfaOnFv = intHelper->computeCoeffAlfai(fi, gi, aList[4], bList[4], 0);
	betaOnFv = intHelper->computeCoeffBetai(gi, hi, aList[4], bList[4], 0);
	alfaTwFv = intHelper->computeCoeffAlfai(fi, gi, aList[4], bList[4], 1);
	betaTwFv = intHelper->computeCoeffBetai(gi, hi, aList[4], bList[4], 1);
	gammaOn = intHelper->ComputeGammai(fi, gi, hi, aList[4], bList[4], 0);
	gammaTw = intHelper->ComputeGammai(fi, gi, hi, aList[4], bList[4], 1);
	if (aList[4] == 1.0 && bList[4] == 0.0)
	{
		try
		{
			if (hi[0] == 0.0)
			{
				throw std::runtime_error("Error: Hi[0] is equal to 0 and we will use it in denominator. ellIntQuarticAllComplexRoots.");
			}
		}
		catch (const std::exception& e)
		{
			std::cerr << "Exception caught: " << e.what() << std::endl;
			return 0;
		}
	}
	bLambda = dOnOn * hi[1] / hi[0];
	bSigTw = pow(m, 2) + bLambda;
	psi = gi[0] * hi[1] - gi[1] * hi[0];
	bix = -1.0 * (ksiOnep * ksi[1] + etaOnep * eta[1]) / (x - y);
	bes = (pow(m, 2) + dOnOn) / 2.0 - pow(u, 2);
	mu = hi[0] / ksi[0] / eta[0];
	bet = mu * bes + 2.0 * hi[0] * hi[1];
	vTw = pow(mu, 2) * (pow(bes, 2) + bLambda * pow(u, 2));
	bTw = (pow(bes / u, 2) + bLambda) * bSigTw * pow(u, 2);
	aTw = bTw + bLambda * (deltaP - bLambda) * (bLambda - deltaAm);
	beh = dOnOn * psi * (cFunctions->rjFunction(pow(m, 2), lmTw, lpTw, bSigTw) / 3.0 + cFunctions->rcFunction(aTw, bTw) / 2.0) / pow(hi[0], 2) - bix * cFunctions->rcFunction(pow(bet, 2), vTw);
	if (pList[4] == -2)
	{
		return -2.0 * (bList[4] * beh + betaOnFv * ffr / gammaOn);
	}
	else if (pList[4] = -4)
	{
		omega = gie - deltaP * ffr + ksiOnep * ksi[0] - etaOnep * eta[0];
		result = bList[4] * (betaOnFv / gammaOn + betaTwFv / gammaTw) * beh + pow(betaOnFv / gammaOn, 2) * ffr;
		result = result + pow(bList[4], 2) * (omega - bList[4] * intHelper->computeAPOnePnDoubleQuadratic(pList, aList, bList, fghOneList, fghTwoList, y, x)) / gammaOn / gammaTw;
		return result;
	}
	else
	{
		bLambda = dOnOn * gammaTw / gammaOn;
		bSigTw = pow(m, 2) + bLambda;
		psi = (alfaOnFv * betaTwFv - alfaTwFv * betaOnFv) / 2.0;
		ksiFv = aList[4] + bList[4] * x;
		etaFv = aList[4] + bList[4] * y;
		bix = (ksiFv * (alfaOnFv + betaOnFv * y) * eta[1] / eta[0] + etaFv * (alfaOnFv + betaOnFv * x) * ksi[1] / ksi[0]) / 2.0 / (x - y);
		bes = (pow(m, 2) + dOnTw) / 2.0 - pow(u, 2);
		mu = gammaOn * ksiFv * etaFv / ksi[0] / eta[0];
		bet = mu * bes + 2.0 * gammaOn * gammaTw;
		vTw = pow(mu, 2) * (pow(bes, 2) + bLambda * pow(u, 2));
		bTw = (pow(bes / u, 2) + bLambda) * pow(bSigTw, 2);
		aTw = bTw + bLambda * (deltaP - bLambda) * (bLambda - deltaAm);
		beh = dOnOn * psi * (cFunctions->rjFunction(pow(m, 2), lmTw, lpTw, bSigTw) / 3.0 + cFunctions->rcFunction(aTw, bTw) / 2.0) / pow(gammaOn, 2) - bix * cFunctions->rcFunction(pow(bet, 2), vTw);
		if (pList[4] == -2)
		{
			return -2.0 * (bList[4] * beh + betaOnFv * ffr / gammaOn);
		}
		else if (pList[4] == -4)
		{
			omega = gie - deltaP * ffr + ksiOnep * ksi[0] - etaOnep * eta[0];
			result = bList[4] * (betaOnFv / gammaOn + betaTwFv / gammaTw) * beh + pow(betaOnFv / gammaOn, 2) * ffr;
			result = result + pow(bList[4], 2) * (omega - bList[4] * intHelper->computeAPOnePnDoubleQuadratic(pList, aList, bList, fghOneList, fghTwoList, y, x)) / gammaOn / gammaTw;
		}
	}
}
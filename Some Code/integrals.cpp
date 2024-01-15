#include "integrals.h"
#include "integralshelper.h"
#include "carlfunctions.h"

Integrals::Integrals()
{
	std::unique_ptr<CarlsonFunction> cFunctions = std::make_unique<CarlsonFunction>();
	std::unique_ptr<IntegralsHelper> intHelper = std::make_unique<IntegralsHelper>();
}

Integrals::~Integrals()
{

}

double_t Integrals::ellIntCubicAllRootsReal(std::vector<int>& pList, std::vector<double_t>& aList, std::vector<double_t>& bList, double_t& ffr, double_t& y, double_t& x)
{
	double_t dOnTw, dOnTr, dOnFo, dTwFo, dTrFo;
	double_t uOne, uTwo, uThree, wTwoTwo, qTwoTwo, pTwoTwo;
	double_t handlerX, handlerY;
	double_t iTwC, iTrC, rOnTw, rOnTr, rTwFo, rTrFo, kTwC;
	std::vector<double_t> xi, yi;
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
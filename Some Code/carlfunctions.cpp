#include "carlfunctions.h"

CarlsonFunction::CarlsonFunction()
{

}

CarlsonFunction::~CarlsonFunction()
{

}

void CarlsonFunction::iterationForRCFunction(double_t& x, double_t& y, double_t& ave, double_t& sol)
{
	x = 0.25 * (x + 2.0 * sqrt(x * y) + y);
	y = 0.25 * (y + 2.0 * sqrt(x * y) + y);
	ave = (x + 2.0 * y) / 3.0;
	sol = (y - ave) / ave;
}

double_t CarlsonFunction::rcFunction(const double_t& x, const double_t& y)
{
	double_t newX, newY, val, ave, so, result;
	try
	{
		if (x < 0)
		{
			throw std::runtime_error("Error: x is lower than 0. rcFunction.");
		}
		else if (y == 0)
		{
			throw std::runtime_error("Error: y cannot be equal to 0. rcFunction.");
		}
		else if (x + fabs(y) < 1.69e-38)
		{
			throw std::runtime_error("Error: Both x and y are to close to 0. rcFunction.");
		}
		else if (x + fabs(y) > 3e37)
		{
			throw std::runtime_error("Error: Sum of x and modile y is to large. rcFunction.");
		}
		else if (y < -1.7e19 || 0.0 < x < 0.01)
		{
			throw std::runtime_error("Error: y is a large negative number and x is a small number. rcFunction.");
		}
	}
	catch (const std::exception& e)
	{
		std::cerr << "Exception caught: " << e.what() << std::endl;
		return 0;
	}
	result = 0;
	ave = 0;
	so = 0;
	if (y > 0.0)
	{
		val = 1.0;
		newX = x;
		newY = y;
	}
	else
	{
		newX = x - y;
		newY = -1.0 * y;
		val = sqrt(x / (x - y));
	}
	// First iteration
	iterationForRCFunction(newX, newY, ave, so);
	while (fabs(so) > 0.04)
	{
		iterationForRCFunction(newX, newY, ave, so);
	}
	// Compute rcFunction
	result = 0.375 + 9.0 * so / 22.0;
	result = 1.0 / 7.0 + so * result;
	result = 0.3 + so * result;
	result = 1.0 * pow(so, 2) * result;
	result = val * result / sqrt(ave);
	return result;
}

void CarlsonFunction::iterationForRDFunction(double_t& x, double_t& y, double_t& z, double_t& dX, double_t& dY, double_t& dZ, double_t& ave, double_t& sol, double_t& f)
{
	double_t alamb = 0;
	alamb = sqrt(x * y) + sqrt(x * z) + sqrt(y * z);
	sol = sol + f / (sqrt(z) * (z + alamb));
	f = 0.25 * f;
	x = 0.25 * (x + alamb);
	y = 0.25 * (y + alamb);
	z = 0.25 * (z + alamb);
	ave = 0.2 * (x + y + 3.0 * z);
	dX = (ave - x) / ave;
	dY = (ave - y) / ave;
	dZ = (ave - z) / ave;
}


double_t CarlsonFunction::rdFunction(double_t x, double_t y, double_t z)
{
	double_t dX, dY, dZ, ave, so, f, result;
	double_t c2, c3, c4, c5;
	try
	{
		if (std::min(x, y) < 0)
		{
			throw std::runtime_error("Error: One or more argument is lower than 0. rdFunction.");
		}
		else if (std::min(x + y, z) < 1.5e-25)
		{
			throw std::runtime_error("Error: Sum of x + y or value of z was to close to 0. rdFunction.");
		}
		else if (std::max(x, y, z) > 4.5e25)
		{
			throw std::runtime_error("Error: One of the argument is to large. rdFunction.");
		}
	}
	catch (const std::exception& e)
	{
		std::cerr << "Exception caught: " << e.what() << std::endl;
		return 0;
	}
	dX = dY = dZ = ave = so = result = 0;
	f = 1.0;
	// First iteration
	iterationForRDFunction(x, y, z, dX, dY, dZ, ave, so, f);
	while (std::max(fabs(x), fabs(y), fabs(z)) > 0.05)
	{
		iterationForRDFunction(x, y, z, dX, dY, dZ, ave, so, f);
	}
	c2 = dX * dY;
	c3 = dX * dY - pow(dZ, 2);
	c4 = dX * dY - 6.0 * pow(dZ, 2);
	c5 = c4 + 2.0 * c3;
	// Compute rcFunction
	result = -9.0 * c3 / 22.0 + 3.0 * dZ * c2 / 26.0;
	result = c5 / 6.0 + dZ * result;
	result = 1.0 + dZ * result;
	result = c4 * ( -.30 / 14.0 + 9.0 * c4 - 45.0 * c5 * dZ / 260.0) + result;
	result = f * result;
	result = 3.0 * so + result / (ave * sqrt(ave));
	return result;
}
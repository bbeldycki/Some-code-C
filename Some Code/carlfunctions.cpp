#include "carlfunctions.h"

CarlsonFunction::CarlsonFunction()
{

}

CarlsonFunction::~CarlsonFunction()
{

}

void CarlsonFunction::iterationForRCFunction(double& x, double& y, double& ave, double& sol)
{
	x = 0.25 * (x + 2.0 * sqrt(x * y) + y);
	y = 0.25 * (y + 2.0 * sqrt(x * y) + y);
	ave = (x + 2.0 * y) / 3.0;
	sol = (y - ave) / ave;
}

double CarlsonFunction::rcFunction(double x, double y)
{
	double newX, newY, val, ave, so, result;
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
		else if (y < -1.7e19 || (0.0 < x && x < 0.1))
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

void CarlsonFunction::iterationForRDFunction(double& x, double& y, double& z, double& dX, double& dY, double& dZ, double& ave, double& sol, double& f)
{
	double alamb = 0;
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


double CarlsonFunction::rdFunction(double x, double y, double z)
{
	double dX, dY, dZ, ave, so, f, result;
	double c2, c3, c4, c5;
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
		else if (std::max(x, std::max(y, z)) > 4.5e25)
		{
			throw std::runtime_error("Error: One of the arguments is to large. rdFunction.");
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
	while (std::max(fabs(x), std::max(fabs(y), fabs(z))) > 0.05)
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

void CarlsonFunction::iterationForRFFunction(double& x, double& y, double& z, double& dX, double& dY, double& dZ, double& ave)
{
	double alamb = 0;
	alamb = sqrt(x * y) + sqrt(x * z) + sqrt(y * z);
	x = 0.25 * (x + alamb);
	y = 0.25 * (y + alamb);
	z = 0.25 * (z + alamb);
	ave = (x + y + z) / 3.0;
	dX = (ave - x) / ave;
	dY = (ave - y) / ave;
	dZ = (ave - z) / ave;
}

double CarlsonFunction::rfFunction(double x, double y, double z)
{
	double dX, dY, dZ, ave, result;
	double c2, c3;
	try
	{
		if (std::min(x, std::min(y, z)) < 0)
		{
			throw std::runtime_error("Error: One or more argument is lower than 0. rfFunction.");
		}
		else if (std::min(x + y, std::min(x + z, y + z)) < 1.5e-38)
		{
			throw std::runtime_error("Error: Sum of any pair ex (x + y) is to close to 0. rfFunction.");
		}
		else if (std::max(x, std::max(y, z)) > 3e37)
		{
			throw std::runtime_error("Error: One of the arguments is to large. rfFunction.");
		}
	}
	catch (const std::exception& e)
	{
		std::cerr << "Exception caught: " << e.what() << std::endl;
		return 0;
	}
	dX = dY = dZ = ave = result = 0;
	// First iteration
	iterationForRFFunction(x, y, z, dX, dY, dZ, ave);
	while (std::max(fabs(x), std::max(fabs(y), fabs(z))) > 0.08)
	{
		iterationForRFFunction(x, y, z, dX, dY, dZ, ave);
	}
	c2 = dX * dY - pow(dZ, 2);
	c3 = dX * dY * dZ;
	// Compute rcFunction
	result = (1.0 + (c2 / 24.0 - 0.1 - 3.0 * c3 / 44.0) * (c2 + c3 / 14.0)) / sqrt(ave);
	return result;
}

void CarlsonFunction::iterationForRJFunction(double& x, double& y, double& z, double& dX, double& dY, double& dZ, double& dP, double& ave, double& so, double& f, double& p)
{
	double alamb, alfa, beta = 0;
	alamb = sqrt(x * y) + sqrt(x * z) + sqrt(y * z);
	alfa = pow((p * (sqrt(x) + sqrt(y) + sqrt(z)) + sqrt(x) * sqrt(y) * sqrt(z)), 2);
	beta = p * pow(p * alamb, 2);
	so = so + f * rcFunction(alfa, beta);
	f = 0.25 * f;
	x = 0.25 * (x + alamb);
	y = 0.25 * (y + alamb);
	z = 0.25 * (z + alamb);
	p = 0.25 * (p + alamb);
	ave = 0.2 * (x + y + z + 2.0 * p);
	dX = (ave - x) / ave;
	dY = (ave - y) / ave;
	dZ = (ave - z) / ave;
	dP = (ave - p) / ave;
}

double CarlsonFunction::rjFunction(double x, double y, double z, double p)
{
	double newX, newY, newZ, dX, dY, dZ, dP, ave, so, f, newP, result;
	double a1, b1, rcx;
	double c2, c3, c4, c5, c6;
	try
	{
		if (std::min(x, std::min(y, z)) < 0)
		{
			throw std::runtime_error("Error: One or more argument is lower than 0. rjFunction.");
		}
		else if (std::min(x + y, std::min(x + z, std::min(y + z, fabs(p)))) < 2.5e-13)
		{
			throw std::runtime_error("Error: Sum of any pair ex (x + y) or abs(p) is to close to 0. rjFunction.");
		}
		else if (std::max(x, std::max(y, std::max(z, fabs(p)))) > 9e11)
		{
			throw std::runtime_error("Error: One of the arguments is to large. rjFunction.");
		}
	}
	catch (const std::exception& e)
	{
		std::cerr << "Exception caught: " << e.what() << std::endl;
		return 0;
	}
	dX = dY = dZ = dP = ave = so = result = 0;
	f = 1.0;
	if (p > 0.0)
	{
		newX = x;
		newY = y;
		newZ = z;
		newP = p;
	}
	else
	{
		newX = std::min(x, std::min(y, z));
		newZ = std::max(x, std::max(y, z));
		newY = x + y + z - newX - newZ;
		a1 = 1.0 / (newY - p);
		b1 = a1 * (newZ - newY) * (newY - newX);
		newP = newY + b1;
		rcx = rcFunction(newX * newZ / newY, p * newP / newY);
	}
	// First iteration
	iterationForRJFunction(newX, newY, newZ, dX, dY, dZ, dP, ave, so, f, newP);
	while (std::max(fabs(x), std::max(fabs(y), std::max(fabs(z), fabs(p)))) > 0.05)
	{
		iterationForRJFunction(newX, newY, newZ, dX, dY, dZ, dP, ave, so, f, newP);
	}
	c2 = dX * dY + dX * dZ + dY * dZ;
	c3 = dX * dY * dZ;
	c4 = pow(dP, 2);
	c5 = c2 - 3.0 * c4;
	c6 = c3 + 2.0 * dP * (c2 - c4);
	// Compute rcFunction
	result = 3.0 * dP / 26.0 - 3.0 / 11.0;
	result = (1.0 / 6.0 + dP * result) * c3;
	result = result - dP * c4 / 3.0;
	result = dP * c2 * (1.0 / 3.0 - 3.0 * dP / 22.0) + result;
	result = 1.0 + c5 * (9.0 * c5 - 4.5 * c6 / 26.0 - 3.0 / 14.0) + result;
	result = f * result;
	result = 3.0 * so + result / (ave * sqrt(ave));
	if (p <= 0.0)
	{
		result = a1 * (b1 * result + 3.0 * (rcx - rfFunction(newX, newY, newZ)));
	}
	return result;
}
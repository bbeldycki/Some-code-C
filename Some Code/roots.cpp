#include "roots.h"

Roots::Roots()
{
	
}

Roots::~Roots()
{

}

void Roots::laguer(std::vector<std::complex<double>>& a, std::complex<double>& x, int& its)
{
	const int mr = 8, mt = 10, maxit = mt * mr;
	const double esp = 1.e-15;
	const std::vector<double> frac = { 0.0, 0.5, 0.25, 0.75, 0.13, 0.38, 0.62, 0.88, 1.0 };
	std::complex<double> dx, xOn, b, d, f, g, h, sq, gp, gm;
	double err;
	int m = a.size() - 1;
	for (int i = 1; i <= maxit; i++)
	{
		its = i;
		b = a[m];
		d = f = 0.0;
		err = std::abs(b);
		for (int j = m - 1; j >= 0; j--)
		{
			f = x * f + d;
			d = x * d + b;
			b = x * b + a[j];
			err = std::abs(b) + std::abs(x) * err;
		}
		err *= esp;
		if (std::abs(b) <= err)
		{
			return;
		}
		g = b / d;
		h = std::pow(g, 2) - 2.0 * f / b;
		sq = std::sqrt(double(m - 1) * (double(m) * h - std::pow(g, 2)));
		gp = g + sq;
		gm = g - sq;
		if (std::abs(gp) < std::abs(gm))
		{
			gp = gm;
		}
	}
}
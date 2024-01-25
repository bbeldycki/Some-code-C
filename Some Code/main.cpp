#include "integralshelper.h"
#include "roots.h"

int main(int argc, char* argv[])
{
	std::vector<int> exOne = {1, 1, 1};
	std::vector<double> exTwo = { 1.f, 1.f, 1.f };
	std::vector<double> exThree = { 2.f, 1.f, 4.f };
	std::complex<double> conv;
	int ajaj = 1;
	double x = 7., y = 6., z = 1.;
	//double y = 6.;
	//double z = 1.;
	IntegralsHelper exercise;
	Roots root;
	z = exercise.computeAPOnePn(exOne, exTwo, exThree, y, x);
	std::cout << z << std::endl;
	// need to convert vector<double> -> vector<complex<double>>
	root.laguer(exTwo, x, ajaj);
	return 0;
}
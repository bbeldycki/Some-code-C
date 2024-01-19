#include "integralshelper.h"

int main(int argc, char* argv[])
{
	std::vector<int> exOne = {1, 1, 1};
	std::vector<double> exTwo = { 1.f, 1.f, 1.f };
	std::vector<double> exThree = { 2.f, 1.f, 4.f };
	double x = 7;
	double y = 6;
	double z = 1;
	std::unique_ptr<IntegralsHelper> exercise = std::make_unique<IntegralsHelper>();
	z = exercise->computeAPOnePn(exOne, exTwo, exThree, y, x);
	std::cout << z << std::endl;
	return 0;
}
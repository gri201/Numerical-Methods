#include <iostream>
#include <cmath>

using namespace std;
const double E = 2.71828;

double f(double x)
{
	if (x == 0) return 0;
	else if (x == 1) return 3 * pow(sin(1), 3) / (E * log(2));
	else return ((pow(x,3) - 1) * pow(sin(x), 3) * exp(-pow(x,3))) / (sqrt(x) * log(x) * log(1+x));
}

double Integrate(double a, double b, int N)
{
	double len = (b - a) / N;
	double S = 0;

	for(int i = 1; i <= N; ++i)
	{
		S += len * f(a + len * (i - 1/2));
	}
	return S;
}

double Runge(double a, double b, int &N, double eps)
{
	N = 1; 
	double S_N = pow(10, 10);
	double S_2N = (b - a) * f(b);
	double len = (b - a);

	while (fabs(S_2N - S_N) / 3 >= eps)
	{
		//cout << fabs(S_2N - S_N) / 3 << endl;
		N *= 2;
		S_N = S_2N;
		S_2N /= 2;
		len /= 2;

		for(int i = 1; i <= N; i += 2)
		{
			S_2N += len * f(a + len * i);
		}
	}

	return S_2N;
}


int main()
{
	double eps = 0.002;
	double delta_1 = 0.146;
	double delta_2 = 0.0005;
	double C = 1.884;
	double a_c = 0;
	double b_c = 0;

	int N_1 = 2827567;
	int N_2 = 1389714;
	int N_1R = 0;
	int N_2R = 0;
	int N_c = 0;
	
	double S_1 = Integrate(delta_1, 1 - delta_2, N_1);
	double S_2 = Integrate(1 + delta_2, C, N_2);
	double S_1R = Runge(delta_1, 1 - delta_2, N_1R, eps);
	double S_2R = Runge(1 + delta_2, C, N_2R, eps);

	cout << fixed;
	cout.precision(5);

	cout << "Results for usual calculation:" << endl;
	cout << "Integral of f(x) in [0.146, 0.9995]: " << S_1 << endl;
	cout << "Integral of f(x) in [1.0005, 1.884]: " << S_2 << endl;
	cout << "Total: " << S_1 + S_2 << endl;
	cout << endl;
	cout << "Results fot Runge calculation:" << endl;
	cout << "Integral of f(x) in [0.146, 0.9995]: " << S_1R << " N_1 = " << N_1R << endl;
	cout << "Integral of f(x) in [1.0005, 1.884]: " << S_2R << " N_2 = " << N_2R << endl;
	cout << "Total: " << S_1R + S_2R << " N = " << N_1R + N_2R << endl;	

	cout << endl;

	cout << "Custom request" << endl;
	cout << "Input a: ";
	cin >> a_c;
	cout << endl;
	cout << "Input b: ";
	cin >> b_c;
	cout << endl;
	cin >> N_c;
	cout << endl;
	cout << "Integral f(x) from " << a_c << " to " << b_c << " is: " << Integrate(a_c, b_c, N_c) << endl;;

	

}
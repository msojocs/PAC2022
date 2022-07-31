#include <iostream>
#include <complex>
#include <cmath>

using namespace std;
int main(void)
{
    cout << "start" << endl;
    complex<double> a(1, 2);
    a.real(3);
    cout << "a: " << a.real() << ", b:" << a.imag() << endl;
    return 0;
}
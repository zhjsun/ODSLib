#include <iostream>

#include <Eigen/Dense>
#include "ODSLib/ODSLib.h"

using namespace std;
using namespace Eigen;
using namespace ODS;

int main( void ) {

    // Test Cartesian-Element transformation
    VectorXd state(6);
    ArrayXd Elem(6);

    state << 6778137, 0, 0, 1, 7000, 1;
    Cart2Elem(state, Elem);
    Elem2Cart(Elem, state);
    cout << state << endl;

    // Test ICS2VVLH transformation
    VectorXd tar(6), cha(6), rel(6);
    tar << 6778137, 0, 0, 1, 7000, 1;
    cha << 6778037, 100, 100, 1, 7000, 1;
    ICS2VVLH(tar, cha, rel);
    cout << rel << endl;

    return 0;
}


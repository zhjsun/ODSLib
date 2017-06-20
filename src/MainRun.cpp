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

    return 0;
}


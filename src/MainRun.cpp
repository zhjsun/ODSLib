#include <Eigen/Dense>
#include "ODSLib.h"

using namespace std;
using namespace Eigen;
using namespace ODS;

int main() {

    VectorXd state(6);
    ArrayXd Elem(6);

    state << 6778137, 0, 0, 0, 7000, 0;
    Cart2Elem(state, Elem);
    cout << Elem << endl;

    return 0;
}


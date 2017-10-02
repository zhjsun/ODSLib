#include <iostream>

#include <Eigen/Dense>
#include "ODSLib/ODSLib.h"

using namespace std;
using namespace Eigen;
using namespace ODS;

int main( void ) {

    // Test coordinates transformation
    VectorXd tar(6), cha(6), rel(6);
    ArrayXd Elem(6);
    tar << 6778137, 1000, 000, 1, 7000, 0;
    cout << "Target state:" << endl << tar << endl;

    Cart2Elem(tar, Elem);
    cout << "Target element:" << endl << Elem << endl;

    Elem << 6778137, 0.999, 0.0, atan(1.0)*2, atan(1.0)*2, atan(1.0)*2;
    Elem2Cart(Elem, tar);
    cout << "Target state:" << endl << tar << endl;

    /*
    rel << -100, -1000, 50, 1, 1, 1;
    cout << "Chaser in LVLH: " << endl << rel << endl;

    LVLH2ICS(tar, rel, cha);
    cout << "Chaser in ICS:" << endl << cha << endl;

    ICS2VVLH(tar, cha, rel);
    cout << "Chaser in VVLH: " << endl << rel << endl;

    VVLH2ICS(tar, rel, cha);
    cout << "Chaser in ICS:" << endl << cha << endl;

    ICS2LVLH(tar, cha, rel);
    cout << "Chaser in LVLH: " << endl << rel << endl;
    */

    return 0;
}


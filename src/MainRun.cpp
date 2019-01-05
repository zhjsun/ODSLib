#include <iostream>

#include <eigen3/Eigen/Dense>
#include <dace/dace.h>
#include "ODSLib/ODSLib.h"

using namespace std;
using namespace Eigen;
using namespace DACE;

int main( void ) {

    // Test coordinates transformation
    VectorXd tar(6), cha(6), rel(6);
    VectorXd Elem(6);

    // tar << 6778137, 1000, 2000, 1, 7000, 1;
    // Elem << 6778137, 0.0, 0.1, 2.0*atan(1.0), 2*atan(1), 0;
    // cout << "Target state:" << endl << tar << endl;

    Elem << 6778137, 0.1, 0.2, ODS::Pi/2, ODS::Pi, atan(1.0)*2;
    cout << "Original element:" << endl << Elem << endl;

    ODS::Elem2Cart(Elem, tar);
    cout << "Target state:" << endl << tar << endl;

    ODS::Cart2Elem(tar, Elem);
    cout << "Target element:" << endl << Elem << endl;

    DA::init(1, 6);
    Matrix<DA, -1, 1> DA_elem(6), DA_state(6);
    DA_elem(0) = 6778137.0 + DA(1);
    DA_elem(1) = 0.1       + DA(2);
    DA_elem(2) = 0.2       + DA(3);
    DA_elem(3) = ODS::Pi/2 + DA(4);
    DA_elem(4) = ODS::Pi/2 + DA(5);
    DA_elem(5) = 1.0*ODS::Pi/2 + DA(6);
    cout << "Orbital element:" << endl;
    for(int i = 0; i < 6; ++i) {
        cout << DA_elem(i) << endl;
    }

    ODS::Elem2Cart(DA_elem, DA_state);
    cout << "Orbital state:" << endl;
    for(int i = 0; i < 6; ++i) {
        cout << DA_state(i) << endl;
    }

    ODS::Cart2Elem(DA_state, DA_elem);
    cout << "Orbital element:" << endl;
    for(int i = 0; i < 6; ++i) {
        cout << DA_elem(i) << endl;
    }

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


#include <iostream>
#include "abstractbasis.h"
#include "jelliumbasis.h"

using namespace std;

int main()
{

    int NumberOfShellsStochastic = 2;
    int NumberOfElectrons = 2;
    double r_s = 0.5;

    /* Single particle basis
     * for the jellium model using std::vector/cartesian classes
     * need to be rewritten via 1-D arrays
     */
    abstractbasis * SPbasis = new jelliumbasis(NumberOfShellsStochastic, NumberOfElectrons, r_s);
    SPbasis->getQuantumStates();

    delete SPbasis;
    return 0;
}

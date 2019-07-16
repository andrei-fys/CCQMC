#include <iostream>
#include "abstractbasis.h"
#include "jelliumbasis.h"

using namespace std;

int main()
{

    int NumberOfShellsStochastic = 2;
    int NumberOfElectrons = 2;
    double r_s = 0.5;

    abstractbasis * SPbasis = new jelliumbasis(NumberOfShellsStochastic, NumberOfElectrons, r_s);
    SPbasis->getQuantumStates();

    return 0;
}

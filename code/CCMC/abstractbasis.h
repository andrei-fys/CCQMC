#ifndef ABSTRACTBASIS_H
#define ABSTRACTBASIS_H

#include "cartesianstate.h"
#include <vector>
#include <armadillo>
#include <iostream>

class abstractbasis
{
public:
    abstractbasis(){}
    ~abstractbasis(){}

    virtual void getQuantumStates() = 0;
    virtual void getQuantumStatesNumber() = 0;
    virtual double TBME(int, int, int, int, const std::vector<CartesianState> &shells) = 0;
    virtual arma::vec computeSPenergies(const std::vector<CartesianState> &shells) = 0;


    virtual std::vector<CartesianState> getStateVec () = 0;
    virtual int getShellsExact      () = 0;
    virtual int getShellsStochastic () = 0;
    virtual int getFermiLevel       () = 0;
    virtual int getStatesExact      () = 0;
    virtual int getStatesStochastic () = 0;
    virtual int getnMax () = 0;
    virtual CartesianState oneState(int) = 0;
    virtual CartesianState sumState(int, int) = 0;
    virtual CartesianState substractState(int, int) = 0;
    virtual CartesianState sumSubstractState(int, int, int) = 0;
    virtual bool isEqual(CartesianState, CartesianState) = 0;
    virtual double calculateReferenceEnergy() = 0;
    virtual arma::vec computeSPenergiesCCMC() = 0;




};

#endif // ABSTRACTBASIS_H

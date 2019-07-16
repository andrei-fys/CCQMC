#ifndef JELLIUMBASIS_H
#define JELLIUMBASIS_H
#include "abstractbasis.h"

class jelliumbasis : public abstractbasis
{
public:

    //constr
    jelliumbasis(int , int , double );
    ~jelliumbasis();
    //vars
    int m_ShellsExact, m_ShellsStochastic, m_FermiLevel, m_StatesExact, m_StatesStochastic, m_nMax;
    double m_L, m_L2 , m_L3;


    //methods
    double TBME(int, int, int, int, const std::vector<CartesianState> &shells);
    arma::vec computeSPenergies(const std::vector<CartesianState> &shells);
    virtual arma::vec computeSPenergiesCCMC();
    void getQuantumStates();
    void getQuantumStatesNumber();
    virtual std::vector<CartesianState> getStateVec () {return this->m_shells;}

    virtual int getShellsExact      () {return this->m_ShellsExact;}
    virtual int getShellsStochastic () {return this->m_ShellsStochastic;}
    virtual int getFermiLevel       () {return this->m_FermiLevel;}
    virtual int getStatesExact      () {return this->m_StatesExact;}
    virtual int getStatesStochastic () {return this->m_StatesStochastic;}
    virtual int getnMax             () {return this->m_nMax;}
    virtual double calculateReferenceEnergy();


    virtual CartesianState oneState(int);
    virtual CartesianState sumState(int, int);
    virtual CartesianState substractState(int, int);
    virtual CartesianState sumSubstractState(int, int, int);
    virtual bool isEqual(CartesianState, CartesianState);


private:
    // vars
    const double pi = M_PI;
    const double eVsHartree = 27.21138505; // eV
    const double FineStruc = 1.0/137.035999139; // a.u.
    double Massc2 = 0.5109989461; // MeV electron mass
    double hbarc = 0.1973269788; // eV micron
    double rBohr;
    double m_prefactor;
    //double m_L;
    double Esq;


    // methods
    void setUpStatesCartesian();
    double VectorAbsoluteDifferenceSquared(arma::vec, arma::vec);
    int KroneckerDelta(arma::vec, arma::vec);
    void setUpUnits(int, double);
    void computeSPenergiesHF(const std::vector<CartesianState> &, arma::vec &);


protected:
    std::vector<CartesianState> m_shells;

};

#endif // JELLIUMBASIS_H

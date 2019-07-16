#include "jelliumbasis.h"

jelliumbasis::jelliumbasis(int NumberOfShellsStochastic, int ParticlesNumber, double r_s){

    setUpUnits(ParticlesNumber, r_s);
    //this->m_L3 = 4.0*pi*ParticlesNumber*r_s*r_s*r_s/3.0;
    //this->m_L = pow(m_L3, 1.0/3.0);
    this->m_L2 = m_L*m_L;
    this->m_L3 = m_L*m_L*m_L;
    this->m_FermiLevel = ParticlesNumber;
    //this->m_ShellsExact = NumberOfShellsExact;
    this->m_ShellsStochastic = NumberOfShellsStochastic;
    setUpStatesCartesian();
}

void jelliumbasis::setUpUnits(int FermiLevel, double r_s){

    this->hbarc = hbarc*10000/eVsHartree; // Hartree Angstrom
    this->Massc2 = Massc2*1000000/eVsHartree; // Hartree
    this->rBohr = hbarc/(Massc2*FineStruc); // Angstrom
    double R0 = r_s*rBohr;
    this->m_prefactor = hbarc*hbarc/(2.*Massc2);
    this->m_L = pow( (4./3.)*pi*FermiLevel, 1./3.)*R0; // Angstrom
    this->Esq = FineStruc*hbarc; // Hartree Angstrom

}

void jelliumbasis::setUpStatesCartesian() {
    // set ups states as vector of statevectors
    // statevector has folloving format (nx, ny, nz, s)

    int s = 1;
    for (int E = 0; E < this->m_ShellsStochastic; E++){
        int nMax = sqrt(E);
        for (int nx = -nMax ; nx <= nMax+1; nx++){
            for (int ny = -nMax; ny <= nMax+1; ny++){
                for (int nz = -nMax; nz <= nMax+1; nz++){
                    if (nx*nx + ny*ny + nz*nz == E){
                        m_shells.emplace_back(CartesianState());
                        m_shells.back().set(nx, ny, nz, s);
                        m_shells.emplace_back(CartesianState());
                        m_shells.back().set(nx, ny, nz, -1*s);
                    }
                }
            }
        }
    }
    m_nMax = sqrt(this->m_ShellsStochastic);
    this->m_StatesStochastic = m_shells.size();
    //number of states for _total_ number of shells
    int total = 0;
    //for (int E = 0; E < this->m_ShellsExact; E++){
    for (int E = 0; E < this->m_FermiLevel; E++){
        int nMax = sqrt(E);
        for (int nx = -nMax ; nx <= nMax+1; nx++){
            for (int ny = -nMax; ny <= nMax+1; ny++){
                for (int nz = -nMax; nz <= nMax+1; nz++){
                    if (nx*nx + ny*ny + nz*nz == E){
                        total++;
                    }
                }
            }
        }
    }
    this->m_StatesExact = 2*total;
}

double jelliumbasis::VectorAbsoluteDifferenceSquared(arma::vec A, arma::vec B){
    double C = 0;
    for (int i = 0; i < 3; i++){
        C += (A(i) - B(i))*(A(i) - B(i));
    }
    return C;
}

int jelliumbasis::KroneckerDelta(arma::vec A, arma::vec B){
    int KroneckerDelta = 1;
    for(int i = 0; i < 3; i++){
        KroneckerDelta *= (A(i)==B(i)) ? 1 : 0;
    }
    return KroneckerDelta;
}

void jelliumbasis::getQuantumStates(){
    for(CartesianState quantum_state : m_shells){
        std::cout << "E = " << quantum_state.E() << std::endl;
        std::cout << "nx = " << quantum_state.nx() << std::endl;
        std::cout << "ny = " <<quantum_state.ny() << std::endl;
        std::cout << "nz = " <<quantum_state.nz() << std::endl;
        std::cout << "s = " <<quantum_state.s() << std::endl;
        std::cout << "-------------------------" << std::endl;
    }
}

void jelliumbasis::getQuantumStatesNumber(){
    std::cout << "System : Electon gas " << std::endl;
    std::cout << "Fermi level is " << this->m_FermiLevel << std::endl;
    std::cout << " --- *** ---- " << std::endl;
    std::cout << "Number of exact states is " <<this->m_StatesExact << std::endl;
    std::cout << " --- *** ---- " << std::endl;
    std::cout << "Number of total shells is " << this->m_ShellsStochastic << std::endl;
    std::cout << "Number of total states is " << this->m_StatesStochastic << std::endl;
}


double jelliumbasis::TBME(int p, int q, int r, int s, const std::vector<CartesianState> &shells){

    CartesianState quantum_state_p = shells.at(p);
    double kp_x = quantum_state_p.nx();
    double kp_y = quantum_state_p.ny();
    double kp_z = quantum_state_p.nz();
    int p_s = quantum_state_p.s();

    CartesianState quantum_state_q = shells.at(q);
    double kq_x = quantum_state_q.nx();
    double kq_y = quantum_state_q.ny();
    double kq_z = quantum_state_q.nz();
    int q_s = quantum_state_q.s();

    CartesianState quantum_state_r = shells.at(r);
    double kr_x = quantum_state_r.nx();
    double kr_y = quantum_state_r.ny();
    double kr_z = quantum_state_r.nz();
    int r_s = quantum_state_r.s();

    CartesianState quantum_state_s = shells.at(s);
    double ks_x = quantum_state_s.nx();
    double ks_y = quantum_state_s.ny();
    double ks_z = quantum_state_s.nz();
    int s_s = quantum_state_s.s();

    arma::vec k_p;
    arma::vec k_q;
    arma::vec k_r;
    arma::vec k_s;

    k_p << kp_x << kp_y << kp_z;
    k_q << kq_x << kq_y << kq_z;
    k_r << kr_x << kr_y << kr_z;
    k_s << ks_x << ks_y << ks_z;

    int KD_kpkq_krks = KroneckerDelta(k_p + k_q, k_r + k_s);
    double first_term = 0.0;
    double second_term = 0.0;

    if (KD_kpkq_krks == 0){
       return 0.0;
    } else {
      if ((p_s == r_s) && (q_s == s_s)){
          if  (KroneckerDelta(k_p,k_r) == 0){
              first_term = m_L2/(4.0*pi*pi*VectorAbsoluteDifferenceSquared(k_r, k_p));
          }

      }

      if ((p_s == s_s) && (q_s == r_s) ){
          if (KroneckerDelta(k_p,k_s) == 0) {
              second_term = m_L2/(4.0*pi*pi*VectorAbsoluteDifferenceSquared(k_s, k_p));
          }
      }
      return (Esq*4.0*pi/m_L3)*(first_term - second_term);
    }
}

arma::vec jelliumbasis::computeSPenergies(const std::vector<CartesianState> &shells){
    int NumberOfStates = shells.size();
    arma::vec Epsilon;
    Epsilon.zeros(NumberOfStates);
    for(int x = 0; x < NumberOfStates; x++){
        CartesianState quantum_state_x = shells.at(x);
        Epsilon(x) = quantum_state_x.E()*(4.0*pi*pi*m_prefactor/m_L2);
    }
    computeSPenergiesHF(shells, Epsilon);
    return Epsilon;
}

void jelliumbasis::computeSPenergiesHF(const std::vector<CartesianState> &shells, arma::vec &Epsilon){
    double Ei,Ea;
    int NumberOfStates = shells.size();
    for(int i = 0; i < m_FermiLevel; i++) {
        Ei = Epsilon(i);
        for(int j = 0; j < m_FermiLevel; j++) {
            Ei += TBME(i,j,i,j, shells);
        }
        Epsilon(i) = Ei;
    }
    for(int a = m_FermiLevel; a < NumberOfStates; a++) {
        Ea = Epsilon(a);
        for(int j = 0 ; j < m_FermiLevel; j++) {
            Ea += TBME(a,j,a,j, shells);
        }
        Epsilon(a) = Ea;
    }
}

CartesianState jelliumbasis::oneState(int p){
    int Nx = getStateVec().at(p).nx();
    int Ny = getStateVec().at(p).ny();
    int Nz = getStateVec().at(p).nz();
    int S =  getStateVec().at(p).s();
    CartesianState QuantumState = CartesianState();
    QuantumState.set(Nx, Ny, Nz, S);
    return QuantumState;
}

CartesianState jelliumbasis::sumState(int p, int q){
    int Nx = getStateVec().at(p).nx() + getStateVec().at(q).nx();
    int Ny = getStateVec().at(p).ny() + getStateVec().at(q).ny();
    int Nz = getStateVec().at(p).nz() + getStateVec().at(q).nz();
    int S = getStateVec().at(p).s() + getStateVec().at(q).s();
    CartesianState QuantumState = CartesianState();
    QuantumState.set(Nx, Ny, Nz, S);
    return QuantumState;
}

CartesianState jelliumbasis::substractState(int p, int q){
    int Nx = getStateVec().at(p).nx() - getStateVec().at(q).nx();
    int Ny = getStateVec().at(p).ny() - getStateVec().at(q).ny();
    int Nz = getStateVec().at(p).nz() - getStateVec().at(q).nz();
    int S = getStateVec().at(p).s() - getStateVec().at(q).s();
    CartesianState QuantumState = CartesianState();
    QuantumState.set(Nx, Ny, Nz, S);
    return QuantumState;
}

CartesianState jelliumbasis::sumSubstractState(int p, int q, int r){
    int Nx = getStateVec().at(p).nx() + getStateVec().at(q).nx() - getStateVec().at(r).nx();
    int Ny = getStateVec().at(p).ny() + getStateVec().at(q).ny() - getStateVec().at(r).ny();
    int Nz = getStateVec().at(p).nz() + getStateVec().at(q).nz() - getStateVec().at(r).nz();
    int S = getStateVec().at(p).s() + getStateVec().at(q).s() - getStateVec().at(r).s();
    CartesianState QuantumState =  CartesianState();
    QuantumState.set(Nx, Ny, Nz, S);
    return QuantumState;
}

bool jelliumbasis::isEqual(CartesianState state1, CartesianState state2){
    if (   state1.nx() == state2.nx()
        && state1.ny() == state2.ny()
        && state1.nz() == state2.nz()
        && state1.s()  == state2.s()
            ) {
        return true;
    } else {
        return false;
    }
}

double jelliumbasis::calculateReferenceEnergy(){
    //returns the reference energy in the current basis
    double reference_energy = 0.0;

    int NumberOfStates = m_shells.size();
    arma::vec Epsilon;
    Epsilon.zeros(NumberOfStates);

    for(int x = 0; x < NumberOfStates; x++){
        CartesianState quantum_state_x = m_shells.at(x);
        Epsilon(x) = quantum_state_x.E()*(4.0*pi*pi*m_prefactor/m_L2);
    }

    for(int i=0; i < getFermiLevel(); i++){
        reference_energy += Epsilon(i);
        for(int j=0; j < getFermiLevel(); j++){
            if(i!=j){
                reference_energy += 0.5*TBME(i,j, i,j, m_shells);
                //cout << calculateTwoBodyElements(i,j, i,j) << endl;
            }
        }
    }
    return reference_energy;
}


arma::vec jelliumbasis::computeSPenergiesCCMC(){

    int NumberOfStates = m_shells.size();
    arma::vec Epsilon;
    Epsilon.zeros(NumberOfStates);

    for(int x = 0; x < NumberOfStates; x++){
        CartesianState quantum_state_x = m_shells.at(x);
        Epsilon(x) = quantum_state_x.E()*(4.0*pi*pi*m_prefactor/m_L2);
    }
    return Epsilon;
}

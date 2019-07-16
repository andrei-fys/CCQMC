#include "cartesianstate.h"

CartesianState::CartesianState()
{

}

void CartesianState::set(int nx, int ny, int nz, int s){
    m_nx = nx;
    m_ny = ny;
    m_nz = nz;
    m_s = s;

}

void CartesianState::flipSpin(){
    m_s = m_s*(-1);
}

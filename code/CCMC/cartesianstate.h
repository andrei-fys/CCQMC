#ifndef CARTESIANSTATE_H
#define CARTESIANSTATE_H


class CartesianState
{
private:
    int m_nx;     // nx
    int m_ny;     // ny
    int m_nz;     // nz
    int m_s;      // spin


public:
    CartesianState();
    double s() { return m_s; }
    int nx() { return m_nx; }
    int ny() { return m_ny; }
    int nz() { return m_nz; }
    double E() {return (m_nx*m_nx + m_ny*m_ny + m_nz*m_nz);} //*0.5
    void set(int, int, int, int);
    void flipSpin();
};

#endif // CARTESIANSTATE_H

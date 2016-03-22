#include <arrayfire.h>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
// armazon de la clase encargada de las simulaciones
using namespace af;
class LatticeBolztmannD3Q15{
 private:
  int q;
  double beta;
  double UMAX;
  int Lx,Ly,Lz; // variables del lattice
  float V[3][15]; // velocidades del lattice
  array fs; // un arreglo de 4 dimensiones  
 
  array w; // funciones de peso
  array vel_x;
  array vel_y;
  array vel_z;
  // condiciones
  array ConditionsU;
 public:
  
  LatticeBolztmannD3Q15(int Lx, int Ly, int Lz); // constructor - reciba el tamaño del lattice
  void Inicie(float r0, float Ux0, float Uy0, float Uz0); // reciba la densidad y velocidades cero
  
  array feq(int i, array &rhos, array &Uxs , array &Uys , array &Uzs);
  array feq2( array &rhos, array &Uxs , array &Uys , array &Uzs);
  array feq3( array &rhos ,array &Us );
  array aux( array &vels, array &Ua );
  void SetConditions( array &Ux , array &Uy, array &Uz);
  array rho();
  array Jx();
  array Jy();
  array Jz();
  
  void Colission(void);
  void Adveccion(void);
  void Run(int steps);
  void Eval();
  
  void SaveFs(std::string route_f );
  void setArrayC(array c);
  void saveVTK(std::string route);
};



// Intento  Modelo 3d para lattice bolztmann utlizando el metodo entropico
/*  Andres Huertas Suarez 2016 
 */ 
// 
#include <arrayfire.h>
#include <LBD3Q15.h> // just a header holding the class definitions 
using namespace af;
LatticeBolztmannD3Q15::LatticeBolztmannD3Q15(int lx,int ly, int lz){
  //constructor
  //this->q = 15;
  this->q = 15;
  this->Lx=lx; // tamaño del lattice
  this->Ly=ly;
  this->Lz=lz;
  
  fs = af::randu( q , Lx , Ly , Lz); // the distribution funtions are 15 for each space site.
  
  w = af::constant( 0 , 15 , f32);
  //
  // the model work with 15 discrete velocity vector as described here
  int vx[] = {0,1,0,0,-1,0,0,1,-1,1,1,-1,1,-1,-1  }; 
  int vy[] = {0,0,1,0,0,-1,0,1,1,-1,1,-1,-1,1,-1  };
  int vz[] = {0,0,0,1,0,0,-1,1,1,1,-1,1,-1,-1,-1 };
  this->ConditionsU = af::constant(1, Lx,Ly,Lz,f32);
  
  // create a set of 3 ArrayFire arrays holding the component of this vectors. 

  vel_x = array( q , vx);
  vel_y = array (q , vy);
  vel_z = array(q , vz);
 
  // the method also needs 15 weigths, one for each velocity, those are the next:
  w(0) = 16; w(1)=w(2)=w(3)=w(4)=w(5)=w(6) = 8;
  w(7)=w(8)=w(9)=w(10)=w(11)=w(12)=w(13)=w(14)=1;
  w = 1/72.0*w;
  
  // this vectors are the same as vel_x vel_y and vel_z,  i keep this repeated information, because the streming step  uses int to shift the arrays in each of the 15 directions
  V[0][0] = 0; V[0][1] = 1; V[0][2] = 0; V[0][3] = 0 ; V[0][4] = -1;
  V[0][5] = 0; V[0][6] = 0; V[0][7] = 1; V[0][8] = -1 ; V[0][9] = 1;
  V[0][10] = 1; V[0][11] = -1; V[0][12] = 1; V[0][13] = -1 ; V[0][14] = -1;
  
  V[1][0] = 0; V[1][1] = 0; V[1][2] = 1; V[1][3] = 0 ; V[1][4] = 0;
  V[1][5] = -1; V[1][6] = 0; V[1][7] = 1; V[1][8] = 1 ; V[1][9] = -1;
  V[1][10] = 1; V[1][11] = -1; V[1][12] = -1; V[1][13] = 1 ; V[1][14] = -1;

  V[2][0] = 0; V[2][1] = 0; V[2][2] = 0; V[2][3] = 1 ; V[2][4] = 0;
  V[2][5] = 0; V[2][6] = -1; V[2][7] = 1; V[2][8] = 1 ; V[2][9] = 1;
  V[2][10] = -1; V[2][11] = 1; V[2][12] = -1; V[2][13] = -1 ; V[2][14] = -1;

 
}
void LatticeBolztmannD3Q15::Inicie(float r0 , float Ux0 , float Uy0 , float Uz0){
  //here i just starting things up. 
  array rh0 = af::constant( r0 ,1 ,Lx, Ly, Lz , f32);
  array Uxs = af::constant( Ux0,1 ,Lx , Ly , Lz , f32);
  array Uys = af::constant( Uy0,1, Lx , Ly , Lz , f32);
  array Uzs = af::constant( Uz0,1, Lx , Ly , Lz , f32);
  fs = feq( rh0 , Uxs , Uys , Uzs);

}
void LatticeBolztmannD3Q15::SetConditions(array &Ux, array &Uy , array &Uz){
  
  // LID DRIVEN CAVITY  
  // this set the physical conditions over the velocity field so simulate a lid driven cavity.
  // in theory this is call each iteration step. 

  Uy(span , 0 , span) = 0;

  // plano y = Ly , velocidad
  
  Uy(span , Ly-1 , span) = 0;

  // plano x = 0

  Ux(0,span,span ) = 0.0;
  Uy(0,span,span ) = 0;
  Uz(0,span,span ) = 0;
  

   // plano x = Lx
  Ux(Lx-1,span,span ) = 0.0;
  Uy(Lx-1,span,span ) = 0;
  Uz(Lx-1,span,span ) = 0;
  
  // plano z = 0
  Ux(span,span,0 ) = 0.0;
  Uy(span,span, 0) = 0;
  Uz(span,span,0 ) = 0;
 
   // plano z = Lz
  Ux(span,span,Lz-1 ) = 0.03;
  Uy(span,span, Lz-1) = 0;
  Uz(span,span, Lz-1) = 0;
  
}
array LatticeBolztmannD3Q15::feq2( array &rhos, array &Uxs , array &Uys , array &Uzs){ 
  
  // this is tricky to explain, but you can assume this works, or that this is what i need to calculate
  // for those familarized with the feq in the entropic lattice bolztmann,  you may recall that if holds a productoria , this term is calculates with the aux function
  return (  af::tile( rhos ,  q )*af::tile( w , 1 , Lx , Ly , Lz ))*aux(vel_x, Uxs)*aux(vel_y,Uys)*aux(vel_z,Uzs);
  
}

array  LatticeBolztmannD3Q15::aux(array &vels, array &Ua){
  // aux function, tricky to explain but this is what i need
  array f2v = af::tile( (2*Ua+af::sqrt(1+3*(Ua*Ua) ))/(1-Ua) , q );
  array f3v = af::tile( vels , 1 , Lx , Ly , Lz  );
  array res = af::pow( f2v ,  f3v   ); 
  return  res;

}


void LatticeBolztmannD3Q15::Adveccion(){
  // each set of distributions,  has to be shifted to the next positions, acording to the velocity 
  // vectors. for example all the distribution functions corresponding to f(1,span,span,span) has to migrate one position in the x axis.  
  for( int i = 1 ; i < q ; i++) { 
    fs( i , span,span,span) = af::shift( fs(i ,span,span,span) , 0 , V[0][i] , V[1][i] , V[2][i] );
  }
  
}

void LatticeBolztmannD3Q15::Run(int steps){
  
  for( int i = 0 ; i < steps; i++){
    
    std::cout << "Step:" << i << std::endl;
    // here i calculate the macroscopic cuantities, needed to calculate the feq
    // this is calculated acording to the Lattice-BOlztmann method 
    array rhos = af::sum( fs , 0 );
    array jx =  af::sum( fs*af::tile( vel_x , 1 , Lx , Ly , Lz ) , 0) ; 
    array jy =  af::sum( fs*af::tile( vel_y , 1 , Lx , Ly , Lz ) ,0 );
    array jz =  af::sum( fs*af::tile( vel_z , 1 , Lx, Ly ,  Lz ) ,0);
    
    array vs = af::join( 0 , jx , jy , jz )/af::tile( rhos , 3) ; 
    
    //af_print ( vs ); 
    array ux = vs(0,span,span,span); // x component of the velocity
    array uy = vs(1,span,span); // y component of the velocity
    array uz = vs(2,span,span); // z component of the velocity
    array feq =  feq2( rhos ,  ux , uy , uz ); // feq.
    // relaxation step.
    fs = fs + beta*(feq - fs);
    
    Adveccion(); 
    
  }
  return ;
}

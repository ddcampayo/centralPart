//#include"pParticles.h"
#include"linear.h"
#include"fields_enum.h"
#include"simu.h"


// Solve for pressure
// The famous pressure Poisson equation

void linear::s_equation(const FT dt ) {

  cout << "Solving s equation " << endl;

  //  fill_Delta_DD(); // This may be important -- or not

  FT ddt = dt;
  if( dt < 1e-10 ) ddt = 1;  // for debugging, mainly

  //  volumes( T );
  
  VectorXd I0  = field_to_vctr( sfield_list::I0 ) ;
  VectorXd I   = field_to_vctr( sfield_list::I ) ;

  VectorXd DI = I - I0 ;

  FT DI_sigma =  DI.array().square().sum() ;
  FT I_mean  =  I.array().square().sum() ;

  cout << " s field  "
       << " rel DI std dev: " << sqrt( DI_sigma / I_mean )
       << endl;

  VectorXd Ds  =  MM_solver.solve( DI );

  VectorXd s0  = field_to_vctr( sfield_list::s ) ;

  vctr_to_field( s0 + Ds / ( ddt * ddt) , sfield_list::s  ) ;

  return;
}




void linear::u_add_grad_s( const FT dt ) {

  VectorXd vol  = field_to_vctr( sfield_list::vol );
  // perhaps mean vol would be just fine

  VectorXd grad_sx,grad_sy;

  MM_times_sfield( sfield_list::s  ,  grad_sx, grad_sy);
  
  VectorXd Ustar_x, Ustar_y;

  vfield_to_vctrs(  vfield_list::Ustar , Ustar_x, Ustar_y );

  VectorXd U_x, U_y;

  FT ddt = dt;
//  if( dt < 1e-10 ) ddt = 1;  // for debugging, mainly
  
  U_x = Ustar_x.array() - ddt * grad_sx.array() / vol.array()  ;
  U_y = Ustar_y.array() - ddt * grad_sy.array() / vol.array() ;

  vctrs_to_vfield( U_x, U_y , vfield_list::U );

}

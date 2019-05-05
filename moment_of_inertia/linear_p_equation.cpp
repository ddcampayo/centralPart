//#include"pParticles.h"
#include"linear.h"
#include"fields_enum.h"
#include"simu.h"


// Solve for pressure
// The famous pressure Poisson equation

void linear::p_equation(const FT dt ) {

  cout << "Solving pressure equation " << endl;

  //  fill_Delta_DD(); // This may be important -- or not

  FT ddt = dt;
  if( dt < 1e-10 ) ddt = 1;  // for debugging, mainly

  //  volumes( T );

  VectorXd vol0  = field_to_vctr( sfield_list::vol0 ) ;
  VectorXd vol   = field_to_vctr( sfield_list::vol ) ;

  //FT target_vol_val =  simu.meanV() ;

  //  FT target_vol_val = vol.array().sum() / FT( vol.size() );
  //VectorXd Dvol = vol.array() - target_vol_val  ;


  VectorXd Dvol = vol - vol0  ;

  FT Dvol_sigma =  Dvol.array().square().sum() ; // / FT( vol.size() );
  FT Dvol_mean  =  vol.array().square().sum() ; // / FT( vol.size() );

  cout << "Pressure  "
       << " rel Dvol std dev: " << sqrt( Dvol_sigma / Dvol_mean )
       << endl;
  

  VectorXd Dp  =  DD_solver.solve( Dvol );

  VectorXd p0  = field_to_vctr( sfield_list::p ) ;

  vctr_to_field( p0 + Dp / ( ddt * ddt) , sfield_list::p  ) ;

  return;
}


void linear::u_add_grad_p( const FT dt ) {

  VectorXd grad_Px,grad_Py;

  DD_times_sfield( sfield_list::p  ,  grad_Px, grad_Py);

  VectorXd vol  = field_to_vctr( sfield_list::vol );
  // perhaps mean vol would be just fine
  
  VectorXd Ustar_x, Ustar_y;

  vfield_to_vctrs(  vfield_list::Ustar , Ustar_x, Ustar_y );

  VectorXd U_x, U_y;

  FT ddt = dt;
//  if( dt < 1e-10 ) ddt = 1;  // for debugging, mainly

  U_x = Ustar_x.array() - ddt * grad_Px.array() / vol.array()  ;
  U_y = Ustar_y.array() - ddt * grad_Py.array() / vol.array() ;

  vctrs_to_vfield( U_x, U_y , vfield_list::U );

}

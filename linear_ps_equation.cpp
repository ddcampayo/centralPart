//#include"pParticles.h"
#include"linear.h"
#include"fields_enum.h"
#include"simu.h"


// Solve for pressure
// The famous pressure Poisson equation

void linear::ps_equation(const FT dt ) {

  cout << "Solving pressure and s equations " << endl;

  //  fill_Delta_DD(); // This may be important -- or not

  FT ddt = dt;
  if( dt < 1e-10 ) ddt = 1;  // for debugging, mainly

  //  volumes( T );

  VectorXd vol  = field_to_vctr( sfield_list::vol ) ;

  //FT target_vol_val =  simu.meanV() ;

  //  FT target_vol_val = vol.array().sum() / FT( vol.size() );
  //VectorXd Dvol = vol.array() - target_vol_val  ;

  VectorXd vol0  = field_to_vctr( sfield_list::vol0 ) ;
  VectorXd Dvol = vol - vol0  ;

  FT Dvol_sigma =  Dvol.array().square().sum() ; // / FT( vol.size() );
  FT Dvol_mean  =  vol.array().square().sum() ; // / FT( vol.size() );

  cout << "Pressure  "
       << " rel Dvol std dev: " << sqrt( Dvol_sigma / Dvol_mean )
       << endl;

  VectorXd I0  = field_to_vctr( sfield_list::I0 ) ;
  VectorXd I  = field_to_vctr( sfield_list::I ) ;

  VectorXd DI = I - I0 ;

  FT DI_sigma =  DI.array().square().sum() ;
  FT I_mean  =  I.array().square().sum() ;

  cout << " s field  "
       << " rel DI std dev: " << sqrt( DI_sigma / I_mean )
       << endl;

  int N = vol0.size(); 
  
  VectorXd DvolDI( 2*N );
  DvolDI << Dvol , DI ;

  VectorXd DpDs  =  DDMM_solver.solve( DvolDI );

  VectorXd Dp = DpDs.head( N );
  VectorXd Ds = DpDs.tail( N );
  
  VectorXd p0  = field_to_vctr( sfield_list::p ) ;

  vctr_to_field( p0 + Dp / ( ddt * ddt) , sfield_list::p  ) ;

  VectorXd s0  = field_to_vctr( sfield_list::s ) ;

  vctr_to_field( s0 + Ds / ( ddt * ddt) , sfield_list::s  ) ;


  // no-no:  //  vctr_to_field( vol , sfield_list::vol0 );

  return;
}




void linear::u_add_grads( const FT dt ) {

  VectorXd gradPx,gradPy;

  DD_times_sfield( sfield_list::p  ,  gradPx, gradPy);

  VectorXd vol  = field_to_vctr( sfield_list::vol );
  // perhaps mean vol would be just fine

  VectorXd gradsx,gradsy;

  MM_times_sfield( sfield_list::s  ,  gradsx, gradsy);
  
  VectorXd Ustar_x, Ustar_y;

  vfield_to_vctrs(  vfield_list::Ustar , Ustar_x, Ustar_y );

  VectorXd U_x, U_y;

  FT ddt = dt;
//  if( dt < 1e-10 ) ddt = 1;  // for debugging, mainly

  VectorXd gradx = gradPx;// + gradsx ;
  VectorXd grady = gradPy;// + gradsy ;
  
  U_x = Ustar_x.array() - ddt * gradx.array() / vol.array()  ;
  U_y = Ustar_y.array() - ddt * grady.array() / vol.array() ;

  vctrs_to_vfield( U_x, U_y , vfield_list::U );

}

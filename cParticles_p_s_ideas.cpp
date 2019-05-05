// pParticles
// no weigths are used, basically
// Pressure is used in order to enforce incompressibility
// An additional field, s, is used to enforce const moments of inertia

#include"cParticles.h"
#include"linear.h"
#include"simu.h"

sim_data simu;

int main() {

  // TODO: read parameter file
  // TODO: read better from parameter file
  
  int init_max_iters; cin >> init_max_iters; // = 40;
  FT  init_tol2 ; cin >> init_tol2 ; // = 1e-3; 

  int inner_max_iters; cin >> inner_max_iters; // = 10; 
  FT  disp_tol; cin >> disp_tol; //  = 1e-6;

  //  int s_iters; cin >> s_iters; //= 10;

  FT total_time =  1/( 2 * 3.14 * 0.2) ;

  const std::string particle_file("particles.dat");
  const std::string diagram_file("diagram.dat");

  Triangulation T;

  cout << "Creating point cloud" << endl;

  //simu.do_perturb(0.01);
  create( T , 1.0 );
  number( T );

  //  set_vels_rotating( T );
  //  set_vels_Lamb_Oseen( T );

  linear algebra( T );

  // Init loop!

  int init_iter=1;
  
  for( ; init_iter <= init_max_iters ; ++init_iter) {
  
    volumes( T ); 

    FT dd = lloyds( T ) ;

    cout << " init loop , iter " << init_iter << " dd = " << dd << endl;
    if( dd < init_tol2) break;

  }

  cout << "Init loop converged in " << init_iter << " steps " << endl;
  
  set_vels_Gresho( T );

  volumes( T ); 
  
  FT d0;
  FT dt=0.001;

  cin >> dt ;

  simu.set_dt( dt );

  // half-step leapfrog
  FT dt2 = dt / 2.0 ;

  // whole step
  //FT dt2 = dt  ;

  draw( T , particle_file     );

  draw_diagram( T , diagram_file );
  
  std::ofstream log_file;
  log_file.open("main.log");

  do {
    simu.next_step();
    simu.advance_time( );

    FT displ;

    int in_iter = 1;// , s_it = 1;

    volumes( T ); 

    backup( T );

    algebra.u_star( );
    
    //    s_it = 1;

    algebra.reset_s();
    algebra.reset_p();

    for ( ; in_iter <= inner_max_iters ; in_iter++) {

      displ = move( T , dt , d0 );

      // whole step, special 1st time
      //      FT ddt = dt2 ;

      cout
	<< "********" << endl
	<< " p and s Iter  " << in_iter
	<< " . Moved from previous (rel.): " << displ <<
	" ; from original (rel.): " << d0
	<< endl ;
  
      volumes( T ); 

      algebra.fill_matrices();

      algebra.ps_equation( dt );

      algebra.u_add_grad_ps( dt2 );
      
      if( displ < disp_tol ) break;


      //      if( simu.current_step() == 1 ) ddt = dt / 2 ;
      
      	// optional .- inner ps, p, and s loops

      //           algebra.u_star( );

      // for ( int p_iter=1 ; p_iter <= inner_max_iters ; p_iter++) {

      //   volumes( T ); 

      // 	algebra.fill_matrices();

      // 	algebra.ps_equation( dt2 );
      // 	algebra.u_add_grad_ps( ddt );
      // 	//	algebra.u_add_grad_ps( ddt );

      // 	displ = move( T , dt2 , d0 );

      // 	cout
      // 	  << "********" << endl
      // 	  << "p - s Iter  " << p_iter
      // 	  << " . Moved from previous (rel.): " << displ <<
      // 	  " ; from original (rel.): " << d0
      // 	  << endl ;
      
      // 	if( displ < disp_tol ) break;

      // }

      
      // algebra.u_star( );
      
      // for ( int p_iter=1 ; p_iter <= inner_max_iters ; p_iter++) {
      // 	volumes( T ); 

      // 	algebra.fill_matrices();

      // 	algebra.p_equation( dt2 );
      // 	algebra.u_add_grad_p( ddt );
      // 	//	algebra.u_add_grad_ps( ddt );

      // 	displ = move( T , dt2 , d0 );

      // 	cout
      // 	  << "********" << endl
      // 	  << "p Iter  " << p_iter
      // 	  << " . Moved from previous (rel.): " << displ <<
      // 	  " ; from original (rel.): " << d0
      // 	  << endl ;
      
      // 	if( displ < disp_tol ) break;

      // }

      // algebra.u_star( );
      
      // for ( int p_iter=1 ; p_iter <= inner_max_iters ; p_iter++) {
      // 	volumes( T ); 

      // 	algebra.fill_matrices();

      // 	algebra.s_equation( dt2 );
      // 	algebra.u_add_grad_s( ddt );
      // 	//	algebra.u_add_grad_ps( ddt );

      // 	displ = move( T , dt2 , d0 );

      // 	cout
      // 	  << "********" << endl
      // 	  << "s Iter  " << p_iter
      // 	  << " . Moved from previous (rel.): " << displ <<
      // 	  " ; from original (rel.): " << d0
      // 	  << endl ;
      
      // 	if( displ < disp_tol ) break;

      // }

      // displ = move( T , dt2 , d0 );

      // volumes( T ); 

      // algebra.fill_matrices();

      // algebra.ps_equation( dt2 );
      // algebra.u_add_grad_ps( ddt );

      // displ = move( T , dt2 , d0 );

      // cout
      //  	<< "********" << endl
      //  	<< "p-s Iter  " << in_iter
      //  	<< " . Moved from previous (rel.): " << displ <<
      //  	" ; from original (rel.): " << d0
      //  	<< endl ;
      
      //      if( displ < disp_tol ) break;

    }

    //    displ = move( T , dt , d0 );
    //    volumes( T ); 

    cout
      << "Whole step  "
      << " : disp " << displ << endl ;

    // half-step:
    //    update_full_vel( T );
    algebra.u_add_grad_ps( dt );

    draw( T , particle_file     );
    draw_diagram( T , diagram_file );

    log_file
      << simu.current_step() << "  "
      << simu.time() << "  "
      << " iters = " << in_iter
      << " L2_vel =  " << L2_vel_Gresho(T)
      << endl ;

  } while ( simu.time() < total_time );

  log_file.close();

  return 0;

}

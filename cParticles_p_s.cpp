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

  int s_iters; cin >> s_iters; //= 10;

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
  //  FT dt2 = dt / 2.0 ;

  // whole step
  FT dt2 = dt  ;

  draw( T , particle_file     );

  draw_diagram( T , diagram_file );
  
  std::ofstream log_file;
  log_file.open("main.log");

  do {
    simu.next_step();
    simu.advance_time( );

    FT displ;

    int in_iter = 1 , s_it = 1;
    
    for ( ; in_iter <= inner_max_iters ; in_iter++) {

      //   moment of inertia (s)   iteration

      volumes( T ); 

      backup( T );

      algebra.u_star( );
    
      s_it = 1;
     
      displ = move( T , dt2 , d0 );

      algebra.reset_s();
      algebra.reset_p();

      for ( ; s_it <= s_iters ; s_it++) {
	volumes( T ); 

	algebra.fill_matrices();

	algebra.ps_equation( dt2 );

	// whole step, special 1st time
	FT ddt = dt2;
	if( simu.current_step() == 1 ) ddt = dt / 2 ;

	algebra.u_add_grads( ddt );

	displ = move( T , dt2 , d0 );

	cout
	  << "********" << endl
	  << "p-s Iter  " << s_it
	  << " . Moved from previous (rel.): " << displ <<
	  " ; from original (rel.): " << d0
	  << endl ;

	//      if( displ < disp_tol ) break;

      }

      volumes( T ); 

 
      if( displ < disp_tol ) break;


    }
    ///////////// end   p iter

    
    cout
      << "Whole step  "
      << " : disp " << displ << endl ;

    volumes( T ); 

    // half-step:
//    update_full_vel( T );

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

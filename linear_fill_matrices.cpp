#include"linear.h"

void linear::fill_matrices(void){ 

  std::cout << " Filling  matrices" << std::endl;

  //  int n=simu.no_of_points();

  std::vector<triplet> dx, dy, mx, my ; // , bb ;            // list of non-zeros coefficients

  // diagonal terms

  typedef std::map<int,FT> diag_map;

  diag_map   diag_dx, diag_dy;
  diag_map   diag_mx, diag_my;

  int N=1;

  F_e_it eit = T.finite_edges_begin();

  for ( ; eit !=T.finite_edges_end(); ++eit) {

    Face_handle f =  eit -> first ;
    int i0 = eit -> second;

    Vertex_handle vi = f->vertex( (i0+1) % 3);
    Vertex_handle vj = f->vertex( (i0+2) % 3);

    int i = vi->idx();
    int j = vj->idx();

    //    if( (i < 0 ) || ( j < 0) ) continue;
    
    CGAL::Object o = T.dual(eit);

    const Segment * Vor_segment = CGAL::object_cast<Segment>( &o );

    if (! Vor_segment ) continue;

    FT Aij = std::sqrt( Vor_segment->squared_length() );

    //    cout << " A0 = " << Aij << endl;
    
    Point pi = vi->point().point();
    Point pj = vj->point().point();

    Vector_2 eij = pj - pi;

    FT l2ij =  eij.squared_length() ;
    FT lij = std::sqrt( l2ij );
    
    Point bij = CGAL::midpoint( Vor_segment->source() , Vor_segment->target() );

    FT voli = vi->vol();
    FT volj = vj->vol();

    Vector_2 r_ij_j = pj - bij;
    Vector_2 r_ij_i = pi - bij;

    Vector_2 DDij = Aij / lij * r_ij_j; // ( pj - bij);
    Vector_2 DDji = Aij / lij * r_ij_i; // ( pi - bij);

    Point bi = vi->centroid.val();
    Point bj = vj->centroid.val();

    Vector_2 di = ( bi - pi ) * voli ;
    Vector_2 dj = ( bj - pj ) * volj ;

    Vector_2 di_perp = ( (di * eij) / l2ij ) * eij ;
    Vector_2 di_para = di - di_perp;

    Vector_2 dj_perp = ( (dj * eij) / l2ij ) * eij ;
    Vector_2 dj_para = dj - dj_perp;
 
    // FT r2_ij_j = rr_ij_j.squared_length();  // (these two are the same on Voronoi)
    // FT r2_ij_i = rr_ij_i.squared_length();
    
    Vector_2 MMij =  -2 * Aij / lij * (
				       (di * r_ij_i ) * r_ij_j
				       + ( Aij*Aij / 12 ) * di_para
				       );

    Vector_2 MMji =  -2 * Aij / lij * (
				       (dj * r_ij_j ) * r_ij_i
				       + ( Aij*Aij / 12 ) * dj_para
				 );

    if( (i >= 0 ) && ( j >= 0) ) {

      dx.push_back( triplet( i, j,  DDij.x() ));
      dy.push_back( triplet( i, j,  DDij.y() ));

      dx.push_back( triplet( j, i,  DDji.x() ));
      dy.push_back( triplet( j, i,  DDji.y() ));

      mx.push_back( triplet( i, j,  MMij.x() ));
      my.push_back( triplet( i, j,  MMij.y() ));

      mx.push_back( triplet( j, i,  MMji.x() ));
      my.push_back( triplet( j, i,  MMji.y() ));

    }
    
    if (i >= 0 ) {
      diag_dx[ i ] -= DDji.x();
      diag_dy[ i ] -= DDji.y();

      //      Vector_2 MMii = - MMij;
      Vector_2 MMii = 2 * Aij / lij * (
				       (di * r_ij_i ) * r_ij_i
				       + ( Aij*Aij / 12 ) * di_para
				       );

      diag_mx[ i ] += MMii.x();
      diag_my[ i ] += MMii.y();

//      di_m_x[ i ] -= MMij.x();
//      di_m_y[ i ] -= MMij.y();
 
    }
    if (j >= 0 ) {

      diag_dx[ j ] -= DDij.x();
      diag_dy[ j ] -= DDij.y();

      //      Vector_2 MMjj = - MMji;

      Vector_2 MMjj =  2 * Aij / lij * (
      					(dj * r_ij_j ) * r_ij_j
      					+ ( Aij*Aij / 12 ) * dj_para
      					);

      diag_mx[ j ] += MMjj.x();
      diag_my[ j ] += MMjj.y();

      //      di_d_x[ j ] -= DDji.x();
      //      di_d_y[ j ] -= DDji.y();

      //      di_m_x[ j ] -= MMji.x();
      //      di_m_y[ j ] -= MMji.y();


    }

    // if( (i!=0) && (j!=0) ) {
    //   aa.push_back( triplet(i - 1 , j -1 ,  ddelta ));
    //   aa.push_back( triplet(j - 1 , i -1 ,  ddelta ));

    //   ++N;
    // }

    if( i+1 > N ) { N = i+1 ; } // keep maximum

    //    cout << i << "  " << j << "  " << ddelta << endl;
  }



  // include "spring" term in M
  for(F_v_it fv=T.finite_vertices_begin();
      fv!=T.finite_vertices_end();
      fv++)  {

    int idx = fv->idx.val();

    if( idx < 0 ) continue;

    FT voli = fv->vol();
    Point ri = fv->point().point();
    Point bi = fv->centroid.val();

    Vector_2 dd = 2 * voli * voli * ( bi - ri ) ;
    diag_mx[ idx ] -= dd.x();
    diag_my[ idx ] -= dd.y();

  }
  
  // Add diagonal terms .-
  
  for( diag_map::const_iterator it = diag_dx.begin(); it != diag_dx.end(); ++it ) {
    int i    = it->first;
    FT diagx = it->second;
    dx.push_back( triplet( i, i,  diagx ));
  }

  
  for( diag_map::const_iterator it = diag_dy.begin(); it != diag_dy.end(); ++it ) {
    int i    = it->first;
    FT diagy = it->second;
    dy.push_back( triplet( i, i,  diagy ));
  }


  for( diag_map::const_iterator it = diag_mx.begin(); it != diag_mx.end(); ++it ) {
    int i    = it->first;
    FT diagx = it->second;
    mx.push_back( triplet( i, i,  diagx ));
  }

  
  for( diag_map::const_iterator it = diag_my.begin(); it != diag_my.end(); ++it ) {
    int i    = it->first;
    FT diagy = it->second;
    my.push_back( triplet( i, i,  diagy ));
  }
    
  Dx.resize( N , N );
  Dx.setFromTriplets(dx.begin(), dx.end());
  // std::cout << " Filled DDx  matrix" << std::endl;
  // cout << "matrix size " << DDx.rows() << " times " << DDx.cols() << endl;

  Dy.resize( N , N );
  Dy.setFromTriplets(dy.begin(), dy.end());
  // std::cout << " Filled DDy  matrix" << std::endl;
  // cout << "matrix size " << DDy.rows() << " times " << DDy.cols() << endl;

  Mx.resize( N , N );
  Mx.setFromTriplets(mx.begin(), mx.end());

  My.resize( N , N );
  My.setFromTriplets(my.begin(), my.end());
  
  
  // set up solvers .-
   
  VectorXd vol  = field_to_vctr( sfield_list::vol ) ;
  VectorXd inv_vol  = 1.0 / vol.array() ;

  //  SpMat
  DD =
    - Dx * inv_vol.asDiagonal() * Dx.transpose()
    - Dy * inv_vol.asDiagonal() * Dy.transpose();

  //SpMat
  MM =
    - Mx * inv_vol.asDiagonal() * Mx.transpose()
    - My * inv_vol.asDiagonal() * My.transpose();

  //SpMat
  DM =
    - Dx *( inv_vol.asDiagonal() * Mx.transpose() )
    - Dy *( inv_vol.asDiagonal() * My.transpose() );

  //SpMat
  MD =
    - Mx *( inv_vol.asDiagonal() * Dx.transpose() )
    - My *( inv_vol.asDiagonal() * Dy.transpose() );


  DD_solver.compute( DD );

  if(  DD_solver.info() != Eigen::Success ) {
    std::cout << "Failure decomposing DD matrix " << endl;
  }

  MM_solver.compute( MM );

  if(  MM_solver.info() != Eigen::Success ) {
    std::cout << "Failure decomposing MM matrix " << endl;
  }
  
  // https://stackoverflow.com/questions/28685877/convert-an-eigen-matrix-to-triplet-form-c
  std::vector<triplet> ddmm ;
  for (int k=0; k < DD.outerSize(); ++k)
    for (SpMat::InnerIterator it(DD,k); it; ++it) {
	  int i = it.row();
	  int j = it.col();
	  
	  FT val= it.value();
 
	  ddmm.push_back( triplet( i, j,  val ));
    }
  
 // for (int k=0; k < DM.outerSize(); ++k)
 //   for (SpMat::InnerIterator it(DM,k); it; ++it) {
 //     int i = it.row();
 //     int j = it.col();
	  
 //     FT val= it.value();
 
 //     ddmm.push_back( triplet( i, j + N ,  val ));
 //   }

 // for (int k=0; k < MD.outerSize(); ++k)
 //   for (SpMat::InnerIterator it(MD,k); it; ++it) {
 //     int i = it.row();
 //     int j = it.col();
	  
 //     FT val= it.value();
 
 //     ddmm.push_back( triplet( i + N , j,  val ));
 //   }

  for (int k=0; k < MM.outerSize(); ++k)
    for (SpMat::InnerIterator it(MM,k); it; ++it) {
      int i = it.row();
      int j = it.col();
	 
      FT val= it.value();
 
      ddmm.push_back( triplet( i + N , j + N ,  val ));
    }

  DDMM.resize( 2*N , 2*N );
  DDMM.setFromTriplets( ddmm.begin() , ddmm.end() );

  DDMM_solver.compute( DDMM );

  if(  DDMM_solver.info() != Eigen::Success ) {
    std::cout << "Failure decomposing DDMM matrix " << endl;
  }

  
  return;

}


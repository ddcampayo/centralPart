#include"linear.h"
#include"fields_enum.h"

// NOTE: these two are the "divergence" (but for a 1/V factor)
void linear::DD_scalar_vfield(const vfield_list::take from , const sfield_list::take to )
{

  VectorXd vx, vy;
  vfield_to_vctrs( from , vx, vy );

  VectorXd div = Dx * vx + Dy*vy;

  vctr_to_field( div , to );

  return;
}

VectorXd linear::DD_scalar_vfield(const vfield_list::take from )
{

  VectorXd vx, vy;

  vfield_to_vctrs( from , vx, vy );

  // cout << "vx cols " << vx.cols() << endl;
  // cout << "vx rows " << vx.rows() << endl;

  // cout << Dx << endl;
  // cout << "vx " << endl;
  // cout << vx << endl;


  return Dx * vx  + Dy * vy ;
}


// NOTE: this is the "gradient" (but for a 1/V factor)
// it features a minus sign, and transposition !!
void linear::DD_times_sfield(const sfield_list::take from ,
			     VectorXd& Dx,VectorXd& Dy)
{

  VectorXd p = field_to_vctr( from );

  Dx = -Dx.transpose() * p;
  Dy = -Dy.transpose() * p;

  return;
}

void linear::MM_times_sfield(const sfield_list::take from ,
			     VectorXd& Dx,VectorXd& Dy)
{

  VectorXd s = field_to_vctr( from );

  Dx = -Mx.transpose() * s;
  Dy = -My.transpose() * s;

  return;
}


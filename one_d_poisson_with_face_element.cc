//kruemelmonster
//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC//    Version 1.0; svn revision $LastChangedRevision: 1307 $
//LIC//
//LIC// $LastChangedDate: 2018-01-18 11:30:14 +0000 (Thu, 18 Jan 2018) $
//LIC// 
//LIC// Copyright (C) 2006-2016 Matthias Heil and Andrew Hazel
//LIC// 
//LIC// This library is free software; you can redistribute it and/or
//LIC// modify it under the terms of the GNU Lesser General Public
//LIC// License as published by the Free Software Foundation; either
//LIC// version 2.1 of the License, or (at your option) any later version.
//LIC// 
//LIC// This library is distributed in the hope that it will be useful,
//LIC// but WITHOUT ANY WARRANTY; without even the implied warranty of
//LIC// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//LIC// Lesser General Public License for more details.
//LIC// 
//LIC// You should have received a copy of the GNU Lesser General Public
//LIC// License along with this library; if not, write to the Free Software
//LIC// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
//LIC// 02110-1301  USA.
//LIC// 
//LIC// The authors may be contacted at oomph-lib@maths.man.ac.uk.
//LIC// 
//LIC//====================================================================
//Driver for a simple 1D poisson problem

// Generic oomph-lib routines
#include "generic.h"

// Include Poisson elements/equations
#include "poisson.h"

// Include the mesh
#include "meshes/one_d_mesh.h"

#include "poisson_sing_face_element.h"

using namespace std;

using namespace oomph;

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

namespace oomph
{

// //====================================================================
// /// New class. Mainly overloads output-related functions to add
// /// "singular function" (which is assumed to satisfy the Laplace
// /// equation; therefore no change to the governing (bulk) equations) 
// /// to the FE solution. 
// //====================================================================
// template<unsigned DIM, unsigned NNODE_1D>
// class MyQPoissonElement : public virtual QPoissonElement<DIM,NNODE_1D>
// {

// public:

//  /// Constructor
//  MyQPoissonElement() : Poisson_sing_el_pt(0)
//   {}

//  /// \short Return FE representation of function value u_poisson(s) 
//  /// plus scaled singular fct (if provided) at local coordinate s
//  inline double interpolated_u_poisson(const Vector<double> &s) const
//   {
//    double u_fe=QPoissonElement<DIM,NNODE_1D>::interpolated_u_poisson(s);
//    if (Poisson_sing_el_pt!=0)
//     {     
//      Vector<double> x(DIM);
//      for(unsigned i=0;i<DIM;i++) 
//       {
//        x[i]=this->interpolated_x(s,i);
//       }
//      u_fe+=Poisson_sing_el_pt->singular_fct(x);
//     }
//    return u_fe;
//   } 


//  /// Output with various contributions
//  void  output_with_various_contributions(std::ostream &outfile, 
//                                          const unsigned &nplot)
//   {
//    //Vector of local coordinates
//    Vector<double> s(DIM);
   
//    // Tecplot header info
//    outfile << this->tecplot_zone_string(nplot);
   
//    // Loop over plot points
//    unsigned num_plot_points=this->nplot_points(nplot);
//    for (unsigned iplot=0;iplot<num_plot_points;iplot++)
//     {
//      // Get local coordinates of plot point
//      this->get_s_plot(iplot,nplot,s);
     
//      Vector<double> x(DIM);
//      for(unsigned i=0;i<DIM;i++) 
//       {
//        x[i]=this->interpolated_x(s,i);
//        outfile << x[i] << " ";
//       }
//      double u_sing=0.0;
//      if (Poisson_sing_el_pt!=0)
//       {
//        u_sing=Poisson_sing_el_pt->singular_fct(x);
//       }
//      outfile << this->interpolated_u_poisson(s) << " "
//              << QPoissonElement<DIM,NNODE_1D>::interpolated_u_poisson(s) << " "
//              << u_sing << " "
//              << std::endl;   
//     }
   
//    // Write tecplot footer (e.g. FE connectivity lists)
//    this->write_tecplot_zone_footer(outfile,nplot);
   
//   }

//  /// Pointer to element that stores singular fct
//  TemplateFreePoissonWithSingularityBCFaceElementBase*& poisson_sing_el_pt()
//   {
//    return Poisson_sing_el_pt;
//   }
 
// private:

//  /// hierher Pointer to element that stores singular fct
//  TemplateFreePoissonWithSingularityBCFaceElementBase* Poisson_sing_el_pt;
 
// };


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////



// //=======================================================================
// /// Face geometry for the MyQPoissonElement elements: The spatial 
// /// dimension of the face elements is one lower than that of the
// /// bulk element but they have the same number of points
// /// along their 1D edges.
// //=======================================================================
// template<unsigned DIM, unsigned NNODE_1D>
// class FaceGeometry<MyQPoissonElement<DIM,NNODE_1D> >: 
//  public virtual QElement<DIM-1,NNODE_1D>
// {

//   public:
 
//  /// \short Constructor: Call the constructor for the
//  /// appropriate lower-dimensional QElement
//  FaceGeometry() : QElement<DIM-1,NNODE_1D>() {}

// };

// ////////////////////////////////////////////////////////////////////////
// ////////////////////////////////////////////////////////////////////////
// ////////////////////////////////////////////////////////////////////////


// //=======================================================================
// /// Face geometry for the 1D MyQPoissonElement elements: Point elements
// //=======================================================================
// template<unsigned NNODE_1D>
// class FaceGeometry<MyQPoissonElement<1,NNODE_1D> >: 
//  public virtual PointElement
// {

//   public:
 
//  /// \short Constructor: Call the constructor for the
//  /// appropriate lower-dimensional QElement
//  FaceGeometry() : PointElement() {}

// };

}

//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////



//==start_of_namespace================================================
/// Namespace for fish-shaped solution of 1D Poisson equation
//====================================================================
namespace FishSolnOneDPoisson
{

 /// \short Sign of the source function 
 /// (- gives the upper half of the fish, + the lower half)
 int Sign=-1;


 /// Exact, fish-shaped solution as a 1D vector
 void get_exact_u(const Vector<double>& x, Vector<double>& u)
 {
  u[0] = double(Sign)*((sin(sqrt(30.0))-1.0)*x[0]-sin(sqrt(30.0)*x[0]));
 }


 /// Source function required to make the fish shape an exact solution 
 void source_function(const Vector<double>& x, double& source)
 {
  source = double(Sign)*30.0*sin(sqrt(30.0)*x[0]);
 }



 /// \short "Singular" function and gradient
 void singular_fct_and_gradient(const Vector<double>& x,
                                double& u, Vector<double>& du_dx)
 {
  // a linear fct which solves the Laplace eqn
  u=-3.5+15.3*x[0];
  du_dx[0]=15.3;
 }


 /// \short "Singular" function
 double singular_fct(const Vector<double>& x)
 {
  double u=0.0;
  Vector<double> du_dx(1);
  singular_fct_and_gradient(x,u,du_dx);
  return u;
 }

 /// \short Gradient of "Singular" function
 Vector<double> gradient_of_singular_fct(const Vector<double>& x)
 {
  double u=0.0;
  Vector<double> du_dx(1);
  singular_fct_and_gradient(x,u,du_dx);
  return du_dx;
 }



} // end of namespace







//==start_of_problem_class============================================
/// 1D Poisson problem in unit interval.
//====================================================================
template<class ELEMENT> 
class OneDPoissonProblem : public Problem
{

public:

 /// Constructor: Pass number of elements and pointer to source function
 OneDPoissonProblem(const unsigned& n_element, 
                    PoissonEquations<1>::PoissonSourceFctPt source_fct_pt);

 /// Destructor (empty)
 ~OneDPoissonProblem(){}

 /// Update the problem specs before solve: (Re)set boundary conditions
 void actions_before_newton_solve();

 /// Update the problem specs after solve (empty)
 void actions_after_newton_solve(){}

 /// \short Doc the solution, pass the number of the case considered,
 /// so that output files can be distinguished.
 void doc_solution(const unsigned& label);

private:

 /// Pointer to source function
 PoissonEquations<1>::PoissonSourceFctPt Source_fct_pt;

 /// Bulk mesh
 Mesh* Bulk_mesh_pt;

 /// Face element mesh
 Mesh* Face_mesh_for_flux_jump_pt;

 /// Face element mesh for BC
 Mesh* Face_mesh_for_bc_pt;

}; // end of problem class





//=====start_of_constructor===============================================
/// \short Constructor for 1D Poisson problem in unit interval.
/// Discretise the 1D domain with n_element elements of type ELEMENT.
/// Specify function pointer to source function. 
//========================================================================
template<class ELEMENT>
OneDPoissonProblem<ELEMENT>::OneDPoissonProblem(const unsigned& n_element,
 PoissonEquations<1>::PoissonSourceFctPt source_fct_pt) : 
 Source_fct_pt(source_fct_pt)
{ 
 Problem::Sparse_assembly_method = Perform_assembly_using_two_arrays;

// Problem::Problem_is_nonlinear = false;
 // Set domain length 
 double L=1.0;

 // Build mesh and store pointer in Problem
 Bulk_mesh_pt = new OneDMesh<ELEMENT>(n_element,L);


 // Left element where we have jump between augmented and non-augmented
 //--------------------------------------------------------------------
 // elements
 //---------

 // Map to keep track of mapping between old and duplicated nodes
 // (not a big deal here, but needed in higher dimensions)
 std::map<Node*,Node*> existing_duplicate_node_pt;

 unsigned left_element_number=unsigned(double(n_element)*0.5);
 oomph_info << "left element number: " << left_element_number << std::endl;

 ELEMENT* left_el_pt=dynamic_cast<ELEMENT*>(Bulk_mesh_pt->
                                            element_pt(left_element_number));

 // Build the corresponding face element
 int face_index=1;
 PoissonWithSingularityFluxJumpFaceElement<ELEMENT>* face_element_pt = new 
  PoissonWithSingularityFluxJumpFaceElement<ELEMENT>
  (left_el_pt,face_index,existing_duplicate_node_pt);
 
 //Add the prescribed-flux element to the mesh
 Face_mesh_for_flux_jump_pt=new Mesh;
 Face_mesh_for_flux_jump_pt->add_element_pt(face_element_pt);

 // Only one node to be added here
 Face_mesh_for_flux_jump_pt->add_node_pt(face_element_pt->node_pt(0));


 
 // Leftmost element where Dirichlet condition is to be applied
 //------------------------------------------------------------
 unsigned left_most_element_number=0;
 ELEMENT* left_most_el_pt=dynamic_cast<ELEMENT*>(
  Bulk_mesh_pt->element_pt(left_most_element_number));
 
 // Build the corresponding face element
 face_index=-1;
 PoissonWithSingularityBCFaceElement<ELEMENT>* bc_face_element_pt = new 
  PoissonWithSingularityBCFaceElement<ELEMENT>(left_most_el_pt,face_index);
 
 //Add the face element to the mesh
 Face_mesh_for_bc_pt=new Mesh;
 Face_mesh_for_bc_pt->add_element_pt(bc_face_element_pt);
 

 // Build global mesh
 add_sub_mesh(Bulk_mesh_pt);
 add_sub_mesh(Face_mesh_for_flux_jump_pt);
 add_sub_mesh(Face_mesh_for_bc_pt);
 build_global_mesh();

 // Set the boundary conditions for this problem: By default, all nodal
 // values are free -- we only need to pin the ones that have 
 // Dirichlet conditions. 

 // Pin the single nodal value at the single node on mesh 
 // boundary 0 (= the left domain boundary at x=0)
 // hierher if this... Bulk_mesh_pt->boundary_node_pt(0,0)->pin(0);
 
 // Pin the single nodal value at the single node on mesh 
 // boundary 1 (= the right domain boundary at x=1)
 Bulk_mesh_pt->boundary_node_pt(1,0)->pin(0);

 // Complete the setup of the 1D Poisson problem:

 // Loop over elements and set pointers to source function
 for(unsigned i=0;i<n_element;i++)
  {
   // Upcast from GeneralisedElement to the present element
   ELEMENT *elem_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(i));
   
   //Set the source function pointer
   elem_pt->source_fct_pt() = Source_fct_pt;
   
   // Pass "singular" fct via pointer to element
   if (i<=left_element_number)
    {
     elem_pt->poisson_sing_el_pt()=bc_face_element_pt;
    }

  }


 // Specify singular fct and its gradient for all elements
 // that need it
 bc_face_element_pt->unscaled_singular_fct_pt()=
  &FishSolnOneDPoisson::singular_fct;
 bc_face_element_pt->gradient_of_unscaled_singular_fct_pt()=
  &FishSolnOneDPoisson::gradient_of_singular_fct;

 face_element_pt->set_poisson_sing_el_pt(bc_face_element_pt);

 // Setup equation numbering scheme
 assign_eqn_numbers();

} // end of constructor




//===start_of_actions_before_newton_solve========================================
/// \short Update the problem specs before solve: (Re)set boundary values
/// from the exact solution. 
//========================================================================
template<class ELEMENT>
void OneDPoissonProblem<ELEMENT>::actions_before_newton_solve()
{
 
 // Assign boundary values for this problem by reading them out
 // from the exact solution.

 // Left boundary is node 0 in the mesh:
 Node* left_node_pt=Bulk_mesh_pt->boundary_node_pt(0,0);

 // Determine the position of the boundary node (the exact solution
 // requires the coordinate in a 1D vector!)
 Vector<double> x(1);
 x[0]=left_node_pt->x(0);
 
 // Boundary value (read in from exact solution which returns
 // the solution in a 1D vector)
 Vector<double> u(1);
 FishSolnOneDPoisson::get_exact_u(x,u);
 
 Vector<double> nodal_boundary_value(1);
 nodal_boundary_value[0]=u[0];
 dynamic_cast<PoissonWithSingularityBCFaceElement<ELEMENT>*>(
  Face_mesh_for_bc_pt->element_pt(0))->set_nodal_boundary_values
  (nodal_boundary_value);

 // Right boundary is last node in the mesh:
 Node* right_node_pt=Bulk_mesh_pt->boundary_node_pt(1,0);

 // Determine the position of the boundary node
 x[0]=right_node_pt->x(0);
 
 // Boundary value (read in from exact solution which returns
 // the solution in a 1D vector)
 FishSolnOneDPoisson::get_exact_u(x,u);
 
 // Assign the boundary condition to one (and only) nodal value
 right_node_pt->set_value(0,u[0]);

 
} // end of actions before solve



//===start_of_doc=========================================================
/// Doc the solution in tecplot format. Label files with label.
//========================================================================
template<class ELEMENT>
void OneDPoissonProblem<ELEMENT>::doc_solution(const unsigned& label)
{ 
 using namespace StringConversion;

 // Number of plot points
 unsigned npts;
 npts=5; 


 oomph_info << "Value of amplitude: "
            << dynamic_cast<PoissonWithSingularityBCFaceElement<ELEMENT>*>(
             Face_mesh_for_bc_pt->element_pt(0))->amplitude_of_singular_fct()
            << std::endl;

 ofstream extended_solution_file(("extended_soln" 
                                  + to_string(label) + ".dat").c_str());
 unsigned nel=Bulk_mesh_pt->nelement();
 for (unsigned e=0;e<nel;e++)
  {
   dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e))->
                          output_with_various_contributions(
                           extended_solution_file,npts);
  }
 extended_solution_file.close();



 // Output solution with specified number of plot points per element
 ofstream solution_file(("soln" + to_string(label) + ".dat").c_str());
 Bulk_mesh_pt->output(solution_file,npts);
 solution_file.close();

 // Output exact solution at much higher resolution (so we can
 // see how well the solutions agree between nodal points)
 ofstream exact_file(("exact_soln" + to_string(label) + ".dat").c_str());
 Bulk_mesh_pt->output_fct(exact_file,20*npts,FishSolnOneDPoisson::get_exact_u); 
 exact_file.close();

 // Doc pointwise error and compute norm of error and of the solution
 double error,norm;
 ofstream error_file(("error" + to_string(label) + ".dat").c_str());
 Bulk_mesh_pt->compute_error(error_file,FishSolnOneDPoisson::get_exact_u,
                          error,norm); 
 error_file.close();

 // Doc error norm:
 cout << "\nNorm of error    : " << sqrt(error) << std::endl; 
 cout << "Norm of solution : " << sqrt(norm) << std::endl << std::endl;
 cout << std::endl;

} // end of doc

 

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


//======start_of_main==================================================
/// Driver for 1D Poisson problem
//=====================================================================
int main()
{

 // Set up the problem: 
 // Solve a 1D Poisson problem using a source function that generates
 // a fish shaped exact solution
 unsigned n_element=40; //Number of elements
 OneDPoissonProblem<MyQPoissonElement<1,4> > 
  problem(n_element,FishSolnOneDPoisson::source_function);

 // Set the sign of the source function:
 cout << "\n\n\nSolving with negative sign:\n" << std::endl;
 FishSolnOneDPoisson::Sign=-1;

 // Solve the problem with this Sign
 problem.newton_solve();

 //Output solution for this case (label output files with "0")
 problem.doc_solution(0);


 // Change the sign of the source function:
 cout << "\n\n\nSolving with positive sign:\n" << std::endl;
 FishSolnOneDPoisson::Sign=1;

 // Re-solve the problem with this Sign (boundary conditions get
 // updated automatically when Problem::actions_before_newton_solve() is
 // called.
 problem.newton_solve();

 //Output solution for this case (label output files with "1")
 problem.doc_solution(1);


} // end of main










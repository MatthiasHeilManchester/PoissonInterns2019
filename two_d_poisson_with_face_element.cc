//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented,
//LIC// multi-physics finite-element library, available
//LIC// at http://www.oomph-lib.org.
//LIC//
//LIC//    Version 1.0; svn revision $LastChangedRevision: 1097 $
//LIC//
//LIC// $LastChangedDate: 2015-12-17 11:53:17 +0000 (Thu, 17 Dec 2015) $
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
//Driver for Poisson in backward step domain -- meshed with triangle

//Generic includes
#include "generic.h"
#include "poisson.h"

// The mesh
#include "meshes/triangle_mesh.h"


#include "poisson_sing_face_element.h"

using namespace std;

using namespace oomph;

//==start_of_namespace==============================
/// Namespace for physical parameters
//==================================================
namespace Global_Physical_Variables
{

 // Dimensionless domain values
 // ---------------------------

 /// Dimless width of inflow channel
 double H_up = 1.0; 

 /// Dimless length of channel upstream of step
 double L_up = 2.0; 

 /// \short Dimless length of channel downstream of step
 double L_down = 2.0; 

 /// Dimless width of outflow channel
 double H_down = 2.0; 

 /// Radius of internal boundary (surrounding the singularity)
 double Radius_of_internal_boundary=0.5;

// Bit of a hack but it facilitates reuse...
#include "unstructured_backward_step_mesh.h"

 /// \short "Singular" function and gradient
 void singular_fct_and_gradient(const Vector<double>& x,
                                double& u, Vector<double>& du_dx)
 {
  // Radius & polar angle
  double r=sqrt((x[0]-L_up)*(x[0]-L_up)+x[1]*x[1]);

  // Little hack to make sure atan2 doesn't overshoot
  double y=x[1];
  double tol_y=-1.0e-12;
  if ((y<0.0)&&(y>tol_y)) y=0.0;
  double phi=atan2(y,(x[0]-L_up));

  // A singular fct that solves the Laplace eqn
  u=1.0-pow(r,2.0/3.0)*sin(2.0/3.0*(phi+MathematicalConstants::Pi/2.0));
  double dudr=-2.0/3.0*pow(r,-1.0/3.0)*
   sin(2.0/3.0*(phi+MathematicalConstants::Pi/2.0));
  double dudphi=pow(r,2.0/3.0)*2.0/3.0*
   cos(2.0/3.0*(phi+MathematicalConstants::Pi/2.0));
  
  du_dx[0]=dudr*cos(phi)-1.0/r*dudphi*sin(phi);
  du_dx[1]=dudr*sin(phi)+1.0/r*dudphi*cos(phi);

 }


 /// \short "Singular" function
 double singular_fct(const Vector<double>& x)
 {
  double u=0.0;
  Vector<double> du_dx(2);
  singular_fct_and_gradient(x,u,du_dx);
  return u;
 }

 /// \short Gradient of "Singular" function
 Vector<double> gradient_of_singular_fct(const Vector<double>& x)
 {
  double u=0.0;
  Vector<double> du_dx(2);
  singular_fct_and_gradient(x,u,du_dx);
  return du_dx;
 }

 /// Exact solution
 double u_exact(const Vector<double>& x)
 {
  double u=0.0;
  if (!CommandLineArgs::command_line_flag_has_been_set
       ("--suppress_sing_in_exact_soln"))
   {
    u+=singular_fct(x);
   }
  if (CommandLineArgs::command_line_flag_has_been_set
       ("--add_sin_cosh_to_exact_soln"))
   {
    u+=sin(x[0])*cosh(x[1]);
   }
  return u;
 }

 /// Exact solution as vector
 void u_exact_as_vector(const Vector<double>& x, Vector<double>& u)
 {
  u[0]=u_exact(x);
 }

} // end_of_namespace



////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////



//==start_of_problem_class============================================
/// Poisson in backward-facing step. Dirichlet boundary conditions along
/// all boundaries.
//====================================================================
template<class ELEMENT>
class StepProblem : public Problem
{

public:


  /// Constructor
  StepProblem();

  /// Destructor 
 ~StepProblem()
  {
   // hierher: at some point delete things properly and do memory leak
   // check
   delete Bulk_mesh_pt->spatial_error_estimator_pt();
   delete Bulk_mesh_pt;
  }
 
 /// Update the after solve (empty)
 void actions_after_newton_solve(){}
 
 /// \short Update the problem specs before solve (empty)
 void actions_before_newton_solve() {}
 
 // Perform actions after mesh adaptation
 void actions_after_adapt()
  {
   // Recreate face elements
   create_face_elements();
   
   // Complete problem setup
   complete_problem_setup();

   // Rebuild global mesh
   rebuild_global_mesh();
  
  }
 
 /// Perform actions after mesh adaptation (empty)
 void actions_before_adapt()
  {
   // Kill face elements
   delete_face_elements();

   // Rebuild global mesh
   rebuild_global_mesh();
  }
 
 /// Access function for the specific mesh
 RefineableTriangleMesh<ELEMENT>* mesh_pt()
  {
   return dynamic_cast<RefineableTriangleMesh<ELEMENT>*>(Problem::mesh_pt());
  }
 
 /// Doc the solution
 void doc_solution(DocInfo& doc_info);
 
private:

 /// Do what it says
 void complete_problem_setup();
 
 /// Helper function to apply boundary conditions
 void apply_boundary_conditions();

 /// hierher Delete face elements and flush meshes
 void delete_face_elements()
  {
   if (CommandLineArgs::command_line_flag_has_been_set
       ("--dont_use_singularity"))
    {
     return;
    }

   // Loop over the flux jump elements
   unsigned n_element = Face_mesh_for_flux_jump_pt->nelement();
   for(unsigned e=0;e<n_element;e++)
    {
     delete Face_mesh_for_flux_jump_pt->element_pt(e);
    }
   
   // hierher: actually kill nodes too because they've been duplicated

   // Wipe the mesh
   Face_mesh_for_flux_jump_pt->flush_element_and_node_storage();

   // Loop over the bc elements
   n_element = Face_mesh_for_bc_pt->nelement();
   for(unsigned e=0;e<n_element;e++)
    {
     // Kill
     delete Face_mesh_for_bc_pt->element_pt(e);
    }
   
   // Wipe the mesh
   Face_mesh_for_bc_pt->flush_element_and_node_storage();

  }

 /// Create face elements
 void create_face_elements()
  {
   
   if (CommandLineArgs::command_line_flag_has_been_set
       ("--dont_use_singularity"))
    {
     return;
    }

   /// Map keeps a running count of duplicate nodes already created;
   /// existing_duplicate_node_pt[orig_node_pt]=new_node_pt.
   std::map<Node*,Node*> existing_duplicate_node_pt;

   // Flux jump elements on boundary 6 of region 1:
   //----------------------------------------------
   // NOTE: Since these duplicate nodes, these elements must be
   //----------------------------------------------------------
   //       constructed first!
   //       ------------------
   {
    // Where are we? Inside region 1 on boundary 6
    unsigned b=6;
    unsigned region_id=1;
    unsigned nel=Bulk_mesh_pt->nboundary_element_in_region(b,region_id);
    for (unsigned e=0;e<nel;e++)
     {
      FiniteElement* el_pt=
       Bulk_mesh_pt->boundary_element_in_region_pt(b,region_id,e);
      
      // What is the index of the face of the bulk element at the boundary
      int face_index = Bulk_mesh_pt->
       face_index_at_boundary_in_region(b,region_id,e);
      
      // Build the corresponding flux jump element
      PoissonWithSingularityFluxJumpFaceElement<ELEMENT>* flux_jump_element_pt 
       = new PoissonWithSingularityFluxJumpFaceElement<ELEMENT>
       (el_pt,face_index,existing_duplicate_node_pt,Flux_jump_el_id);

      //Add the flux jump element to the mesh
      Face_mesh_for_flux_jump_pt->add_element_pt(flux_jump_element_pt);
     }
   }
   
   // Now add all new (duplicated) nodes to mesh
   for (std::map<Node*,Node*>::iterator it=
         existing_duplicate_node_pt.begin();
        it!=existing_duplicate_node_pt.end();it++)
    {
     Face_mesh_for_flux_jump_pt->add_node_pt((*it).second);
    }
   

   // Now loop over bulk elements in region 1 ("torus" around singularity)
   //---------------------------------------------------------------------
   // and swap over any of their nodes have been replaced
   //----------------------------------------------------
   unsigned region_id=1;
   unsigned n_el=Bulk_mesh_pt->nregion_element(region_id);
   for (unsigned e=0;e<n_el;e++)
    {
     ELEMENT* bulk_el_pt=dynamic_cast<ELEMENT*>(
      Bulk_mesh_pt->region_element_pt(region_id,e));
   
     // Loop over all nodes and check if they're amongst the replaced
     // ones
     unsigned nnod=bulk_el_pt->nnode();
     for (unsigned j=0;j<nnod;j++)
      {
       Node* nod_pt=bulk_el_pt->node_pt(j);
     
       // Find this original node in the map; if we find it
       // it's already been duplicated
       std::map<Node*,Node*>::iterator it=
        existing_duplicate_node_pt.find(nod_pt);
       if (it!=existing_duplicate_node_pt.end())
        {
         // Use the existing duplicate node
         bulk_el_pt->node_pt(j)=(*it).second;
        }
      }   
    }


   // BC elements live on boundary 4:
   //--------------------------------
   {
    unsigned b=4;
    unsigned n_element = Bulk_mesh_pt->nboundary_element(b);
    
    // Loop over the bulk elements adjacent to boundary b
    for(unsigned e=0;e<n_element;e++)
     {
      // Get pointer to the bulk element that is adjacent to boundary b
      ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
       Bulk_mesh_pt->boundary_element_pt(b,e));

      //Find the index of the face of element e along boundary b 
      int face_index = Bulk_mesh_pt->face_index_at_boundary(b,e);
      
      // Build the corresponding bc element
      PoissonWithSingularityBCFaceElement<ELEMENT>* bc_element_pt = new 
       PoissonWithSingularityBCFaceElement<ELEMENT>(bulk_elem_pt,face_index,
                                                    BC_el_id);
      
      //Add the bc element to the surface mesh
      Face_mesh_for_bc_pt->add_element_pt(bc_element_pt);
      
     }
   }

  } 


 /// Pointer to the bulk mesh
 RefineableTriangleMesh<ELEMENT> *Bulk_mesh_pt;
 
 /// Face element mesh for jump in interior of domain
 Mesh* Face_mesh_for_flux_jump_pt;

 /// Face element mesh for BC/singularity
 Mesh* Face_mesh_for_bc_pt;

 /// Mesh for (single) element containing singular fct
 Mesh* Singular_fct_element_mesh_pt;

 /// \short Enumeration for IDs of FaceElements (used to figure out
 /// who's added what additional nodal data...)
 enum{Flux_jump_el_id,BC_el_id};
 
}; // end_of_problem_class


//==start_of_constructor==================================================
/// Constructor for StepProblem problem
//========================================================================
template<class ELEMENT>
StepProblem<ELEMENT>::StepProblem()
{

  // Build the mesh
  double uniform_element_area=0.1;
  Bulk_mesh_pt=Global_Physical_Variables::build_the_mesh<ELEMENT>
   (uniform_element_area);

  // Let's have a look at the boundary enumeration
  Bulk_mesh_pt->output_boundaries("boundaries.dat");

  // Set error estimator for bulk mesh
  Z2ErrorEstimator* error_estimator_pt=new Z2ErrorEstimator;
  Bulk_mesh_pt->spatial_error_estimator_pt()=error_estimator_pt;
  
  // Set element size limits
  Bulk_mesh_pt->max_element_size()=0.1;
  Bulk_mesh_pt->min_element_size()=1e-30;
  Bulk_mesh_pt->max_permitted_error()=0.005;
  Bulk_mesh_pt->min_permitted_error()=0.0;
  
  
  // Add sub-mesh
  add_sub_mesh(Bulk_mesh_pt);


  if (!CommandLineArgs::command_line_flag_has_been_set
      ("--dont_use_singularity"))
   {
    
    // Create element that stores the singular fct and its amplitude
    //---------------------------------------------------------------
    ScalableSingularityForPoissonElement<ELEMENT>* el_pt=
     new ScalableSingularityForPoissonElement<ELEMENT>;
    
    // Pass fct pointers:
    el_pt->unscaled_singular_fct_pt()
     =&Global_Physical_Variables::singular_fct;
    el_pt->gradient_of_unscaled_singular_fct_pt()=
     &Global_Physical_Variables::gradient_of_singular_fct;
    
    // Add to mesh
    Singular_fct_element_mesh_pt=new Mesh;
    Singular_fct_element_mesh_pt->add_element_pt(el_pt);
    add_sub_mesh(Singular_fct_element_mesh_pt);
    
    
    // Create face elements for impositition of BC and flux jump
    //----------------------------------------------------------
    Face_mesh_for_bc_pt=new Mesh;
    Face_mesh_for_flux_jump_pt=new Mesh;
    create_face_elements();
    
    // Add to mesh
    add_sub_mesh(Face_mesh_for_bc_pt);
    add_sub_mesh(Face_mesh_for_flux_jump_pt);
    
   }

  // Build global mesh
  build_global_mesh();
  
  // Complete problem setup
  complete_problem_setup();
    
  // Setup equation numbering scheme
  oomph_info <<"Number of equations: " 
             << this->assign_eqn_numbers() 
             << std::endl;

} // end_of_constructor


//==start_of_complete======================================================
/// Set boundary condition, and complete the build of
/// all elements
//========================================================================
template<class ELEMENT>
void StepProblem<ELEMENT>::complete_problem_setup()
{
 
 if (!CommandLineArgs::command_line_flag_has_been_set
       ("--dont_use_singularity"))
  {
   
   // Loop over the elements to set up element-specific
   // things that cannot be handled by constructor
   
   // Bulk elements in region 1
   unsigned region_id=1;
   unsigned n_el=Bulk_mesh_pt->nregion_element(region_id);
   for (unsigned e=0;e<n_el;e++)
    {
     ELEMENT* bulk_el_pt=dynamic_cast<ELEMENT*>(
      Bulk_mesh_pt->region_element_pt(region_id,e));
     
     // Tell the bulk element about the singular fct
     bulk_el_pt->poisson_sing_el_pt()=
      dynamic_cast<TemplateFreeScalableSingularityForPoissonElement*>(
       Singular_fct_element_mesh_pt->element_pt(0));
    }
   
   
   // Flux jump elements
   unsigned n_element = Face_mesh_for_flux_jump_pt->nelement();
   for(unsigned e=0;e<n_element;e++)
    {
     // Upcast from GeneralisedElement to the present element
     PoissonWithSingularityFluxJumpFaceElement<ELEMENT>* el_pt = dynamic_cast<
      PoissonWithSingularityFluxJumpFaceElement<ELEMENT>*>(
       Face_mesh_for_flux_jump_pt->element_pt(e));
     
     // Tell the element about the singular fct
     el_pt->set_poisson_sing_el_pt(
      dynamic_cast<ScalableSingularityForPoissonElement<ELEMENT>*>(
       Singular_fct_element_mesh_pt->element_pt(0)));
    }
   
   
   
   // BC elements
   n_element =  Face_mesh_for_bc_pt->nelement();
   for(unsigned e=0;e<n_element;e++)
    {
     // Upcast from GeneralisedElement to the present element
     PoissonWithSingularityBCFaceElement<ELEMENT>* el_pt = dynamic_cast<
      PoissonWithSingularityBCFaceElement<ELEMENT>*>(
       Face_mesh_for_bc_pt->element_pt(e));
     
     // Tell the element about the singular fct
     el_pt->set_poisson_sing_el_pt(
      dynamic_cast<ScalableSingularityForPoissonElement<ELEMENT>*>(
       Singular_fct_element_mesh_pt->element_pt(0)));
    }
   
   
   // Find regularised bulk element
   ELEMENT* regularised_bulk_el_pt=0;
   Vector<double> s_reg(2);
   
   // We're looking for a bulk element that has a node on the
   // the corner of the backward step; all candidate elements live on
   // boundary 4
   unsigned b=4;
   n_element = Bulk_mesh_pt->nboundary_element(b);
   
   // Loop over the bulk elements adjacent to boundary b
   for(unsigned e=0;e<n_element;e++)
    {
     // Get pointer to the bulk element that is adjacent to boundary b
     ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
      Bulk_mesh_pt->boundary_element_pt(b,e));
     
     unsigned nnod=bulk_elem_pt->nnode();
     for (unsigned j=0;j<nnod;j++)
      {
       Node* nod_pt=bulk_elem_pt->node_pt(j);
       
       // Are we at the corner?
       double distance=sqrt((nod_pt->x(0)-Global_Physical_Variables::L_up)*
                            (nod_pt->x(0)-Global_Physical_Variables::L_up)+
                            (nod_pt->x(1)-0.0)*(nod_pt->x(1)-0.0));
       double tol=1.0e-10;
       if (distance<tol)
        {
         
         oomph_info << "hierher: better check that some of its nodes are"
                    << " above and some below the normal\n";
         
         regularised_bulk_el_pt=bulk_elem_pt;
         bulk_elem_pt->local_coordinate_of_node(j,s_reg);
         oomph_info 
          << "Found regularised bulk element; s_reg ="
          << s_reg[0] << " " << s_reg[1] << " coords: "
          << regularised_bulk_el_pt->interpolated_x(s_reg,0) << " "
          << regularised_bulk_el_pt->interpolated_x(s_reg,1) << " "
          << std::endl;
         break;
        }
      }
     if (regularised_bulk_el_pt!=0) break;
    }
   
   // Set pointer to bulk element and local coordinate
   // in it to identify where regularity (zero flux in specified
   // direction) is imposed on the FE solution. Needs to be done
   // after the node pointers have been updated otherwise
   // the classification of external data goes wrong!
   Vector<double> n(2);
   n[0]=1.0/sqrt(2.0);
   n[1]=1.0/sqrt(2.0);
   dynamic_cast<ScalableSingularityForPoissonElement<ELEMENT>*>(
    Singular_fct_element_mesh_pt->element_pt(0))->
    specify_regularisation(regularised_bulk_el_pt,s_reg,n);
   
  }
 
 // Apply bcs
 apply_boundary_conditions();
 
}

//==start_of_apply_bc=====================================================
/// Helper function to apply boundary conditions
//========================================================================
template<class ELEMENT>
void StepProblem<ELEMENT>::apply_boundary_conditions()
{
 
 // Set the boundary conditions for this problem: All nodes are
 // free by default -- just pin the ones that have Dirichlet conditions
 // here.
 unsigned num_bound = Bulk_mesh_pt->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {   
   // Leave internal boundary alone
   if (ibound!=6)
    { 
     unsigned num_nod=Bulk_mesh_pt->nboundary_node(ibound);
     for (unsigned inod=0;inod<num_nod;inod++)
      {
       Bulk_mesh_pt->boundary_node_pt(ibound,inod)->pin(0);
      }
    }
  } // end loop over boundaries
 
 
 // Now set boundary values
 for (unsigned ibound=0;ibound<num_bound;ibound++)
  {
   // Leave internal boundary alone
   if (ibound!=6)
    {
     unsigned num_nod= Bulk_mesh_pt->nboundary_node(ibound);
     for (unsigned inod=0;inod<num_nod;inod++)
      {
       Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(ibound,inod);
       Vector<double> x(2);
       x[0]=nod_pt->x(0);
       x[1]=nod_pt->x(1);
       double u=Global_Physical_Variables::u_exact(x);
       nod_pt->set_value(0,u);
      }
    }
  }

 if (!CommandLineArgs::command_line_flag_has_been_set
       ("--dont_use_singularity"))
  {
   // Now unpin nodal values where the bc conditions are enforceed
   // by Lagrange multiplier to ensure that the sum of fe and singular
   // solution is correct
   unsigned nel=Face_mesh_for_bc_pt->nelement();
   for (unsigned e=0;e<nel;e++)
    {
     // Get element
     PoissonWithSingularityBCFaceElement<ELEMENT>* el_pt=
      dynamic_cast<PoissonWithSingularityBCFaceElement<ELEMENT>*>(
       Face_mesh_for_bc_pt->element_pt(e));
     
     // Specify desired nodal values for compound solution
     unsigned nnod=el_pt->nnode();
     Vector<double> nodal_boundary_value(nnod,0.0);

     // Unpin the fe part of the solution
     for (unsigned j=0;j<nnod;j++)
      {
       el_pt->unpin_u_fe_at_specified_local_node(j);
       
       // Now here's another subtle one! If this node stores three values
       // they come from: original Poisson dof; one Lagrange multiplier
       // to ensure continuity of solution across "jump"; one Lagrange
       // multiplier from enforcing Dirichlet boundary conditions. 
       // The latter two enforce the same thing, so we can only apply
       // the constraint once --> Pin the Lagrange multiplier that tries to 
       // enforce the Dirichlet BC
       unsigned nval=el_pt->node_pt(j)->nvalue();
       if (nval==3)
        {
         el_pt->pin_lagrange_multiplier_at_specified_local_node(j);
        }
       
       Node* nod_pt=el_pt->node_pt(j);
       Vector<double> x(2);
       x[0]=nod_pt->x(0);
       x[1]=nod_pt->x(1);
       double u=Global_Physical_Variables::u_exact(x);
       nodal_boundary_value[j]=u;
      }
     
     // Tell the element about the desired nodal boundary values
     el_pt->set_nodal_boundary_values(nodal_boundary_value);
     
     
    }
  }

} // end set bc


//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void StepProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
{
  ofstream some_file;
  char filename[100];

  // Number of plot points
  unsigned npts=10;

  // Output solution
  sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
  doc_info.number());
  some_file.open(filename);
  Bulk_mesh_pt->output(some_file,npts);
  some_file.close();


  // Output solution just using vertices so we can see the mesh
  sprintf(filename,"%s/coarse_soln%i.dat",doc_info.directory().c_str(),
  doc_info.number());
  some_file.open(filename);
  npts=2;
  Bulk_mesh_pt->output(some_file,npts);
  some_file.close();
  

  // Plot "extended solution" showing contributions; also work out
  // average element size
  double av_el_size=0.0;
  sprintf(filename,"%s/extended_soln%i.dat",doc_info.directory().c_str(),
          doc_info.number());
  some_file.open(filename);
  unsigned nel=Bulk_mesh_pt->nelement();
  for (unsigned e=0;e<nel;e++)
   {
    npts=20;
    ELEMENT* el_pt=
     dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));
    el_pt->output_with_various_contributions(some_file,npts);
    av_el_size+=el_pt->size();
  }
  some_file.close();
  av_el_size/=double(nel);


  // Get error
  double error,norm; 
  sprintf(filename,"%s/error%i.dat",doc_info.directory().c_str(),
          doc_info.number());
  some_file.open(filename);
  Bulk_mesh_pt->compute_error(some_file,
                              Global_Physical_Variables::u_exact_as_vector,
                              error,norm);
  some_file.close();

  // Doc error norm:
  oomph_info << "\n av el size, av h, Ndof, # bulk els, Norm of error    : "   
             << av_el_size << " " 
             << sqrt(av_el_size) << " " 
             << ndof() << " " 
             << Bulk_mesh_pt->nelement() << " " 
             << sqrt(error) << std::endl;
  oomph_info << "Norm of solution : " << sqrt(norm) << std::endl << std::endl;
  oomph_info << std::endl;
  
  // Exact solution
  sprintf(filename,"%s/exact_soln%i.dat",doc_info.directory().c_str(),
          doc_info.number());
  some_file.open(filename);
  unsigned nplot=5;
  Bulk_mesh_pt->output_fct(some_file,nplot,
                           Global_Physical_Variables::u_exact_as_vector);
  some_file.close();
  

  if (!CommandLineArgs::command_line_flag_has_been_set
      ("--dont_use_singularity"))
   {
    oomph_info 
     << "Amplitude of singular function: "
     << dynamic_cast<TemplateFreeScalableSingularityForPoissonElement*>(
      Singular_fct_element_mesh_pt->element_pt(0))->
     amplitude_of_singular_fct() << std::endl;
   }

} // end_of_doc_solution



////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


//==start_of_main======================================================
/// Driver for backward step with impedance outflow bc
//=====================================================================
int main(int argc, char **argv)
{

 // Store command line arguments
 CommandLineArgs::setup(argc,argv);
  
 // Don't subtract off singularity
 CommandLineArgs::specify_command_line_flag(
  "--dont_use_singularity");

 // Use sin cosh in exact solution
 CommandLineArgs::specify_command_line_flag(
  "--add_sin_cosh_to_exact_soln");

 // Suppress singular term in exact solution
 CommandLineArgs::specify_command_line_flag(
  "--suppress_sing_in_exact_soln");

 // Parse command line
 CommandLineArgs::parse_and_assign(); 
 
 // Doc what has actually been specified on the command line
 CommandLineArgs::doc_specified_flags();

 
 // Set up doc info
 // ---------------
 
 // Label for output
 DocInfo doc_info;
 
 // Set output directory
 doc_info.set_directory("RESLT");
 
 // Step number
 doc_info.number()=0;
 
  // Build the problem with 
 StepProblem<ProjectablePoissonElement<MyTPoissonElement<2,3> > >
  problem;
 
 // Solve, refine uniformly and keep going
  unsigned max_adapt = 5;
  for (unsigned i=0;i<max_adapt;i++)
   {
    // Solve the bloody thing
    problem.newton_solve();
    problem.doc_solution(doc_info);
    doc_info.number()++;
    if (i!=(max_adapt-1)) problem.refine_uniformly();
   }

 
}

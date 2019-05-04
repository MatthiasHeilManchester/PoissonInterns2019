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

#include<fenv.h>

//Generic routines
#include "generic.h"

// Poisson
#include "poisson.h"

// Get the mesh
#include "meshes/tetgen_mesh.h" 

// The mesh
#include "meshes/triangle_mesh.h"
 
// The mesh
#include "gmsh_tet_mesh.h"

#include "poisson_sing_face_element.h"

using namespace std;

using namespace oomph;



///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////



//=========================================================================
/// \short "Disk" in 3d: zeta[0]=x; zeta[1]=y. (Parametrisation must not 
/// have coordinate singularities!)
//=========================================================================
class Disk : public GeomObject
{

public:

 /// Constructor
 Disk(const double& epsilon, const double& z_offset=0.0) :
  GeomObject(2,3), Epsilon(epsilon),  N(3), Z_offset(z_offset)
  {}

 /// Broken copy constructor
 Disk(const Disk& dummy) 
  { 
   BrokenCopy::broken_copy("Disk");
  } 
 
 /// Broken assignment operator
 void operator=(const Disk&) 
  {
   BrokenCopy::broken_assign("Disk");
  }

 /// Destructor
 virtual ~Disk()
  {}

 /// Access fct to amplitude of disk warping
 double& epsilon()
  {
   return Epsilon;
  }

 /// \short Position Vector at Lagrangian coordinate zeta 
 void position(const Vector<double>& zeta, Vector<double>& r) const
  {
   // Position Vector
   r[0] = zeta[0];
   r[1] = zeta[1];
   double radius=sqrt(r[0]*r[0]+r[1]*r[1]);
   double phi=atan2(r[1],r[0]);
   r[2]=Z_offset+w(radius,phi);
  }


 /// \short Position Vector at Lagrangian coordinate zeta 
 void outer_normal_on_boundary(const double& phi, 
                               Vector<double>& r,
                               Vector<double>& n,
                               Vector<double>& t) const
  {
   // Position Vector
   r[0] = cos(phi);
   r[1] = sin(phi);
   r[2]=0.0+w(1.0,phi);

   Vector<double> dr_dr(3);
   dr_dr[0]=cos(phi);
   dr_dr[1]=sin(phi);
   dr_dr[2]=dwdr(1.0,phi);
   double inv_norm=1.0/sqrt(dr_dr[0]*dr_dr[0]+
                            dr_dr[1]*dr_dr[1]+
                            dr_dr[2]*dr_dr[2]);

   n[0]=dr_dr[0]*inv_norm;
   n[1]=dr_dr[1]*inv_norm;
   n[2]=dr_dr[2]*inv_norm;

   Vector<double> dr_dphi(3);
   dr_dphi[0]=-sin(phi);
   dr_dphi[1]=cos(phi);
   dr_dphi[2]=dwdphi(1.0,phi);

   inv_norm=1.0/sqrt(dr_dphi[0]*dr_dphi[0]+
                     dr_dphi[1]*dr_dphi[1]+
                     dr_dphi[2]*dr_dphi[2]);
   
   t[0]=dr_dphi[0]*inv_norm;
   t[1]=dr_dphi[1]*inv_norm;
   t[2]=dr_dphi[2]*inv_norm;
  }

private:

 /// Vertical deflection
 double w(const double& r, const double& phi) const
  {
   return Epsilon*cos(double(N)*phi)*pow(r,2);
  }

 /// Deriv of vertical deflection w.r.t. radius
 double dwdr(const double& r, const double& phi) const
  {
   return Epsilon*cos(double(N)*phi)*2.0*r;
  }

 /// Deriv of vertical deflection w.r.t. angle
 double dwdphi(const double& r, const double& phi) const
  {
   return -Epsilon*double(N)*sin(double(N)*phi)*pow(r,2);
  }

 /// Amplitude of non-axisymmetric deformation
 double Epsilon;

 /// Wavenumber of non-axisymmetric deformation
 unsigned N;

 /// Vertical offset
 double Z_offset;

};



///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////


//=============================================================
/// Namespace for problem parameters
//=============================================================
namespace Global_Parameters
{
 /// (Half-)width of the box
 double Box_half_width = 1.5;

 /// (Half)height of the box
 double Box_half_length = 1.0;

 /// Amplitude of disk warping
 double Epsilon=0.0;

 /// GeomObject representing the disk:
 Disk* Disk_pt=new Disk(Epsilon);
 
 /// Asymptotic solution near perimeter of disk
 double asymptotic_solution(const double& rho, 
                            const double& phi)
 {
  return 1.0-(2.0/MathematicalConstants::Pi)*sqrt((cos(phi)+1.0)*rho); 
 }

 /// \short "Singular" function and gradient
 void singular_fct_and_gradient(const Vector<double>& x,
                                double& u, Vector<double>& du_dx)
 {
  // inverse transformation 
  // https://en.wikipedia.org/wiki/Oblate_spheroidal_coordinates#Inverse_transformation
  double radius=1.0;
  double rho=sqrt(x[0]*x[0]+x[1]*x[1]);
  double d1=sqrt((rho+radius)*(rho+radius)+x[2]*x[2]);
  double d2=sqrt((rho-radius)*(rho-radius)+x[2]*x[2]);
  double mu=acosh((d1+d2)/(2.0*radius));
  //double nu=acos ((d1-d2)/(2.0*radius));

  double eta=sinh(mu);
  //double xi = sin(nu);

  // https://colalg.math.csusb.edu/~devel/IT/main/m06_inverse/src/s02_tanflip.html
  double arccot_eta=MathematicalConstants::Pi/2.0-atan(eta);
  u=2.0/MathematicalConstants::Pi*arccot_eta; 


  double X=x[0];
  double Y=x[1];
  double Z=x[2];
  double MapleGenVar1=0.0;
  double MapleGenVar2=0.0;
  double t0=0.0;

      MapleGenVar1 = -2.0*X*(sqrt(X*X+Y*Y+2.0*sqrt(X*X+Y*Y)*radius+Z*Z+radius*
radius)*sqrt(X*X+Y*Y)+sqrt(X*X+Y*Y-2.0*sqrt(X*X+Y*Y)*radius+Z*Z+radius*radius)*
sqrt(X*X+Y*Y)-sqrt(X*X+Y*Y+2.0*sqrt(X*X+Y*Y)*radius+Z*Z+radius*radius)*radius+
sqrt(X*X+Y*Y-2.0*sqrt(X*X+Y*Y)*radius+Z*Z+radius*radius)*radius)*(sqrt(X*X+Y*Y+
2.0*sqrt(X*X+Y*Y)*radius+Z*Z+radius*radius)+sqrt(X*X+Y*Y-2.0*sqrt(X*X+Y*Y)*
radius+Z*Z+radius*radius))/sqrt((sqrt(X*X+Y*Y+2.0*sqrt(X*X+Y*Y)*radius+Z*Z+
radius*radius)+sqrt(X*X+Y*Y-2.0*sqrt(X*X+Y*Y)*radius+Z*Z+radius*radius)-2.0*
radius)/radius);
      MapleGenVar2 = 1/(sqrt(X*X+Y*Y+2.0*sqrt(X*X+Y*Y)*radius+Z*Z+radius*radius
))/sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y-2.0*sqrt(X*X+Y*Y)*radius+Z*Z+radius*radius)/sqrt(
(sqrt(X*X+Y*Y+2.0*sqrt(X*X+Y*Y)*radius+Z*Z+radius*radius)+sqrt(X*X+Y*Y-2.0*sqrt
(X*X+Y*Y)*radius+Z*Z+radius*radius)+2.0*radius)/radius)/0.3141592653589793E1/(
sqrt(X*X+Y*Y+2.0*sqrt(X*X+Y*Y)*radius+Z*Z+radius*radius)*sqrt(X*X+Y*Y-2.0*sqrt(
X*X+Y*Y)*radius+Z*Z+radius*radius)+X*X+Y*Y+Z*Z+radius*radius);
      t0 = MapleGenVar1*MapleGenVar2;

  du_dx[0]=t0;



      MapleGenVar1 = -2.0*Y*(sqrt(X*X+Y*Y+2.0*sqrt(X*X+Y*Y)*radius+Z*Z+radius*
radius)*sqrt(X*X+Y*Y)+sqrt(X*X+Y*Y-2.0*sqrt(X*X+Y*Y)*radius+Z*Z+radius*radius)*
sqrt(X*X+Y*Y)-sqrt(X*X+Y*Y+2.0*sqrt(X*X+Y*Y)*radius+Z*Z+radius*radius)*radius+
sqrt(X*X+Y*Y-2.0*sqrt(X*X+Y*Y)*radius+Z*Z+radius*radius)*radius)*(sqrt(X*X+Y*Y+
2.0*sqrt(X*X+Y*Y)*radius+Z*Z+radius*radius)+sqrt(X*X+Y*Y-2.0*sqrt(X*X+Y*Y)*
radius+Z*Z+radius*radius))/sqrt((sqrt(X*X+Y*Y+2.0*sqrt(X*X+Y*Y)*radius+Z*Z+
radius*radius)+sqrt(X*X+Y*Y-2.0*sqrt(X*X+Y*Y)*radius+Z*Z+radius*radius)-2.0*
radius)/radius);
      MapleGenVar2 = 1/(sqrt(X*X+Y*Y+2.0*sqrt(X*X+Y*Y)*radius+Z*Z+radius*radius
))/sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y-2.0*sqrt(X*X+Y*Y)*radius+Z*Z+radius*radius)/sqrt(
(sqrt(X*X+Y*Y+2.0*sqrt(X*X+Y*Y)*radius+Z*Z+radius*radius)+sqrt(X*X+Y*Y-2.0*sqrt
(X*X+Y*Y)*radius+Z*Z+radius*radius)+2.0*radius)/radius)/0.3141592653589793E1/(
sqrt(X*X+Y*Y+2.0*sqrt(X*X+Y*Y)*radius+Z*Z+radius*radius)*sqrt(X*X+Y*Y-2.0*sqrt(
X*X+Y*Y)*radius+Z*Z+radius*radius)+X*X+Y*Y+Z*Z+radius*radius);
      t0 = MapleGenVar1*MapleGenVar2;

  du_dx[1]=t0;



      t0 = -4.0*Z/0.3141592653589793E1/sqrt((sqrt(X*X+Y*Y+2.0*sqrt(X*X+Y*Y)*
radius+Z*Z+radius*radius)+sqrt(X*X+Y*Y-2.0*sqrt(X*X+Y*Y)*radius+Z*Z+radius*
radius)-2.0*radius)/radius)/sqrt(X*X+Y*Y+2.0*sqrt(X*X+Y*Y)*radius+Z*Z+radius*
radius)/sqrt(X*X+Y*Y-2.0*sqrt(X*X+Y*Y)*radius+Z*Z+radius*radius)/sqrt((sqrt(X*X
+Y*Y+2.0*sqrt(X*X+Y*Y)*radius+Z*Z+radius*radius)+sqrt(X*X+Y*Y-2.0*sqrt(X*X+Y*Y)
*radius+Z*Z+radius*radius)+2.0*radius)/radius);


  du_dx[2]=t0;

  // du_dx[0]=0.0;
  // du_dx[1]=0.0;
  // du_dx[2]=0.0;

 }


 /// \short "Singular" function
 double singular_fct(const Vector<double>& x)
 {
  double u=0.0;
  Vector<double> du_dx(3);
  singular_fct_and_gradient(x,u,du_dx);
  return u;
 }

 /// \short Gradient of "Singular" function
 Vector<double> gradient_of_singular_fct(const Vector<double>& x)
 {
  double u=0.0;
  Vector<double> du_dx(3);
  singular_fct_and_gradient(x,u,du_dx);
  return du_dx;
 }


 /// Exact solution as a Vector: Gupta solution 
 void get_exact_u(const Vector<double>& x, Vector<double>& u)
 {
  u[0]=singular_fct(x);
 }

}

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////


#include "tetmesh_faceted_surfaces.h"

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////



//====================================================================
/// Disk in container
//====================================================================
template<class ELEMENT> 
class DiskInContainerProblem : public Problem
{

public:

 /// Constructor
 DiskInContainerProblem();
  
 /// Destructor (empty)
 ~DiskInContainerProblem()
  {
   //Delete the objects
   unsigned nh = Inner_boundary_pt.size();
   for(unsigned h=0;h<nh;++h)
    {
     delete Inner_boundary_pt[h];
    }
   delete Outer_boundary_pt;
  }

      
 /// Actions before adapt. Kill stuff
 void actions_before_adapt()
  {
   delete LV_pt;
   LV_pt=0;

   // Geom object representations need to be (re)built
   Geom_objects_are_out_of_date=true;

   // hierher
   // // Kill face elements
   // delete_face_elements();

   // // Rebuild global mesh
   // rebuild_global_mesh();
  }


 /// Totally new mesh; build elements and apply boundary conditions
 void actions_after_adapt()
  {
   // hierher
   // Kill face elements
   delete_face_elements();

   // Recreate face elements
   create_face_elements();

   // Complete problem setup
   complete_problem_setup();

   // Rebuild global mesh
   rebuild_global_mesh();
  }
 

 /// Build line visualiser
 void build_line_visualiser()
  {

   LV_pt=0;
   oomph_info << "hierher Not building line visualiser\n";
   return;

   // How many points to you want?
   unsigned int  npt=1000;
   Vector<Vector<double> > coord_vec(npt);
   for (unsigned j=0;j<npt;j++)
    {
     coord_vec[j].resize(3);
     double coord=-Global_Parameters::Box_half_width+
      2.0*Global_Parameters::Box_half_width*
      double(j)/double(npt-1);
     coord_vec[j][0]=coord;
     coord_vec[j][1]=coord;
     coord_vec[j][2]=1.0e-2;
    }

   // Setup line visualiser
   LV_pt=new LineVisualiser(Bulk_mesh_pt,
                            coord_vec);
   
  }



 /// Update the problem specs before solve: (empty)
 void actions_before_newton_solve(){}

 /// Update the problem specs before solve (empty)
 void actions_after_newton_solve(){}
 
 /// Doc the solution
 void doc_solution(const unsigned& nplot, DocInfo& doc_info);

private:
 
 /// Apply BCs and make elements functional
 void complete_problem_setup();

 /// Helper function to apply boundary conditions
 void apply_boundary_conditions();
 
 /// Setup disk on disk plots
 void setup_disk_on_disk_plots();


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
     // Kill
     delete Face_mesh_for_flux_jump_pt->element_pt(e);
    }
   
   // Kill (duplicated!) nodes
   unsigned nnod=Face_mesh_for_flux_jump_pt->nnode();
   for (unsigned j=0;j<nnod;j++)
    {
     delete Face_mesh_for_flux_jump_pt->node_pt(j);
    }

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


   // Map to keep track of mapping between old and duplicated nodes
   std::map<Node*,Node*> existing_duplicate_node_pt;

   // Flux jump elements on boundary of torus
   //----------------------------------------
   // NOTE: Since these duplicate nodes, these elements must be
   //----------------------------------------------------------
   //       constructed first!
   //       ------------------
   {
    // hierher
    ofstream some_file;
    some_file.open("flux_jump_elements.dat");

    // Where are we?
    unsigned region_id=One_based_torus_region_id-1;
    for (unsigned b=First_boundary_for_torus_id;
         b<=Last_boundary_for_torus_id;b++)
     {
      unsigned nel=Bulk_mesh_pt->nboundary_element_in_region(b,region_id);
      for (unsigned e=0;e<nel;e++)
       {
        FiniteElement* el_pt=
         Bulk_mesh_pt->boundary_element_in_region_pt(b,region_id,e);
        
        // What is the index of the face of the bulk element at the boundary
        int face_index = Bulk_mesh_pt->
         face_index_at_boundary_in_region(b,region_id,e);
        
        // Build the corresponding flux jump element
        PoissonWithSingularityFluxJumpFaceElement<ELEMENT>* 
         flux_jump_element_pt 
         = new PoissonWithSingularityFluxJumpFaceElement<ELEMENT>
         (el_pt,face_index,existing_duplicate_node_pt,Flux_jump_el_id);
        
        //Add the flux jump element to the mesh
        Face_mesh_for_flux_jump_pt->add_element_pt(flux_jump_element_pt);

        // hierher
        flux_jump_element_pt->output(some_file);
        
       }
     }

    // hierher
    some_file.close();
    
   }
   

   {
    ofstream some_file;
    some_file.open("bulk_nodes.dat");
    unsigned nnod=Bulk_mesh_pt->nnode();
    for (unsigned j=0;j<nnod;j++)
     {
      Node* nod_pt=Bulk_mesh_pt->node_pt(j);
      some_file << nod_pt << " " 
                << nod_pt->x(0) << " "  
                << nod_pt->x(1) << " "  
                << nod_pt->x(2) << " "  
                << std::endl;
     }
    some_file.close();
    some_file.open("bulk_nodes_from_elements.dat");
    unsigned nel=Bulk_mesh_pt->nelement();
    for (unsigned e=0;e<nel;e++)
     {
      FiniteElement* el_pt=Bulk_mesh_pt->finite_element_pt(e);
      unsigned nnod=el_pt->nnode();
      for (unsigned j=0;j<nnod;j++)
       {
        Node* nod_pt=el_pt->node_pt(j);
        some_file<< nod_pt << " " 
                 << nod_pt->x(0) << " "  
                 << nod_pt->x(1) << " "  
                 << nod_pt->x(2) << " "  
                 << std::endl;
       }
     }
    some_file.close();
   }

   // Now add all duplicated nodes to mesh
   ofstream some_file;
   some_file.open("duplicated_nodes.dat");
   for (std::map<Node*,Node*>::iterator it=
         existing_duplicate_node_pt.begin();
        it!=existing_duplicate_node_pt.end();it++)
    {
     Face_mesh_for_flux_jump_pt->add_node_pt((*it).second);
     some_file << (*it).second->x(0) << " " 
               << (*it).second->x(1) << " " 
               << (*it).second->x(2) << " "
               << std::endl;
    }
   some_file.close();
   //pause("done duplicated_nodes.dat");
   
   
   // Now loop over bulk elements in torus region ("torus" around singularity)
   //-------------------------------------------------------------------------
   // and swap over any of their nodes that have been replaced
   //---------------------------------------------------------
   unsigned region_id=One_based_torus_region_id-1;
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
  


   // BC elements live on disk inside torus
   //--------------------------------------
   {
    // hierher
    ofstream some_file;
    some_file.open("bc_elements.dat");

    unsigned n=One_based_boundary_id_for_disk_within_torus.size();
    for (unsigned i=0;i<n;i++)
     {
      unsigned b=One_based_boundary_id_for_disk_within_torus[i]-1;
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

        // hierher
        bc_element_pt->output(some_file);
        
       }
     }
    
    // hierher
    some_file.close();
   }

  }

 /// Bulk mesh
 RefineableGmshTetMesh<ELEMENT>* Bulk_mesh_pt;

 /// Face element mesh for jump in interior of domain
 Mesh* Face_mesh_for_flux_jump_pt;

 /// Face element mesh for BC/singularity
 Mesh* Face_mesh_for_bc_pt;

 /// Mesh for (single) element containing singular fct
 Mesh* Singular_fct_element_mesh_pt;

 /// \short Enumeration for IDs of FaceElements (used to figure out
 /// who's added what additional nodal data...)
 enum{bla_hierher,Flux_jump_el_id,BC_el_id};


 /// Flag to force update on geom object representations 
 bool Geom_objects_are_out_of_date;

 /// Storage for the outer boundary object
 TetMeshFacetedClosedSurface* Outer_boundary_pt;

 /// Inner boundary
 Vector<TetMeshFacetedSurface*> Inner_boundary_pt;

 /// First boundary ID for outer boundary
 unsigned First_boundary_id_for_outer_boundary;

 /// Boundary ID: top
 unsigned Top_boundary_id;

 /// Boundary ID: bottom
 unsigned Bottom_boundary_id;

 /// Boundary ID: first for sheet
 unsigned First_sheet_boundary_id;

 /// Boundary ID: last for sheet
 unsigned Last_sheet_boundary_id;

// hierher reinstate these
 /// Region ID for elements above disk
 unsigned Above_disk_region_id;

 /// Region ID for elements below disk
 unsigned Below_disk_region_id;

 // One-based region ID for torus
 unsigned One_based_torus_region_id;
  
 /// First boundary ID for torus
 unsigned First_boundary_for_torus_id;
 
 /// Last boundary ID for torus
 unsigned Last_boundary_for_torus_id;

 /// \short One-based boundary IDs of boundaries on disk
 /// and within the torus region
 Vector<unsigned> One_based_boundary_id_for_disk_within_torus;

 /// \short One-based boundary IDs of boundaries on disk
 /// and outside the torus region
 Vector<unsigned> One_based_boundary_id_for_disk_outside_torus;
 
 /// The Line Visualiser.
 LineVisualiser* LV_pt;

 /// \short Number of "disks on disk" around the edge where solution is
 /// to be visualised
 unsigned Ndisk_on_disk_plot;

 /// Mesh as geom object representation of mesh
 MeshAsGeomObject* Mesh_as_geom_object_pt;

 /// \short Number of azimuthal plot points in "disks on disk" plots 
 /// around the edge where solution is to be visualised
 unsigned Nphi_disk_on_disk_plot;

 /// \short Number of radial plot points in "disks on disk" plots 
 /// around the edge where solution is to be visualised
 unsigned Nrho_disk_on_disk_plot;

 /// Disk_on_disk_plot_point[k_disk][i_rho][i_phi]
 Vector<Vector<Vector<std::pair<
        Vector<double>,std::pair<GeomObject*,Vector<double> > > > > >
 Disk_on_disk_plot_point;

};



//========================================================================
/// Constructor for DiskInContainer problem
//========================================================================
template<class ELEMENT>
DiskInContainerProblem<ELEMENT>::DiskInContainerProblem()
{ 

 //Add a steady time stepper
 this->add_time_stepper_pt(new Steady<0>);
 

 // OUTER BOUNDARY
 //===============

 // Start boundary  IDs for outer boundary from some crazy offset
 // (just for testing). By default the one-based boundary IDs go from
 // 1 to 6.
 unsigned outer_boundary_id_offset=1000;

 //Make the outer boundary object
 Outer_boundary_pt = new CubicTetMeshFacetedSurface(
  Global_Parameters::Box_half_width,
  Global_Parameters::Box_half_length,
  outer_boundary_id_offset);

 // Look, we can visualise the faceted surface!
 Outer_boundary_pt->output("outer_facet_thing.dat");

 // oomph-lib boundary id for top
 Top_boundary_id=outer_boundary_id_offset+2;
 
 // oomph-lib boundary id for bottom
 Bottom_boundary_id=outer_boundary_id_offset+5;

 // First boundary ID for outer boundary
 First_boundary_id_for_outer_boundary=outer_boundary_id_offset;
 
 oomph_info << "Top/bottom boundary ids for outer boundary: "
            << Top_boundary_id << " "
            << Bottom_boundary_id << " " 
            << std::endl;
 
 // INTERNAL BOUNDARIES
 //====================
 
 // One-based region ID for torus
 One_based_torus_region_id=1453;
 
 // Disk with torus around edge
 //----------------------------

 // (Half) number of segments used to represent the disk perimeter
 unsigned half_nsegment=30; // hierher 15; 
 
 // Number of vertices around perimeter of torus
 unsigned nvertex_torus=20;
 
 // Start enumeration from here
 unsigned first_one_based_boundary_for_disk_id=19001; 
 
 // These get returned
 unsigned last_one_based_boundary_for_disk_id=0;
 unsigned first_one_based_boundary_for_torus_id=0;
 unsigned last_one_based_boundary_for_torus_id=0;
 
 // Created faceted surface 
 double z_offset=0.0; 
 double r_torus=0.1;
 DiskWithTorusAroundEdgeTetMeshFacetedSurface* srf_pt=
  new DiskWithTorusAroundEdgeTetMeshFacetedSurface(
   z_offset,
   half_nsegment,
   r_torus,
   nvertex_torus,
   first_one_based_boundary_for_disk_id,
   One_based_torus_region_id,
   last_one_based_boundary_for_disk_id,
   first_one_based_boundary_for_torus_id,
   last_one_based_boundary_for_torus_id,
   One_based_boundary_id_for_disk_within_torus,
   One_based_boundary_id_for_disk_outside_torus);
 
 Inner_boundary_pt.push_back(srf_pt);
 
 // Look, we can visualise the faceted surface!
 srf_pt->output("inner_facet_thing_with_torus.dat");
 
 First_sheet_boundary_id=first_one_based_boundary_for_disk_id-1;
 Last_sheet_boundary_id=last_one_based_boundary_for_disk_id-1;
 
 First_boundary_for_torus_id=first_one_based_boundary_for_torus_id-1;
 Last_boundary_for_torus_id=last_one_based_boundary_for_torus_id-1;
 
 

 // Build the mesh
 //--------------- 

 // Initial element volume
 double initial_element_volume=1.0; // was 5.0; reducing to 1.0 for seg fault checking 5.0; // 0.5;

 // How to call gmsh from the command line
 std::string gmsh_command_line_invocation="/home/mheil/gmesh/bin/bin/gmsh";

 // Setup parameters for gmsh
 GmshParameters* gmsh_parameters_pt=
  new GmshParameters(Outer_boundary_pt,
                     gmsh_command_line_invocation);

 gmsh_parameters_pt->element_volume()=initial_element_volume;
 gmsh_parameters_pt->internal_surface_pt()=Inner_boundary_pt;

 gmsh_parameters_pt->stem_for_filename_gmsh_size_transfer()=
  "target_size_on_grid";
 gmsh_parameters_pt->counter_for_filename_gmsh_size_transfer()=0;

 // Problem is linear so we don't need to transfer the solution to the
 // new mesh
 gmsh_parameters_pt->disable_projection();

 // Redirect gmsh on-screen output
 // hierher delete existing file...
 gmsh_parameters_pt->gmsh_onscreen_output_file_name()=
  "RESLT/gmsh_on_screen_output.dat";

 Bulk_mesh_pt = new RefineableGmshTetMesh<ELEMENT>
  (gmsh_parameters_pt,this->time_stepper_pt());


  // Add sub-mesh
  add_sub_mesh(Bulk_mesh_pt);

 // hierher uncomment these again

 // Bulk_mesh_pt->output("mesh.dat");
 // std::ofstream quality_file;
 // quality_file.open("quality.dat");
 // Bulk_mesh_pt->assess_mesh_quality(quality_file);
 // quality_file.close();
 // exit(0);



 // Set error estimator for bulk mesh
 Z2ErrorEstimator* error_estimator_pt=new Z2ErrorEstimator;
 Bulk_mesh_pt->spatial_error_estimator_pt()=error_estimator_pt;

 // Set targets for spatial adaptivity
 Bulk_mesh_pt->max_permitted_error()=0.0005; 
 Bulk_mesh_pt->min_permitted_error()=0.00001;


 // hierher renable
 // Bulk_mesh_pt->max_element_size()=1.0; // hierher use these for uniform element area in adaptation
 // Bulk_mesh_pt->min_element_size()=0.0; 

 /// Mesh as geom object representation of mesh
 Mesh_as_geom_object_pt=0;

 // Number of "disks on disk" around the edge where solution is
 // to be visualised
 Ndisk_on_disk_plot=4;

 // Number of azimuthal plot points in "disks on disk" plots 
 // around the edge where solution is to be visualised
 Nphi_disk_on_disk_plot=30;
 
 // Number of radial plot points in "disks on disk" plots 
 // around the edge where solution is to be visualised
 Nrho_disk_on_disk_plot=50;
 
 // Geom object representations need to be (re)built
 Geom_objects_are_out_of_date=true;
 
 
 if (!CommandLineArgs::command_line_flag_has_been_set
     ("--dont_use_singularity"))
  {
   
   // Create element that stores the singular fct and its amplitude
   //---------------------------------------------------------------
   ScalableSingularityForPoissonElement<ELEMENT>* el_pt=
    new ScalableSingularityForPoissonElement<ELEMENT>;
   
   // Pass fct pointers:
   el_pt->unscaled_singular_fct_pt()
    =&Global_Parameters::singular_fct;
   el_pt->gradient_of_unscaled_singular_fct_pt()=
    &Global_Parameters::gradient_of_singular_fct;
   
   // Add to mesh
   Singular_fct_element_mesh_pt=new Mesh;
   Singular_fct_element_mesh_pt->add_element_pt(el_pt);
   add_sub_mesh(Singular_fct_element_mesh_pt);
   

   // // Create face elements for impositition of BC and flux jump
   // //----------------------------------------------------------
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
 
 oomph_info << "hierher hypre only if no singularity!" << std::endl;

 // Set the linear solver for problem
 linear_solver_pt() = new MumpsSolver;
 
   
 // // Build linear solver
 // linear_solver_pt() = new HSL_MA42;
 
 // /// Switch on full doc
 // static_cast<HSL_MA42*>(linear_solver_pt())->enable_doc_stats();

 // /// Switch on re-ordering
 // static_cast<HSL_MA42*>(linear_solver_pt())->enable_reordering();

 // // Write to disk...
 // static_cast<HSL_MA42*>(linear_solver_pt())->enable_direct_access_files();

 // // Create a new Hypre linear solver
 // HypreSolver* hypre_linear_solver_pt = new HypreSolver;
 
 // // Set the linear solver for problem
 // linear_solver_pt() = hypre_linear_solver_pt;
 
 // // Set some solver parameters
 // hypre_linear_solver_pt->max_iter() = 100;
 // hypre_linear_solver_pt->tolerance() = 1e-10;
 // hypre_linear_solver_pt->amg_simple_smoother() = 1;
 // hypre_linear_solver_pt->disable_doc_time();
 // hypre_linear_solver_pt->enable_hypre_error_messages();
 // hypre_linear_solver_pt->amg_print_level() = 0;
 // hypre_linear_solver_pt->krylov_print_level() = 0;
 // hypre_linear_solver_pt->hypre_method() = HypreSolver::BoomerAMG;
   
 // Setup equation numbering scheme
 oomph_info <<"Number of equations: " << assign_eqn_numbers() << std::endl; 

}


//========================================================================
/// Setup disk on disk plots
//========================================================================
template<class ELEMENT>
void DiskInContainerProblem<ELEMENT>::setup_disk_on_disk_plots()
{
 oomph_info << "Not making geom object" << std::endl;
 return;

 oomph_info << "Starting make geom object" << std::endl;
 double t_start=TimingHelpers::timer();
      
 // Make new geom object
 delete Mesh_as_geom_object_pt;
 Mesh_as_geom_object_pt=new MeshAsGeomObject(Bulk_mesh_pt);

 // Make space for plot points: Disk_on_disk_plot_point[k_disk][i_rho][i_phi]
 Disk_on_disk_plot_point.resize(Ndisk_on_disk_plot);
 for (unsigned i=0;i<Ndisk_on_disk_plot;i++)
  {
   Disk_on_disk_plot_point[i].resize(Nrho_disk_on_disk_plot);
   for (unsigned j=0;j<Nrho_disk_on_disk_plot;j++)
    {
     Disk_on_disk_plot_point[i][j].resize(Nphi_disk_on_disk_plot);
    }
  }

  Vector<double> r_edge(3);
  Vector<double> normal(3);  
  Vector<double> tangent(3); 
  Vector<double> normal_normal(3);   
  Vector<double> x(3);     
  Vector<double> s(3);  
  Vector<double> rho_and_phi(2);
  GeomObject* geom_object_pt=0;
  for (unsigned k=0;k<Ndisk_on_disk_plot;k++)
   {
    double theta=double(k)/double( Ndisk_on_disk_plot)*
     2.0*MathematicalConstants::Pi;
    Global_Parameters::Disk_pt->
     outer_normal_on_boundary(theta,r_edge,normal,tangent);
    normal_normal[0]=normal[1]*tangent[2]-normal[2]*tangent[1];
    normal_normal[1]=normal[2]*tangent[0]-normal[0]*tangent[2];
    normal_normal[2]=normal[0]*tangent[1]-normal[1]*tangent[0];
    for (unsigned i=0;i<Nrho_disk_on_disk_plot;i++)
     {
      double rho_min=0.0;
      double rho_max=0.1;
      double rho=rho_min+(rho_max-rho_min)*double(i)/
       double(Nrho_disk_on_disk_plot-1);
      rho_and_phi[0]=rho;
      for (unsigned j=0;j<Nphi_disk_on_disk_plot;j++)
       {
        double phi=double(j)/double(Nphi_disk_on_disk_plot-1)*
         2.0*MathematicalConstants::Pi;
        rho_and_phi[1]=phi;
        x[0]=r_edge[0]+rho*cos(phi)*normal[0]+rho*sin(phi)*normal_normal[0];
        x[1]=r_edge[1]+rho*cos(phi)*normal[1]+rho*sin(phi)*normal_normal[1];
        x[2]=r_edge[2]+rho*cos(phi)*normal[2]+rho*sin(phi)*normal_normal[2];

        Mesh_as_geom_object_pt->locate_zeta(x,geom_object_pt,s); 
        if (geom_object_pt==0)
         {
          oomph_info << "Point : " 
                     << x[0] << " " 
                     << x[1] << " " 
                     << x[2] << " "
                     << " not found in setup of disk on disk plots" 
                     << std::endl;
         }        
        Disk_on_disk_plot_point[k][i][j]=
         std::make_pair(rho_and_phi,
                        std::make_pair(geom_object_pt,s));
       }
     }
   }

  oomph_info << "Completed setup of disk on disk plots. This took " 
             << TimingHelpers::timer()-t_start << " sec"
             << std::endl;

}

//========================================================================
/// Complete problem setup: No slip everywhere apart from top where
/// we impose parallel, axially traction free outflow
//========================================================================
template<class ELEMENT>
void DiskInContainerProblem<ELEMENT>::complete_problem_setup()
{

 if (!CommandLineArgs::command_line_flag_has_been_set
       ("--dont_use_singularity"))
  {
   
   // Loop over the elements to set up element-specific
   // things that cannot be handled by constructor

   // Bulk elements in torus region
   unsigned region_id=One_based_torus_region_id-1;;
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
   //------------------------------
   ELEMENT* regularised_bulk_el_pt=0;
   Vector<double> s_reg(3);
   
   // Loop over bulk elements in torus region
   double x_max_on_boundary=-DBL_MAX;
   region_id=One_based_torus_region_id-1;
   n_el=Bulk_mesh_pt->nregion_element(region_id);
   for (unsigned e=0;e<n_el;e++)
    {
     ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
      Bulk_mesh_pt->region_element_pt(region_id,e));

     // Now loop over nodes
     double x_max_el_on_boundary=-DBL_MAX;
     double x_max_el=-DBL_MAX;
     bool has_node_on_disk=false;
     unsigned j_right_most_boundary_node=0;
     unsigned nnod=bulk_elem_pt->nnode();
     for (unsigned j=0;j<nnod;j++)
      {
       Node* nod_pt=bulk_elem_pt->node_pt(j);

       // Check if node is on disk, and if so check if it's
       // the rightmost node of all such nodes in this element
       for (unsigned ibound=First_sheet_boundary_id;
            ibound<=Last_sheet_boundary_id;ibound++)
        {
         if (nod_pt->is_on_boundary(ibound))
          {
           has_node_on_disk=true;
           if (nod_pt->x(0)>x_max_el_on_boundary)
            {
             x_max_el_on_boundary=nod_pt->x(0);
             j_right_most_boundary_node=j;
            }
           break;
          }
        }
       
       // Find the rightmost node anywhere in the element
       if (nod_pt->x(0)>x_max_el) x_max_el=nod_pt->x(0);
      }

     // Now: if the element has a node on the disk and its
     // rightmost node (anywhere) is further to the left than
     // its rightmost on the disk it's a candidate for the
     // rightmost element that's attached to the disk:
     if ((has_node_on_disk)&&(x_max_el>x_max_el_on_boundary))
      {
       // Is the rightmost node on the disk further to the right
       // than any previously located one?
       if (x_max_el_on_boundary>x_max_on_boundary)
        {
         x_max_on_boundary=x_max_el_on_boundary;
         regularised_bulk_el_pt=bulk_elem_pt;
         bulk_elem_pt->local_coordinate_of_node(j_right_most_boundary_node,
                                                s_reg);

         oomph_info 
          << "Found new potential regularised bulk element; s_reg ="
          << s_reg[0] << " " << s_reg[1] << " " << s_reg[2] << " coords: "
          << regularised_bulk_el_pt->interpolated_x(s_reg,0) << " "
          << regularised_bulk_el_pt->interpolated_x(s_reg,1) << " "
          << regularised_bulk_el_pt->interpolated_x(s_reg,2) << " "
          << std::endl;
        }
      }
    }

   
   oomph_info 
    << "Found final regularised bulk element; s_reg ="
    << s_reg[0] << " " << s_reg[1] << " " << s_reg[2] << " coords: "
    << regularised_bulk_el_pt->interpolated_x(s_reg,0) << " "
    << regularised_bulk_el_pt->interpolated_x(s_reg,1) << " "
    << regularised_bulk_el_pt->interpolated_x(s_reg,2) << " "
    << std::endl;

   
   // Set pointer to bulk element and local coordinate
   // in it to identify where regularity (zero flux in specified
   // direction) is imposed on the FE solution. Needs to be done
   // after the node pointers have been updated otherwise
   // the classification of external data goes wrong!
   Vector<double> n(3);
   n[0]=1.0;
   n[1]=0.0;
   n[2]=0.0;
   dynamic_cast<ScalableSingularityForPoissonElement<ELEMENT>*>(
    Singular_fct_element_mesh_pt->element_pt(0))->
    specify_regularisation(regularised_bulk_el_pt,s_reg,n);
   

   // oomph_info << "hierher for now pin amplitude of singular fct"
   //            << std::endl;

   // // hierher for now pin: 
   // dynamic_cast<ScalableSingularityForPoissonElement<ELEMENT>*>(
   //  Singular_fct_element_mesh_pt->element_pt(0))->
   //  pin_amplitude_of_singular_fct();

   // double value=1.0;
   // dynamic_cast<ScalableSingularityForPoissonElement<ELEMENT>*>(
   //  Singular_fct_element_mesh_pt->element_pt(0))->
   //  set_amplitude_of_singular_fct(value);


  }

 // Apply bcs
 apply_boundary_conditions();
 
 // Build line visualiser
 build_line_visualiser();
 
 
}

//==start_of_apply_bc=====================================================
/// Helper function to apply boundary conditions
//========================================================================
template<class ELEMENT>
void DiskInContainerProblem<ELEMENT>::apply_boundary_conditions()
{
 
 // Doc pinned nodes
 ofstream pin_file;
 pin_file.open("pinned_nodes.dat");

 // Identify boundary ids of pinned nodes
 Vector<unsigned> pinned_boundary_id;
 for (unsigned ibound=First_sheet_boundary_id;
      ibound<=Last_sheet_boundary_id;ibound++)
  {
   pinned_boundary_id.push_back(ibound);
  }
 for (unsigned ibound=First_boundary_id_for_outer_boundary;
      ibound<=First_boundary_id_for_outer_boundary+6;ibound++)
  {
   pinned_boundary_id.push_back(ibound);
  }

 unsigned num_pin_bnd=pinned_boundary_id.size();
 for (unsigned bnd=0;bnd<num_pin_bnd;bnd++)
  {
   // Pin whatever needs to be pinned
   unsigned ibound=pinned_boundary_id[bnd];
   unsigned num_nod= Bulk_mesh_pt->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     Bulk_mesh_pt->boundary_node_pt(ibound,inod)->pin(0);
     Vector<double> x(3);
     x[0]=Bulk_mesh_pt->boundary_node_pt(ibound,inod)->x(0);
     x[1]=Bulk_mesh_pt->boundary_node_pt(ibound,inod)->x(1);
     x[2]=Bulk_mesh_pt->boundary_node_pt(ibound,inod)->x(2);
     Vector<double> u(1);
     Global_Parameters::get_exact_u(x,u);
     Bulk_mesh_pt->boundary_node_pt(ibound,inod)->set_value(0,u[0]);
     pin_file << x[0] << " " 
              << x[1] << " " 
              << x[2] << " " 
              << std::endl;
    }
  }
 
 pin_file.close();


 if (!CommandLineArgs::command_line_flag_has_been_set
       ("--dont_use_singularity"))
  {
   // Now unpin nodal values where the bc conditions are enforced
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
     double u_bc=1.0; // hierher
     Vector<double> nodal_boundary_value(nnod,u_bc);
     
     // Unpin the fe part of the solution
     for (unsigned j=0;j<nnod;j++)
      {
       el_pt->unpin_u_fe_at_specified_local_node(j);
       
       // Now here's another subtle one! If this node stores three values
       // they come from: original Poisson dof; one Lagrange multiplier
       // to ensure continutuity of solution across "jump"; one Lagrange
       // multiplier from enforcing Dirichlet boundary conditions. 
       // The latter two enforce the same thing, so we can only apply
       // the constraint once --> Pin the Lagrange multiplier that tries to 
       // enforces the Dirichlet BC
       unsigned nval=el_pt->node_pt(j)->nvalue();
       if (nval==3) // hierher
        {
         el_pt->pin_lagrange_multiplier_at_specified_local_node(j);
        }
       
      }
     
     // Tell the element about the desired nodal boundary values
     el_pt->set_nodal_boundary_values(nodal_boundary_value);
     
     
    }
  }

} // end set bc



//========================================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void DiskInContainerProblem<ELEMENT>::doc_solution(const unsigned& nplot,
                                                DocInfo& doc_info)
{ 

 ofstream some_file;
 ofstream some_file2;
 ofstream face_some_file;
 ofstream coarse_some_file;
 char filename[100];


 // Plot normals on disk
 {
  sprintf(filename,"%s/disk_normal%i.dat",
          doc_info.directory().c_str(),
          doc_info.number());
  some_file.open(filename);
  sprintf(filename,"%s/disk_tangent%i.dat",
          doc_info.directory().c_str(),
          doc_info.number());
  some_file2.open(filename);
  unsigned n=100;
  Vector<double> r(3);
  Vector<double> normal(3);  
  Vector<double> tangent(3);
  for (unsigned j=0;j<n;j++)
   {
    double phi=double(j)/double(n-1)*2.0*MathematicalConstants::Pi;
    Global_Parameters::Disk_pt->outer_normal_on_boundary(phi,r,normal,tangent);
    some_file << r[0] << " " 
              << r[1] << " " 
              << r[2] << " " 
              << normal[0] << " " 
              << normal[1] << " " 
              << normal[2] << " " 
              << std::endl;
    some_file2 << r[0] << " " 
               << r[1] << " " 
               << r[2] << " " 
               << tangent[0] << " " 
               << tangent[1] << " " 
               << tangent[2] << " " 
               << std::endl;
   }
  some_file.close();
  some_file2.close();
 }


 // Plot disk
 {
  sprintf(filename,"%s/disk_%i.dat",
          doc_info.directory().c_str(),
          doc_info.number());
  some_file.open(filename);
  Vector<double> r(3);
  Vector<double> x(2);
  unsigned n_phi=100;
  unsigned n_r=200;
  double r_min=0.5;
  double r_max=1.0;
  some_file << "ZONE I=" << n_r << ", J=" << n_phi << std::endl;
  for (unsigned i=0;i<n_r;i++)
   {
    double radius=r_min+(r_max-r_min)*double(i)/double(n_r-1);
    for (unsigned j=0;j<n_phi;j++)
     {
      double phi=double(j)/double(n_phi-1)*2.0*MathematicalConstants::Pi;
      x[0]=radius*cos(phi);
      x[1]=radius*sin(phi);
      Global_Parameters::Disk_pt->position(x,r);
      some_file << r[0] << " " 
                << r[1] << " " 
                << r[2] << " " 
                << std::endl;
     }
   }
  some_file.close();
 }



 // // Plot disks around the perimeter of the disk...
 // {
 //  sprintf(filename,"%s/disk_on_disk%i.dat",
 //          doc_info.directory().c_str(),
 //          doc_info.number());
 //  some_file.open(filename);
 //  Vector<double> r_edge(3);
 //  Vector<double> normal(3);  
 //  Vector<double> tangent(3); 
 //  Vector<double> normal_normal(3);   
 //  Vector<double> x(3);
 //  for (unsigned k=0;k<Ndisk_on_disk_plot;k++)
 //   {
 //    double theta=double(k)/double(Ndisk_on_disk_plot)*
 //     2.0*MathematicalConstants::Pi;
 //    Global_Parameters::Disk_pt->outer_normal_on_boundary(theta,r_edge,normal,tangent);
 //    normal_normal[0]=normal[1]*tangent[2]-normal[2]*tangent[1];
 //    normal_normal[1]=normal[2]*tangent[0]-normal[0]*tangent[2];
 //    normal_normal[2]=normal[0]*tangent[1]-normal[1]*tangent[0];

 //    double rho_min=0.0;
 //    double rho_max=0.1;
 //    some_file << "ZONE I=" << Nrho_disk_on_disk_plot 
 //              << ", J=" << Nphi_disk_on_disk_plot << std::endl;
 //    for (unsigned i=0;i<Nrho_disk_on_disk_plot;i++)
 //     {
 //      double rho=rho_min+(rho_max-rho_min)*double(i)/double(Nrho_disk_on_disk_plot-1);
 //      for (unsigned j=0;j<Nphi_disk_on_disk_plot;j++)
 //       {
 //        double phi=double(j)/double(Nphi_disk_on_disk_plot-1)*
 //         2.0*MathematicalConstants::Pi;
 //        x[0]=r_edge[0]+rho*cos(phi)*normal[0]+rho*sin(phi)*normal_normal[0];
 //        x[1]=r_edge[1]+rho*cos(phi)*normal[1]+rho*sin(phi)*normal_normal[1];
 //        x[2]=r_edge[2]+rho*cos(phi)*normal[2]+rho*sin(phi)*normal_normal[2];
 //        some_file << x[0] << " " 
 //                  << x[1] << " " 
 //                  << x[2] << " " 
 //                  << std::endl;
 //       }
 //     }
 //   }
 //  some_file.close();
 // }


  oomph_info << "Not doing disk on disk plot" << std::endl;

 // // Plot disks around the perimeter of the disk...
 // {
 //  if (Geom_objects_are_out_of_date)
 //   {
 //    // Setup disk on disk plots
 //    setup_disk_on_disk_plots();
    
 //    // Now they're not...
 //    Geom_objects_are_out_of_date=false;
 //   }

 //  sprintf(filename,"%s/disk_on_disk%i.dat",
 //          doc_info.directory().c_str(),
 //          doc_info.number());
 //  some_file.open(filename);
 //  Vector<double> x(3);
 //  for (unsigned k=0;k<Ndisk_on_disk_plot;k++)
 //   {
 //    some_file << "ZONE I=" << Nphi_disk_on_disk_plot 
 //              << ", J=" << Nrho_disk_on_disk_plot << std::endl;
 //    for (unsigned i=0;i<Nrho_disk_on_disk_plot;i++)
 //     {
 //      for (unsigned j=0;j<Nphi_disk_on_disk_plot;j++)
 //       {
 //        (Disk_on_disk_plot_point[k][i][j].second.first)->
 //         position(Disk_on_disk_plot_point[k][i][j].second.second,x);
 //        double rho=(Disk_on_disk_plot_point[k][i][j].first)[0];
 //        double phi=(Disk_on_disk_plot_point[k][i][j].first)[1];
 //        some_file 
 //         << x[0] << " " 
 //         << x[1] << " " 
 //         << x[2] << " " 
 //         << dynamic_cast<ELEMENT*>(Disk_on_disk_plot_point[k][i][j].second.first)->
 //         interpolated_u_poisson(Disk_on_disk_plot_point[k][i][j].second.second) << " "
 //         << rho << " " 
 //         << phi << " " 
 //         << Global_Parameters::asymptotic_solution(rho,phi) << " " 
 //         << std::endl;
 //       }
 //     }
 //   }
 //  some_file.close();
 // }








 //  // hierher
 // unsigned nb=Bulk_mesh_pt->nboundary();
 // oomph_info << "nb " << nb << std::endl;
 // for (unsigned b=0;b<nb;b++)
 //  {
 //   sprintf(filename,"%s/boundary_coordinate%i_%i.dat",
 //           doc_info.directory().c_str(),
 //           b,
 //           doc_info.number());
 //   some_file.open(filename);
 //   Bulk_mesh_pt->Mesh::template doc_boundary_coordinates<ELEMENT>(b,some_file);
 //   some_file.close();
 //  }   


 // for (unsigned b=0;b<nb;b++)
 //  {
 //   sprintf(filename,"%s/boundary_area_coordinate%i_%i.dat",
 //           doc_info.directory().c_str(),
 //           b,
 //           doc_info.number());
 //   some_file.open(filename);
 //   Bulk_mesh_pt->doc_boundary_area_coordinates(b,some_file);
 //   some_file.close();
 //  }   


 // Output elements adjacent to disk boundary 
 sprintf(filename,"%s/elements_next_to_disk%i.dat",
         doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 sprintf(filename,"%s/coarse_elements_next_to_disk%i.dat",
         doc_info.directory().c_str(),
         doc_info.number());
 coarse_some_file.open(filename);
 for (unsigned ibound=First_sheet_boundary_id;
      ibound<=Last_sheet_boundary_id;ibound++)
  {
   unsigned n_el=Bulk_mesh_pt->nboundary_element(ibound);
   for (unsigned e=0;e<n_el;e++)
    {
     Bulk_mesh_pt->boundary_element_pt(ibound,e)->
      output(some_file,nplot);
     Bulk_mesh_pt->boundary_element_pt(ibound,e)->
      output(coarse_some_file,2);
    }
  }
 some_file.close();
 coarse_some_file.close();

 
 // Output elements adjacent to outer boundary
 //-------------------------------------------
 sprintf(filename,"%s/elements_next_to_outer_boundary%i.dat",
         doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 for (unsigned ibound=First_boundary_id_for_outer_boundary;
      ibound<=First_boundary_id_for_outer_boundary+6;ibound++)
  {
   unsigned n_el=Bulk_mesh_pt->nboundary_element(ibound);
   for (unsigned e=0;e<n_el;e++)
    {
     Bulk_mesh_pt->boundary_element_pt(ibound,e)->
      output(some_file,nplot);
    }
  }
 some_file.close();


 // Output boundaries
 //------------------
 sprintf(filename,"%s/boundaries%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 Bulk_mesh_pt->output_boundaries(some_file);
 some_file.close();


 if (!CommandLineArgs::command_line_flag_has_been_set
     ("--dont_use_singularity"))
  {
   
   // Output flux jump elements
   //--------------------------
   sprintf(filename,"%s/flux_jump_elements%i.dat",doc_info.directory().c_str(),
           doc_info.number());
   some_file.open(filename);
   Face_mesh_for_flux_jump_pt->output(some_file,nplot);
   some_file.close();
   
   // Output flux jump element dofs
   //------------------------------
   {
    sprintf(filename,"%s/flux_jump_element_dofss%i.dat",
            doc_info.directory().c_str(),
            doc_info.number());
    some_file.open(filename);
    unsigned nel=Face_mesh_for_flux_jump_pt->nelement();
    for (unsigned e=0;e<nel;e++)
     {
      FiniteElement* el_pt=Face_mesh_for_flux_jump_pt->finite_element_pt(e);
      unsigned nnod=el_pt->nnode();
      for (unsigned j=0;j<nnod;j++)
       {
        unsigned n_dummy=3;
        Vector<int> eqn_number(n_dummy,-9);
        Node* nod_pt=el_pt->node_pt(j);
        some_file << nod_pt->x(0) << " " 
                  << nod_pt->x(1) << " " 
                  << nod_pt->x(2) << " ";
        unsigned nval=nod_pt->nvalue();
        for (unsigned i=0;i<nval;i++)
         {
          eqn_number[i]=nod_pt->eqn_number(i);
         }
        some_file << nval << " ";
        for (unsigned i=0;i<n_dummy;i++)
         {
          some_file << eqn_number[i] << " ";
         }
        some_file << std::endl;
       }
     }
    some_file.close();
   }
  }

 
 // Doc face elements attached to various disk regions
 //---------------------------------------------------
 
 
 //... within torus region
 {
  sprintf(filename,"%s/bulk_elements_on_disk_within_torus%i.dat",
          doc_info.directory().c_str(),
          doc_info.number());
  some_file.open(filename);
  sprintf(filename,"%s/face_elements_on_disk_within_torus%i.dat",
          doc_info.directory().c_str(),
          doc_info.number());
  face_some_file.open(filename);
  unsigned n=One_based_boundary_id_for_disk_within_torus.size();
  for (unsigned i=0;i<n;i++)
   {
    unsigned b=One_based_boundary_id_for_disk_within_torus[i]-1;
    unsigned nel=Bulk_mesh_pt->nboundary_element(b);
    for (unsigned e=0;e<nel;e++)
     {
      FiniteElement* el_pt=Bulk_mesh_pt->boundary_element_pt(b,e);
      el_pt->output(some_file,nplot);     
      
      // What is the index of the face of the bulk element at the boundary
      int face_index = Bulk_mesh_pt->face_index_at_boundary(b,e);
      
      // Build the corresponding prescribed-flux element
      PoissonFluxElement<ELEMENT>* flux_element_pt = new 
       PoissonFluxElement<ELEMENT>(el_pt,face_index);
      flux_element_pt->output(face_some_file,nplot);
      delete flux_element_pt;
     }
   }
  some_file.close();
  face_some_file.close();
 }

 // ... outside torus region
 {
  sprintf(filename,"%s/bulk_elements_on_disk_outside_torus%i.dat",
          doc_info.directory().c_str(),
          doc_info.number());
  some_file.open(filename);
  sprintf(filename,"%s/face_elements_on_disk_outside_torus%i.dat",
          doc_info.directory().c_str(),
          doc_info.number());
  face_some_file.open(filename);
  unsigned n=One_based_boundary_id_for_disk_outside_torus.size();
  for (unsigned i=0;i<n;i++)
   {
    unsigned b=One_based_boundary_id_for_disk_outside_torus[i]-1;
    unsigned nel=Bulk_mesh_pt->nboundary_element(b);
    for (unsigned e=0;e<nel;e++)
     {
      FiniteElement* el_pt=Bulk_mesh_pt->boundary_element_pt(b,e);
      el_pt->output(some_file,nplot);     
      
      // What is the index of the face of the bulk element at the boundary
      int face_index = Bulk_mesh_pt->face_index_at_boundary(b,e);
      
      // Build the corresponding prescribed-flux element
      PoissonFluxElement<ELEMENT>* flux_element_pt = new 
       PoissonFluxElement<ELEMENT>(el_pt,face_index);
      flux_element_pt->output(face_some_file,nplot);
      delete flux_element_pt;
     }
   }
  some_file.close();
  face_some_file.close();
 }



 // Output bulk elements in torus region
 //---------------------------------------
 sprintf(filename,"%s/soln_in_torus%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 unsigned region_id=One_based_torus_region_id-1;
 unsigned n_el=Bulk_mesh_pt->nregion_element(region_id);
 for (unsigned e=0;e<n_el;e++)
  {
   Bulk_mesh_pt->region_element_pt(region_id,e)->output(some_file,nplot);
  }
 some_file.close();

 
 // Output solution
 //----------------
 sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 Bulk_mesh_pt->output(some_file,nplot);
 some_file.close();

 // Output solution showing element outlines
 //-----------------------------------------
 sprintf(filename,"%s/coarse_soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 Bulk_mesh_pt->output(some_file,2);
 some_file.close();

 // Output solution
 //----------------
 sprintf(filename,"%s/soln%i.vtu",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 Bulk_mesh_pt->output_paraview(some_file,nplot);
 some_file.close();

 // Output solution showing element outlines
 //-----------------------------------------
 sprintf(filename,"%s/coarse_soln%i.vtu",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 Bulk_mesh_pt->output_paraview(some_file,2);
 some_file.close();


 // Output line visualiser solution
 //--------------------------------
 if (LV_pt!=0)
  {
   sprintf(filename,"%s/line_soln%i.dat",doc_info.directory().c_str(),
           doc_info.number());
   some_file.open(filename);
   LV_pt->output(some_file);
   some_file.close();
  }

 // Doc pointwise error and compute norm of error and of the solution
 if (Global_Parameters::Epsilon==0.0)
  {
   double error,norm; 
   sprintf(filename,"%s/error%i.dat",doc_info.directory().c_str(),
           doc_info.number());
   some_file.open(filename);
   Bulk_mesh_pt->compute_error(some_file,Global_Parameters::get_exact_u,
                            error,norm);
   some_file.close();
   // Doc error norm:
   cout << "\nNorm of error    : " << sqrt(error) << std::endl;
   cout << "Norm of solution : " << sqrt(norm) << std::endl << std::endl;
   cout << std::endl;
   
   // Exact solution
   sprintf(filename,"%s/exact_soln%i.dat",doc_info.directory().c_str(),
           doc_info.number());
   some_file.open(filename);
   Bulk_mesh_pt->output_fct_paraview(some_file,nplot,
                                  Global_Parameters::get_exact_u);
   some_file.close();
  }

 // What's the amplitude?
 if (!CommandLineArgs::command_line_flag_has_been_set
     ("--dont_use_singularity"))
  {
   oomph_info 
    << "Amplitude of singular function: "
    << dynamic_cast<TemplateFreeScalableSingularityForPoissonElement*>(
     Singular_fct_element_mesh_pt->element_pt(0))->
    amplitude_of_singular_fct() << std::endl;
  }

} // end of doc




//========================================================================
/// Driver
//========================================================================
int main(int argc, char* argv[])
{

 MPI_Helpers::init(argc,argv);

 // Store command line arguments
 CommandLineArgs::setup(argc,argv);
  
 // length of downstream region occupied by impedance elements
 CommandLineArgs::specify_command_line_flag(
  "--dont_use_singularity");

 // Parse command line
 CommandLineArgs::parse_and_assign(); 
 
 // Doc what has actually been specified on the command line
 CommandLineArgs::doc_specified_flags();

 // Note that this can make tetgen die!
 //feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);


 // Shut up prefix
 oomph_info.output_modifier_pt()=&default_output_modifier;

 // Label for output
 DocInfo doc_info;
 
 // Output directory
 doc_info.set_directory("RESLT");
  
 // Number of output points per edge
 unsigned nplot=5;

 // Build problem
 DiskInContainerProblem<ProjectablePoissonElement<
  MyTPoissonElement<3,3> > > problem;


 //Output initial guess
 problem.doc_solution(nplot,doc_info);
 doc_info.number()++;

 // //hierher 
 // DoubleVector residuals;
 // CRDoubleMatrix jacobian;
 // problem.get_jacobian(residuals,
 //                      jacobian);

 // jacobian.sparse_indexed_output("jac.dat");
 // problem.describe_dofs();

 // exit(0);

 unsigned max_adapt=1; // hierher 1;
 for (unsigned i=0;i<=max_adapt;i++)
  {
   // Solve the bastard!
   problem.newton_solve();

   //Output solution
   problem.doc_solution(nplot,doc_info);
 
   //Increment counter for solutions 
   doc_info.number()++;

   if (i!=max_adapt)
    {
     problem.adapt();
    }
  }

}




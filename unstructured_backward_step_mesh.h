


//=====================================================================
/// Helper function to build the mesh; assumed to live in namespace
/// where dimensions of mesh are defined
//=====================================================================
template<class ELEMENT>
RefineableTriangleMesh<ELEMENT>* build_the_mesh(const double& 
                                                uniform_element_area)
{

 // The boundary is bounded by six distinct boundaries, each
 // represented by its own polyline
 Vector<TriangleMeshCurveSection*> boundary_polyline_pt(6);
 TriangleMeshClosedCurve* closed_curve_pt=0;
 
 // Boundary 0: Inflow boundary
 Vector < Vector < double > > bound_0(2);
 for(unsigned ipoint=0; ipoint<2;ipoint++)
  {
   // Resize the vector
   bound_0[ipoint].resize(2);
  }
 bound_0[0][0]=0.0;
 bound_0[0][1]=0.0;
 bound_0[1][0]=0.0;
 bound_0[1][1]=H_up;
 
 unsigned boundary_id=0;
 boundary_polyline_pt[0]=new TriangleMeshPolyLine(bound_0,boundary_id);
 
 // Boundary 1: Top wall (straight)
 Vector < Vector < double > > bound_1(2);
 for(unsigned ipoint=0; ipoint<2;ipoint++)
  {
   // Resize the vector
   bound_1[ipoint].resize(2);
  }
 bound_1[0][0]=0.0;
 bound_1[0][1]=H_up;
 bound_1[1][0]=L_up+L_down;
 bound_1[1][1]=H_up;
 
 boundary_id=1;
 boundary_polyline_pt[1]=new TriangleMeshPolyLine(bound_1,boundary_id);
 
 // Boundary 2: Outflow
 Vector < Vector < double > > bound_2(2);
 for(unsigned ipoint=0; ipoint<2;ipoint++)
  {
   // Resize the vector
   bound_2[ipoint].resize(2);
  }
 bound_2[0][0]=L_up+L_down;
 bound_2[0][1]=H_up;
 bound_2[1][0]=L_up+L_down;
 bound_2[1][1]=H_up-H_down;
 
 boundary_id=2;
 boundary_polyline_pt[2]=new TriangleMeshPolyLine(bound_2,boundary_id);
 
 // Boundary 3: Bottom wall up to start of region surrounded by circle
 Vector < Vector < double > > bound_3(3);
 for(unsigned ipoint=0; ipoint<3;ipoint++)
  {
   // Resize the vector
   bound_3[ipoint].resize(2);
  }
 bound_3[0][0]=L_up+L_down;
 bound_3[0][1]=H_up-H_down;
 bound_3[1][0]=L_up;
 bound_3[1][1]=H_up-H_down;
 bound_3[2][0]=L_up;
 bound_3[2][1]=-Radius_of_internal_boundary;

 boundary_id=3;
 boundary_polyline_pt[3]=new TriangleMeshPolyLine(bound_3,boundary_id);
 

 // Boundary 4: Corner of step, surrounded by circle
 Vector < Vector < double > > bound_4(3);
 for(unsigned ipoint=0; ipoint<3;ipoint++)
  {
   // Resize the vector
   bound_4[ipoint].resize(2);
  }
 bound_4[0][0]=L_up;
 bound_4[0][1]=-Radius_of_internal_boundary;
 bound_4[1][0]=L_up;
 bound_4[1][1]=0.0;
 bound_4[2][0]=L_up-Radius_of_internal_boundary;;
 bound_4[2][1]=0.0;
 
 boundary_id=4;
 TriangleMeshPolyLine* corner_polyline_pt=
  new TriangleMeshPolyLine(bound_4,boundary_id);
 boundary_polyline_pt[4]=corner_polyline_pt;
 

 // Boundary 5: From region surrounded by circle back to inflow
 Vector < Vector < double > > bound_5(2);
 for(unsigned ipoint=0; ipoint<2;ipoint++)
  {
   // Resize the vector
   bound_5[ipoint].resize(2);
  }
 bound_5[0][0]=L_up-Radius_of_internal_boundary;;
 bound_5[0][1]=0.0;
 bound_5[1][0]=0.0;
 bound_5[1][1]=0.0;
 
 boundary_id=5;
 boundary_polyline_pt[5]=new TriangleMeshPolyLine(bound_5,boundary_id);
 


 // Make GeomObject representing Circle for internal boundary
 // around step
 double x_c=L_up;
 double y_c=0.0;
 Circle* circle_pt=new Circle(x_c,y_c,Radius_of_internal_boundary);

  // Number of segments used for representing the curvlinear internal boundary
 unsigned n_segments = 20;
 
 // The intrinsic coordinates for the beginning and end of the curve
 double s_start = -MathematicalConstants::Pi/2.0; // hierher get these from opening angle
 double s_end   =  MathematicalConstants::Pi;

 // State the vertex number for connection on the destination
 // boundaries
  unsigned vertex_to_connect_initial = 0;
  unsigned vertex_to_connect_final = 2;
  boundary_id = 6;
  TriangleMeshCurviLine *boundary6_pt =  
   new TriangleMeshCurviLine(circle_pt,
                             s_start,
                             s_end,
                             n_segments,
                             boundary_id);

  // Do the connection with the destination boundary
  boundary6_pt->connect_initial_vertex_to_polyline(corner_polyline_pt, 
                                                   vertex_to_connect_initial);
  
  // Do the connection with the destination boundary, in this case
  // the connection is done with the outer boundary
  boundary6_pt->connect_final_vertex_to_polyline(corner_polyline_pt,
                                                 vertex_to_connect_final);


  // The open curve that define this boundary is composed of just one
  // curve section
  Vector<TriangleMeshCurveSection*> internal_curve_section_pt(1);
  internal_curve_section_pt[0] = boundary6_pt;
  Vector<TriangleMeshOpenCurve *> inner_open_boundaries_pt(1);
  inner_open_boundaries_pt[0] =
   new TriangleMeshOpenCurve(internal_curve_section_pt);

  // Point identifying "torus" region near backward facing step
  Vector<double> region1_point(2);

  // Define the coordinates of the regions on the domain
  region1_point[0] = L_down+0.5*Radius_of_internal_boundary;
  region1_point[1] = 0.0   +0.5*Radius_of_internal_boundary;

 // Create the triangle mesh polygon for outer boundary
 //----------------------------------------------------
 TriangleMeshPolygon *outer_polygon = 
  new TriangleMeshPolygon(boundary_polyline_pt);
 
 // Enable redistribution of polylines
 outer_polygon -> enable_redistribution_of_segments_between_polylines();
 
 // Set the pointer
 closed_curve_pt = outer_polygon;
 
 // Now build the mesh
 //===================
 
 // Use the TriangleMeshParameters object for helping on the manage of the
 // TriangleMesh parameters
 TriangleMeshParameters triangle_mesh_parameters(closed_curve_pt);
 
 // Specify the maximum area element
 triangle_mesh_parameters.element_area() = uniform_element_area;
 
 // Specify the internal open boundaries
 triangle_mesh_parameters.internal_open_curves_pt() = inner_open_boundaries_pt;

 // Identify region 1
 triangle_mesh_parameters.add_region_coordinates(1, region1_point);

 // Create the mesh
 RefineableTriangleMesh<ELEMENT>* Bulk_mesh_pt=
  new RefineableTriangleMesh<ELEMENT>(triangle_mesh_parameters);
 
 return Bulk_mesh_pt;
}

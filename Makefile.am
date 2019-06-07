# Name of executables: The prefix "noinst" means the executables don't
# have to be installed anywhere.
noinst_PROGRAMS= \
one_d_poisson_with_face_element \
two_d_poisson_with_face_element \
three_d_poisson_with_face_element \
cube_with_disk


#---------------------------------------------------------------------------

# Local sources that Joe's Poisson code depends on:
cube_with_disk_SOURCES = cube_with_disk.cc

# Required libraries: 
cube_with_disk_LDADD = -L@libdir@  -lpoisson -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

# Use this to get the single disk without the surrounding torus
cube_with_disk_CXXFLAGS=-DSINGLE_DISK


#---------------------------------------------------------------------------

# Local sources that Joe's Poisson code depends on:
three_d_poisson_with_face_element_SOURCES = three_d_poisson_with_face_element.cc poisson_sing_face_element.h  tetmesh_faceted_surfaces.h

# Required libraries: 
three_d_poisson_with_face_element_LDADD = -L@libdir@  -lpoisson -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#---------------------------------------------------------------------------

# Local sources that Joe's Poisson code depends on:
two_d_poisson_with_face_element_SOURCES = two_d_poisson_with_face_element.cc poisson_sing_face_element.h unstructured_backward_step_mesh.h

# Required libraries: 
two_d_poisson_with_face_element_LDADD = -L@libdir@  -lpoisson -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#---------------------------------------------------------------------------

# Local sources that Joe's Poisson code depends on:
one_d_poisson_with_face_element_SOURCES = one_d_poisson_with_face_element.cc poisson_sing_face_element.h

# Required libraries: 
one_d_poisson_with_face_element_LDADD = -L@libdir@  -lpoisson -lgeneric $(EXTERNAL_LIBS) $(FLIBS)


#---------------------------------------------------------------------------


# Include path for library headers: All library headers live in 
# the include directory which we specify with -I
# Automake will replace the variable @includedir@ with the actual
# include directory.
INCLUDES = -I@includedir@  

EXTRA_DIST = boundaries_2d.lay boundaries_2d.png convergence_2d.lay run2d.bash




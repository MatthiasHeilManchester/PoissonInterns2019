#! /bin/bash

make two_d_poisson_with_face_element


# Stuff to move
important_files="two_d_poisson_with_face_element two_d_poisson_with_face_element.cc convergence_2d.lay run2d.bash"


# Setup directories
main_dir=NEW_RUNS
if [ -e $main_dir ]; then
    echo "Directory " $main_dir "already exists -- please move it out of the way."
    exit
fi
mkdir $main_dir

cp $important_files $main_dir
cd $main_dir


mkdir RESLT
./two_d_poisson_with_face_element --dont_use_singularity > RESLT/OUTPUT 
mv RESLT RESLT_with_sing_only_no_sing_in_fe

mkdir RESLT
./two_d_poisson_with_face_element > RESLT/OUTPUT 
mv RESLT RESLT_with_sing_only_with_sing_in_fe


mkdir RESLT
./two_d_poisson_with_face_element --add_sin_cosh_to_exact_soln --dont_use_singularity > RESLT/OUTPUT 
mv RESLT RESLT_with_sing_sin_cosh_no_sing_in_fe

mkdir RESLT
./two_d_poisson_with_face_element --add_sin_cosh_to_exact_soln  > RESLT/OUTPUT 
mv RESLT RESLT_with_sing_sin_cosh_with_sing_in_fe


#------------------
#POST PROCESS
#------------------
echo "VARIABLES=\"h\",\"||e||\",\"h<sup>3</sup>\"" > convergence.dat
dir_list=`ls -d RESLT*`
for dir in `echo $dir_list`; do
    echo "ZONE T=\""$dir"\"" >> convergence.dat
    grep "Norm of error" $dir/OUTPUT | awk '{print $15 " " $18 " " 0.1*$15*$15*$15}' >> convergence.dat
done

exit




mkdir RESLT
./two_d_poisson_with_face_element --add_sin_cosh_to_exact_soln --suppress_sing_in_exact_soln > RESLT/OUTPUT 
mv RESLT RESLT_with_sin_cosh_only_with_sing_in_fe



mkdir RESLT
./two_d_poisson_with_face_element --dont_use_singularity > RESLT/OUTPUT 
mv RESLT RESLT_with_sing_only_no_sing_in_fe

mkdir RESLT
./two_d_poisson_with_face_element > RESLT/OUTPUT 
mv RESLT RESLT_with_sing_only_with_sing_in_fe

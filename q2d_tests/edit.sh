# Automated editing of inputdata files. Danger!
basedir='/home/gareth/storage/codes_post_thesis/single_crosssection_postthesis/with_new_variable_suspendedload/q2d_tests'
for i in $(find $basedir -name inputdata2.modin);
    do echo $i;

       #sed 's/read_geo=.true./read_initial_geo=.true./g' $i > $(dirname $i)'/inputdata3.modin'
       
       sed '/read_initial_geo=.true./a\read_initial_waters=.false. ! Do we read in the initial water elevation
            ' $i > $(dirname $i)'/inputdata3.modin'

       #sed 's/read_initial_waters=.false.//' $i > $(dirname $i)'/inputdata3.modin'
     
       # Save old file 
       mv $i $(dirname $i)'/inputdata_old.modin'
       # Copy output file to inputdata2.modin 
       mv $(dirname $i)'/inputdata3.modin' $(dirname $i)'/inputdata2.modin'

       # OVERWRITE THE NEW FILE WITH THE OLD ONE
       #mv $(dirname $i)'/inputdata_old.modin' $i
    done

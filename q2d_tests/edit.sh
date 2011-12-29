# Automated editing of inputdata files. Danger!
basedir='/home/gareth/storage/codes_post_thesis/single_crosssection_postthesis/with_new_variable_suspendedload/q2d_tests'
for i in $(find $basedir -name inputdata2.modin);
    do echo $i;

        ############################################################
        # Various operations that were useful at one time or another
        ############################################################

        ## Renaming a variable 

        #sed 's/read_geo=.true./read_initial_geo=.true./g' $i > $(dirname $i)'/inputdata3.modin'
      
         
        ## Adding a newline after a known line

        #sed '/read_initial_geo=.true./a\read_initial_waters=.false. ! Do we read in the initial water elevation
        #     ' $i > $(dirname $i)'/inputdata3.modin'
 
        #sed '/eddis1D=/a\eddis1D_constant=0.0 ! This adds a constant to the eddy dispersion. Useful for some analytical cases
        #     ' $i > $(dirname $i)'/inputdata3.modin'

        sed '/Cmouth=/a\Cmouth_file="Cmouth" ! File with timeseries of Cmouth, if Cmouth_read=.True.
             ' $i > $(dirname $i)'/inputdata3.modin'

        ## Add a newline before a known line

        #sed '/Cmouth=/i\Cmouth_read=.false. ! Do we read the concentration of suspended sediment at the upstream boundary?
        #     ' $i > $(dirname $i)'/inputdata3.modin'

        ## Removing a variable
        
        #sed 's/read_initial_waters=.false.//' $i > $(dirname $i)'/inputdata3.modin'
      
        ## Save old file 
        mv $i $(dirname $i)'/inputdata_old.modin'
        ## Copy output file to inputdata2.modin 
        mv $(dirname $i)'/inputdata3.modin' $(dirname $i)'/inputdata2.modin'
 
        ## OVERWRITE THE NEW FILE WITH THE OLD ONE -- if you make a mistake
        #mv $(dirname $i)'/inputdata_old.modin' $i
    done

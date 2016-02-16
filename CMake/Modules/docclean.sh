#!/bin/sh

#echo "Want to clean doc dir: $1"


recursive_remove()
{
	# sphinx docu
	if [ -e scai-* ]
	then
		cd scai-*
		#pwd
		rm -rf doctrees html json
		cd ..
	fi
	
	# doxygen docu
	if [ -e system ]
	then
		cd system
		#pwd
		rm -rf html
		cd ..
	fi
	
    for dir in `ls`
    do
        if [ -d $dir ]
        then
	    	cd $dir
	    	recursive_remove
	    	cd ..
        fi
    done
}


if [ -e $1 ]
then
	cd $1
	recursive_remove
	cd ..
fi
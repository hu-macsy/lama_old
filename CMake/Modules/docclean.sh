#!/bin/sh

#echo "Want to clean doc dir: $1"


recursive_remove()
{
	# user docu
	if [ -e user ]
	then
		cd user
		#pwd
		rm -rf doctrees html json
		cd ..
	fi
	
	# system docu
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
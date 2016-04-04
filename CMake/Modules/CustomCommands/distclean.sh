#!/bin/sh

recursive_remove()
{
	#echo $PWD
	rm -f *.so *.cmake install_manifest.txt
	rm -rf tmp src CMakeFiles
	
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

recursive_remove $PWD
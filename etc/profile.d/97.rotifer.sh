ROTIFER_PATH=$(echo $(readlink -f ${BASH_SOURCE[0]}) | sed 's/rotifer\/.*/rotifer/')
append_directory_to_path $ROTIFER_PATH/bin
add_to_environment_variable PYTHONPATH first $ROTIFER_PATH/lib
unset ROTIFER_PATH

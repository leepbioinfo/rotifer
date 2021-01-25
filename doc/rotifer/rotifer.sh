# Example Bash basic configuration for Rotifer
ROTIFER_PATH=$(echo $(readlink -f ${BASH_SOURCE[0]}) | sed 's/rotifer\/.*/rotifer/')
append_directory_to_path $ROTIFER_PATH/bin
add_to_environment_variable PYTHONPATH first $ROTIFER_PATH/lib
add_to_environment_variable PERL5LIB last $ROTIFER_PATH/perl/lib
#export ROTIFER_DATA=/databases
unset ROTIFER_PATH

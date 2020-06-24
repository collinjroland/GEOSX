#!/bin/bash

# check if -a or --all is provided as 1st arg
METHOD=all
case $1 in
    -a|--all) METHOD=all; shift
    ;;
    -g|--git) METHOD=git; shift
    ;;
    *)
esac

# nothing to do if formatting script not specified
if [ -z "$1" ]; then
    echo "Usage: $0 [-a|--all|-g|--git] <PYTHON_EXECUTABLE> <FORMAT_SCRIPT> [<path>...]"
    exit
fi

PYTHON_EXECUTABLE=$1; shift
FORMAT_SCRIPT=$1; shift
LOGFILE=xml_formatting_results.log

# "-r" in GNU xargs omits the call if input is empty
# OS X xargs does not support it, but does the same by default
if [ "$(uname)" == "Darwin" ]; then
    XARGS="xargs"
else
    XARGS="xargs -r"
fi

# check if python/xml is present
if ! $PYTHON_EXECUTABLE -c "import xml" &> /dev/null; then
    >&2 echo "Error: python xml module is required to format xml files"
    exit
fi

# check if git is present, if needed
if [ "$METHOD" = "git" ] && ! (hash git &> /dev/null); then
    >&2 echo "Error: git is required when -g or --git is specified"
    exit
fi

abs_path ()
{
    if [ "$#" -gt 0 ]; then
        realpath -s "$@"
    fi
}

list_xml_files_all () 
{
    abs_path $(find $1 -name "*.xml" -not -path "*/\.*")
}

# git does not have -C flag prior to 1.8.5, so we emulate it with cd
# this is actually the most reliable way, although not the fastest
# run in a subprocess so as to avoid having to cd back
list_xml_files_git ()
{
    local git_root=$(cd $path; git rev-parse --show-toplevel 2>/dev/null)
    if ! [ $? -eq 0 ]; then
        >&2 echo "Error: $path does not appear to be part of a git repository"
        exit 1
    fi
    local prefix=$(cd $path; git rev-parse --show-prefix 2>/dev/null)
    git --git-dir=$git_root/.git ls-files $prefix | grep -e .*[.]xml$ | sed "s|^|$git_root/|g"
}

# create/nullify the log file
echo -n > $LOGFILE

# validate each path separately and write results in the log
for path in "$@"; do
    for file in $(list_xml_files_$METHOD $path); do
        $PYTHON_EXECUTABLE $FORMAT_SCRIPT $file &>> $LOGFILE
    done
done




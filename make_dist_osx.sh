#!/bin/bash
set -e

#Set Script Name variable
SCRIPT=`basename ${BASH_SOURCE[0]}`

#Set fonts for Help.
NORM=`tput sgr0`
BOLD=`tput bold`
REV=`tput smso`

#Help function
function HELP {
    echo -e \\n"${REV}Basic usage:${NORM} ${BOLD}$SCRIPT Maven.app${NORM}"\\n
    echo -e "${REV}-h${NORM}  --Displays this help message. No further functions are performed."\\n
    echo -e "Example: ${BOLD}$SCRIPT \"${apppath}\"${NORM}"\\n
    exit 1
}

#Check the number of arguments. If none are passed, print help and exit.
NUMARGS=$#
if [ $NUMARGS -eq 0 ]; then
    HELP
fi

# A POSIX variable
OPTIND=1         # Reset in case getopts has been used previously in the shell.

# Initialize our own variables:

while getopts "h?f:" opt; do
    case "$opt" in
    h|\?)
        show_help
        exit 0
        ;;
    esac
done

shift $((OPTIND-1))  #This tells getopts to move on to the next argument.

apppath=$1
basepath="${apppath%.app}"
appfn="${apppath##*/}"
distpath="dist"
basefn="${appfn%.app}"
# Get version
VERSION=$("./get_version.sh")
dmgfn="${basefn}_${VERSION}.dmg"

# Set QT Environment
#source qt-5.env

rm -rf "${distpath}"
mkdir -p "${distpath}"

# TODO Should these be in the repo?
#copy support files
# cp bin/*.csv   "${apppath}/Contents/Resources"
# cp bin/*.model "${apppath}/Contents/Resources"

mkdir "${apppath}/Contents/Resources/methods"
cp bin/methods/* "${apppath}/Contents/Resources/methods"
# TODO Should pathways be in the repo?
# mkdir "${apppath}/Contents/Resources/pathways"
# cp bin/pathways/* "${apppath}/Contents/Resources/pathways"
mkdir "${apppath}/Contents/Resources/scripts"
cp bin/scripts/* "${apppath}/Contents/Resources/scripts"

#fix Qt dynamic library dependancy
macdeployqt "${apppath}" -dmg
# TODO This doesn't appear to get built, is that correct or are we missing someething?
# macdeployqt "bin/peakdetector.app" -dmg
mv "${basepath}.dmg" "${distpath}/${dmgfn}"

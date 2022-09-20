# Install/unInstall package files in LAMMPS
# mode = 0/1/2 for uninstall/install/update

mode=$1

# arg1 = file, arg2 = file it depends on

action () {
    if (test $mode = 0) then
        rm -f ../$1
    elif (! cmp -s $1 ../$1) then
        if (test -z "$2" || test -e ../$2) then
            cp $1 ..
            if (test $mode = 2) then
                echo "  updating src/$1"
            fi
        fi
    elif (test -n "$2") then
        if (test ! -e ../$2) then
            rm -f ../$1
        fi
    fi
}

# define LAMMPS_VERSION_NUMBER in pair_nep.cpp
if (test $mode = 1) then
    LAMMPS_VERSION=`grep "VERSION" ../version.h | awk -F "[\"\"]" '{print $2}'`
    LAMMPS_VERSION_NUMBER=`date -d "$LAMMPS_VERSION" +%Y%m%d`
    sed -i "s/NUMBER\s\+[0-9]\{8\}/NUMBER $LAMMPS_VERSION_NUMBER/g" pair_NEP.cpp
fi

# all package files with no dependencies

for file in *.cpp *.h; do
    action $file
done

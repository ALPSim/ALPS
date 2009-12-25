#!/bin/bash
SRC_PATH="vistrails"
BIN_PATH="Contents/Resources/lib/python2.5"
DIRS="api core db gui packages tests"

if [ -z "$1" ] || [ -z "$2" ]
then
    echo "usage: $0 <src_dir> <bin_dir>"
    exit 65
fi

DEST=$2/../dist/Vistrails.app
mkdir -p $DEST
cp -rp $2/* $DEST


for dir in $DIRS
do
    if [ -e "$DEST/$BIN_PATH/$dir" ]
    then
	rm $DEST/$BIN_PATH/$dir
    fi
    cp -rp $1/$SRC_PATH/$dir $DEST/$BIN_PATH/$dir
done

if [ -e "$DEST/$BIN_PATH/../../vistrails.py" ]
then
    rm $DEST/$BIN_PATH/../../vistrails.py
fi
cp -p $1/$SRC_PATH/vistrails.py $DEST/$BIN_PATH/../../vistrails.py
zip $DEST
exit 0

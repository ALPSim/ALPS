ln -s $PWD/subworkflows/ShowAsHTML.xml ~/.vistrails/subworkflows/ShowAsHTML.xml
ln -s $PWD/subworkflows/PersistentResults.xml ~/.vistrails/subworkflows/PersistentResults.xml
ln -s $PWD/userpackages/alps ~/.vistrails/userpackages/alps
ln -s $PWD/userpackages/test ~/.vistrails/userpackages/test
rmdir ~/.vistrails/persistent_files
ln -s $PWD/persistent_files ~/.vistrails/persistent_files

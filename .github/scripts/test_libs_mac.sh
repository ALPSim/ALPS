for f in /opt/alps/bin/* /opt/alps/lib/*dylib /opt/alps/lib/pyalps
do
  echo $f
  otool -L $f | grep /opt/local
  otool -L $f | grep /usr/local
  otool -L $f | grep libzzip
done

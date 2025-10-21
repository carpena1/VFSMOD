#i!/bin/sh
#read old version, new version, old date, new date
#example: sh version.sh v2.4.4 v2.4.5 '04\/2007' '05\/2007'

sed s/$1/$2/g < finput.f > finput.1
sed s/$1/$2/g < hyetograph.f > hyetograph.1
sed s/$1/$2/g < hydrograph.f > hydrograph.1
sed s/$1/$2/g < io.f > io.1
sed s/$3/$4/g < finput.1 > finput.f
sed s/$3/$4/g < hyetograph.1 > hyetograph.f
sed s/$3/$4/g < hydrograph.1 > hydrograph.f
sed s/$3/$4/g < io.1 > io.f
rm *.1

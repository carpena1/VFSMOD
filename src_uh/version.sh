#i!/bin/sh
#read old version, new version, old date, new date
#example: sh version.sh v2.4.4 v2.4.5 '04\/2007' '05\/2007'

sed s/$1/$2/g < outmass.f > outmass.1
sed s/$1/$2/g < finput.f > finput.1
sed s/$1/$2/g < inputs.f > inputs.1
sed s/$1/$2/g < vfsmod.f > vfsmod.1
sed s/$3/$4/g < finput.1 > finput.f
sed s/$3/$4/g < inputs.1 > inputs.f
sed s/$3/$4/g < outmass.1 > outmass.f
sed s/$3/$4/g < vfsmod.1 > vfsmod.f
rm *.1



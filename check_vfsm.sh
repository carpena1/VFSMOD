vfsm 19220105.prj
vfsm sample.prj
vfsm samplePnoRain.prj 
vfsm sampleWTP.prj
vfsm 841215.prj
vfsm sampleP.prj
vfsm sampleP3.prj
vfsm sampleP3nodeg.prj
vfsm samplePnodeg.prj
vfsm samplePnoRain.prj
vfsm sampleWT.prj
vfsm sampleWTP.prj

grep dQ output/*.owq |dm s2 >dQ.txt
grep dE output/*.owq |dm s2 >dE.txt
grep dP output/*.owq |dm s2 >dP.txt
abut dQ.txt dE.txt dP.txt > owq.txt
rm dQ.txt dE.txt dP.txt
cat owq.txt

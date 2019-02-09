all : emat5.exe polydet.exe polyquotient.exe emat6.exe

emat5.exe : Ematrix_5pt.cc Ematrix_6pt.cc polydet.cc polyquotient.cc sturm.cc
	g++ -DMAIN_5PT -DNO_TARGETJR -O3 Ematrix_5pt.cc Ematrix_6pt.cc polydet.cc polyquotient.cc sturm.cc
	mv a.exe emat5.exe

polydet.exe : polydet.cc polyquotient.cc
	g++ -O3 -DPOLYDET_HAS_MAIN polydet.cc polyquotient.cc
	mv a.exe polydet.exe

polyquotient.exe : polyquotient.cc
	g++ -O3 -DPOLYQUOTIENT_HAS_MAIN polyquotient.cc
	mv a.exe polyquotient.exe

emat6.exe :  Ematrix_6pt.cc polydet.cc polyquotient.cc sturm.cc
	g++ -DMAIN_6PT -DNO_TARGETJR -O3 Ematrix_6pt.cc polydet.cc polyquotient.cc sturm.cc
	mv a.exe emat6.exe



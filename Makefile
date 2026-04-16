LIBS:=`root-config --libs`
INCS:=`root-config --cflags`
MYINC    := -I./headers -I./headers/radInts

OBJ := fit_XS

fit_XS: fit_XS.o
	g++ -O3 fit_XS.o ${INCS} ${LIBS} -o fit_XS	-lMinuit  -lRooFitCore -lMathMore -lgsl -lgslcblas -lm
fit_XS.o: $(OBJ).cxx
	g++ -O3 -c  $(OBJ).cxx ${myClasses} ${INCS} ${LIBS} ${MYINC}  -o fit_XS.o    -lMinuit  -lRooFitCore -lMathMore -lgsl -lgslcblas -lm

clean:
	rm -r fit_XS.o

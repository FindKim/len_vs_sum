all: leRCC

leRCC: main.o ExtractSigOrf.o ExtractMMSeq.o RunLength.o RunSum.o
	g++ main.o ExtractSigOrf.o ExtractMMSeq.o RunLength.o RunSum.o -o leRCC

main.o: main.cpp
	g++ -c main.cpp

ExtractSigOrf.o: ExtractSigOrf.cpp ExtractSigOrf.h
	g++ -c ExtractSigOrf.cpp

ExtractMMSeq.o: ExtractMMSeq.cpp ExtractMMSeq.h
	g++ -c ExtractMMSeq.cpp

RCClust.o: RCClust.cpp RCClust.h
	g++ -c RCClust.cpp

clean:
	rm -f *.o *.pyc *~ leRCC
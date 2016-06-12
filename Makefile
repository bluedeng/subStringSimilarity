#Makefile for similarity search and join in complex structure databases
CC=g++ -Wno-deprecated -DNDEBUG -O3
CFLAGS=-c

OFILES = time.o gram.o list.o filtb.o query.o db.o

all:topkPositionSearch topkSubStringMatch_app topkSubStringMatch_exa rangeSubStringMatch

rebuild:clean all
	
clean:
	rm -f *.o

topkPositionSearch:${OFILES} topkPositionSearch.o
	${CC} -o $@ ${OFILES} topkPositionSearch.o -lpthread

topkSubStringMatch_app:${OFILES} topkSubStringMatch_app.o
	${CC} -o $@ ${OFILES} topkSubStringMatch_app.o -lpthread

topkSubStringMatch_exa:${OFILES} topkSubStringMatch_exa.o
	${CC} -o $@ ${OFILES} topkSubStringMatch_exa.o -lpthread

rangeSubStringMatch:${OFILES} rangeSubStringMatch.o
	${CC} -o $@ ${OFILES} rangeSubStringMatch.o -lpthread

# ------------------------------------------
topkPositionSearch.o:topkPositionSearch.cpp
	${CC} ${CFLAGS} $< -o $@ ${INCLUDE}

topkSubStringMatch_app.o:topkSubStringMatch_app.cpp
	${CC} ${CFLAGS} $< -o $@ ${INCLUDE}

topkSubStringMatch_exa.o:topkSubStringMatch_exa.cpp
	${CC} ${CFLAGS} $< -o $@ ${INCLUDE}

rangeSubStringMatch.o:rangeSubStringMatch.cpp
	${CC} ${CFLAGS} $< -o $@ ${INCLUDE}

# COMMON OBJECT FILES

time.o:Time.cpp Time.h
	${CC} ${CFLAGS} $< -o $@ ${INCLUDE}
	
gram.o:Gram.cpp Gram.h
	${CC} ${CFLAGS} $< -o $@ ${INCLUDE}

list.o:GramList.cpp GramList.h
	${CC} ${CFLAGS} $< -o $@ ${INCLUDE}
	
filtb.o:CountFilter.cpp CountFilter.h
	${CC} ${CFLAGS} $< -o $@ ${INCLUDE}

query.o:Query.cpp Query.h
	${CC} ${CFLAGS} $< -o $@ ${INCLUDE}

db.o:SeqDB.cpp SeqDB.h
	${CC} ${CFLAGS} $< -o $@ ${INCLUDE}
CC=cc
SRC = src
BIN=bin
CURRENT_DIR = \"`pwd`\"
CFLAGS = -DSEQUERY_HOME=${CURRENT_DIR}



install: sequery_home resnum_subs sequery
	/bin/mv sequery.exe ${BIN}/sequery
	/bin/rm *.o

all:  sequery_home resnum_subs sequery
	 /bin/mv sequery.exe ${BIN}/sequery

sequery:${SRC}/sequery.c
	${CC} ${CFLAGS} -o sequery.exe ${SRC}/sequery.c sequery_home.o resnum_subs.o

sequery_home:${SRC}/sequery_home.c
	${CC} ${CFLAGS} -c ${SRC}/sequery_home.c

resnum_subs:${SRC}/resnum_subs.c
	${CC} -c ${SRC}/resnum_subs.c

clean:
	/bin/rm *.o



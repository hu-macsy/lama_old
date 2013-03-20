BOOST_ROOT=/home/brandes/local/boost_1_46_0
LAMA_ROOT=/home/brandes/local/lama

CPP_FLAGS = -I${BOOST_ROOT}/include -I${LAMA_ROOT}/include

LD_FLAGS = -L${LAMA_ROOT}/lib -llama -Wl,-rpath=${LAMA_ROOT}/lib 

NAME = simple

EXECUTABLE: ${NAME}

${NAME}: ${NAME}.o
	g++ -o ${NAME} ${NAME}.o ${LD_FLAGS} 

${NAME}.o: ${NAME}.cpp Makefile
	g++ ${CPP_FLAGS} -c ${NAME}.cpp

clean:
	rm -rf ${NAME} *.o

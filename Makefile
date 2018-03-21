INC	= -I/usr/local/include/eigen3
OBJ	= $(SRC:.cpp=.o)
SRC	= amoeba_try.cpp readMatrix.cpp
TARGET	= amoeba
all : $(TARGET) $(OBJS)

$(TARGET) : $(OBJ)
	c++ -o $@ $(OBJ)
.cpp.o :
	c++ -c $< $(INC)
.PHONY : clean
clean :
	rm -f *~ $(TARGET) $(OBJ)


PROJECT (simpleopt)

SET (CMAKE_C_FLAGS_DEBUG "-Wall -Wextra -W -Wshadow -Wcast-qual -Wwrite-strings -g3 -ggdb -O0")

ADD_EXECUTABLE (basicsample basicsample.c simpleopt.c simpleopt.h)


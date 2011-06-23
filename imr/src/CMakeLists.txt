PROJECT (imr)
SET (imr_SRCS 
  "SimpleGlob.h"
  "SimpleOpt.h"
  "tinystr.cpp"
  "tinystr.h"
  "tinyxml"
  "tinyxml.cpp"
  "tinyxml.h"
  "tinyxmlerror.cpp"
  "tinyxmlparser.cpp"
  "imr.cpp"
  "imrTree.cpp"
  "imrTree.h"
  "imrTreeEdge.cpp"
  "imrTreeEdge.h"
)

SET (CMAKE_C_FLAGS_DEBUG "-Wall -Wextra -W -Wshadow -Wcast-qual -Wwrite-strings -g3 -ggdb -O0")

ADD_EXECUTABLE (imr ${imr_SRCS})


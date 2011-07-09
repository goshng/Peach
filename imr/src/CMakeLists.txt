PROJECT (imr)
SET (imr_SRCS 
  "SimpleGlob.h"
  "SimpleOpt.h"
  "tinystr.cpp"
  "tinystr.h"
  "tinyxml.cpp"
  "tinyxml.h"
  "tinyxmlerror.cpp"
  "tinyxmlparser.cpp"
  "imr.cpp"
  "imrTree.cpp"
  "imrTree.h"
  "imrEdge.cpp"
  "imrEdge.h"
  "imrSequence.cpp"
  "imrSequence.h"
  "imrType.h"
  "imrHelp.cpp"
  "imrHelp.h"
  "imrXmlFastsimcoal.cpp"
  "imrXmlFastsimcoal.h"
  "imrXmlImrInput.cpp"
  "imrXmlImrInput.h"
  "imrXmlImrOutput.cpp"
  "imrXmlImrOutput.h"
  "imrDefault.cpp"
  "imrDefault.h"
  "imrXml.cpp"
  "imrXml.h"
)

SET (CMAKE_C_FLAGS_DEBUG "-Wall -Wextra -W -Wshadow -Wcast-qual -Wwrite-strings -g3 -ggdb -O0")

ADD_EXECUTABLE (imr ${imr_SRCS})


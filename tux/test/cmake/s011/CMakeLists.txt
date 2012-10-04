cmake_minimum_required (VERSION 2.8)
set(CMAKE_VERBOSE_MAKEFILE true)

find_package(Perl)

if(PERL_FOUND)
  message("We use find_package to find perl at " ${PERL_EXECUTABLE} " version " ${PERL_VERSION_STRING})
else(PERL_FOUND)
  message("perl is not found with find_package")
endif(PERL_FOUND)

find_path(PERL perl)
if(${PERL} MATCHES NOTFOUND)
  message("We could not find perl: " ${PERL})
else(${PERL} MATCHES NOTFOUND)
  message("We use find_path to find perl at " ${PERL})
endif(${PERL} MATCHES NOTFOUND)


find_path(FUNNYPERL funnyperl)
if(${FUNNYPERL} MATCHES NOTFOUND)
  message("We shoud not find funnyperl:" ${FUNNYPERL})
else(${FUNNYPERL} MATCHES NOTFOUND)
  message("Do we find funnyperl? No!")
endif(${FUNNYPERL} MATCHES NOTFOUND)

find_program(COPYPROGRAM cp)
if(${COPYPROGRAM} MATCHES NOTFOUND)
  message("We could not find cp: " ${COPYPROGRAM})
else(${COPYPROGRAM} MATCHES NOTFOUND)
  message("We use find_program to find cp at " ${COPYPROGRAM})
endif(${COPYPROGRAM} MATCHES NOTFOUND)

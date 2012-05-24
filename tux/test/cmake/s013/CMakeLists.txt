cmake_minimum_required (VERSION 2.8)

function(F1 arg1 arg2)
  message("F1: " ${arg1})
  message("F2: " ${arg2})
endfunction(F1)

F1(3 5)

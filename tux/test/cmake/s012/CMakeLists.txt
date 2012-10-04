cmake_minimum_required (VERSION 2.8)

foreach(MSG a1 a2)
  message(${MSG})
endforeach(MSG a1 a2)

foreach(MSG RANGE 3)
  message(${MSG})
endforeach(MSG RANGE 3)

foreach(MSG RANGE 105 107)
  message(${MSG})
endforeach(MSG RANGE 105 107)


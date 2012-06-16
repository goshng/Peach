#include <iostream>
#include <stdio.h>                                                              
#include <stdlib.h>                                                             
#include <math.h>                                                               
int main (int argc, char *argv[])                                               
{                                                                               
  if (argc < 2)                                                                 
    {                                                                           
    fprintf(stdout,"Usage: %s number\n",argv[0]);                               
    return 1;                                                                   
    }                                                                           
  double inputValue = atof(argv[1]);                                            
  double outputValue = sqrt(inputValue);                                        
  fprintf(stdout,"The square root of %g is %g\n",                               
          inputValue, outputValue);                                             
  double d1 = 3.7;
  unsigned int ui1 = 7;
  double s1 = d1 + ui1;
  std::cout << "s1: " << s1 << " = " << d1 << " + " << ui1 << std::endl;
  return 0;                                                                     
}

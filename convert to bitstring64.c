/*
 *  convert to bitstring.c
 *  bitwise edit distance alignment
 *
 *  Created by Gary Benson on 7/17/12.
 *  Copyright 2012 Boston University. All rights reserved.
 *
 */

#include "convert to bitstring64.h"
#define wordSize64 64

char *convertToBitString64(unsigned long long int number)
{
  int i;
  unsigned long long int mask;
  static char result[wordSize64+1];
  
  mask = 0x0000000000000001;
  for (i=0; i<wordSize64; i++) {
    if (number&mask) {
      result[i]='1';
    }
    else result[i]='0';
    number>>=1;
  }
  return(result);
}



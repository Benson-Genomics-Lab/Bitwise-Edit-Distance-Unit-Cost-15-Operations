/*
 *  bitwise edit distance alignment single word.h
 *  bitwise edit distance alignment
 *
 *  Created by Gary Benson on 7/10/12.
 *  Copyright 2012 Boston University. All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include "bitwise edit distance alignment single word.h"
#include "convert to bitstring64.h"

#define wordSize 64

int Edit_Distance_single_word(char *string1, char *string2, int n, int m){
  //the procedure computes fast bitwise edit distance alignment
  //it uses 15 logic and addition operations
  //this procedure only works if string1 is <=63 characters 
  //n is the length of string1
  //m is the length of string2
  //n and m are not determined by strlen because they 
  //may be part of longer strings
  //returns -1 if error
  
  int i;
  unsigned long long int matchA;
  unsigned long long int matchC;
  unsigned long long int matchG;
  unsigned long long int matchT;
  unsigned long long int matchN;
  unsigned long long int bitmask;
  unsigned long long int matchString;
  unsigned long long int P_D;
  unsigned long long int P_S_or_P_D;
  unsigned long long int M_m;
  unsigned long long int R_IS;
  unsigned long long int not_M_I;
  unsigned long long int not_M_I_xor_P_D;
  unsigned long long int VC_0;
  unsigned long long int VC_plus_1;
  unsigned long long int VC_0_shift;
  unsigned long long int VC_plus_1_shift;
  unsigned long long int sum;
  unsigned long long int sum_and_R_IS;
  
  unsigned long long int P_DErase;
  unsigned long long int P_S_or_P_DErase;
  int PDScore;
  int PSorPDScore;
  int editDistanceScore;
  int countOneBits;
  
  //length of strings for LCS is n and m
  //this procedure only works if string1 is <=63 characters 
  if (n>=wordSize) {
    printf("\nError, first string is longer than %d characters. n=%d",wordSize-1,n);
    return(-1);
  }
  
  //*************************encode match strings A C G T N for string1 
  //loop through string1 and store bits in matchA, matchC, etc.
  //position zero corresponds to column zero in the score matrix, i.e., no character
  //so we start with i = 0 and bitmask = 2
  bitmask=0x0000000000000002;
  matchA=0x0000000000000000;
  matchC=0x0000000000000000;
  matchG=0x0000000000000000;
  matchT=0x0000000000000000;
  matchN=0x0000000000000000;	
  for(i=0;i<n;i++)
  {
    //printf("\nbitmask %6llu",bitmask);
    switch (string1[i]) {
      case 'A':
        matchA |= bitmask;
        break;
      case 'C':
        matchC |= bitmask;
        break;
      case 'G':
        matchG |= bitmask;
        break;
      case 'T':
        matchT |= bitmask;
        break;
      case 'N':
        matchN |= bitmask;
        break;
      default:
        printf("\nError, non-ACGTN character read in string:%c",string1[i]);
        return(-1);
        break;
    }
    bitmask <<= 1; //bitmask = bitmask<<1
    
  }
  
  //debug
  /*
   printf("\nmatchA: %s",convertToBitString64(matchA));
   printf("\nmatchC: %s",convertToBitString64(matchC));
   printf("\nmatchG: %s",convertToBitString64(matchG));
   printf("\nmatchT: %s",convertToBitString64(matchT));
   printf("\nmatchN: %s",convertToBitString64(matchN));
  */
  
  //initialize PD and PS||PD for row zero
  //all of row zero is increase, except zero position
  //which must be S or D so that an initial first R_IS run is computed correctly
  //here it is set to D for global alignment so that the vertical class value
  //below it will be set to +1 (first column increasing)
  //otherwise, set to S for semi-global alignment so no penalty 
  //for starting row (first column all zeros)
  
  P_D = 1;        //puts 1 in zero position, as required above for global alignment
  P_S_or_P_D = 1; //required for consistency
   
  //debug
  printf("\n\nRow 0");
  printf("\nP_D             %s",convertToBitString64(P_D));
  printf("\nP_S_or_P_D      %s",convertToBitString64(P_S_or_P_D));
 
  //loop for each letter in string2
  for(i=0;i<m;i++)
  {
    switch (string2[i]) {
      case 'A':
        matchString = matchA;
        break;
      case 'C':
        matchString = matchC;
        break;
      case 'G':
        matchString = matchG;
        break;
      case 'T':
        matchString = matchT;
        break;
      case 'N':
        matchString = matchN;
        break;
      default:
        printf("\nError, non-ACGTN character read in string:%c",string2[i]);
        return(-1);
        break;
    }
    
    //compute bitstrings
    //match subsets
    M_m = matchString|P_D;
    R_IS = ~M_m;
    not_M_I = R_IS|P_S_or_P_D;
    not_M_I_xor_P_D = not_M_I^P_D;
    //the add
    sum = not_M_I+P_S_or_P_D;
    sum_and_R_IS = sum&R_IS;
    //the new vertical classes
    VC_0 = sum_and_R_IS^not_M_I_xor_P_D;
    VC_plus_1 = P_D|(sum_and_R_IS&P_S_or_P_D);
    //the vertical class shifts
    VC_0_shift=VC_0<<1;
    VC_plus_1_shift=VC_plus_1<<1;
    //reset first position in VC_plus_1_shift for boundary condition
    VC_plus_1_shift+=0x0000000000000001;  
    
    //new pair classes
    P_D = M_m&VC_plus_1_shift;
    P_S_or_P_D = VC_plus_1_shift|(M_m&VC_0_shift);

    //debug
    /*
    printf("\n\nRow %d",i+1);
    printf("\nmatchString     %s",convertToBitString64(matchString));
    printf("\nM_m             %s",convertToBitString64(M_m));
    printf("\nR_IS            %s",convertToBitString64(R_IS));
    printf("\nnot_M_I         %s",convertToBitString64(not_M_I));
    printf("\nnot_M_I_xor_P_D %s",convertToBitString64(not_M_I_xor_P_D));
    printf("\nsum             %s",convertToBitString64(sum));
    printf("\nsum_and_R_IS    %s",convertToBitString64(sum_and_R_IS));
    printf("\nVC_0            %s",convertToBitString64(VC_0));
    printf("\nVC_0_shift      %s",convertToBitString64(VC_0_shift));
    printf("\nVC_plus_1       %s",convertToBitString64(VC_plus_1));
    printf("\nVC_plus_1_shift %s",convertToBitString64(VC_plus_1_shift));
    printf("\nP_D             %s",convertToBitString64(P_D));
    printf("\nP_S_or_P_D      %s",convertToBitString64(P_S_or_P_D));
    */
    
    
  }
  
  //debug
  printf("\n\nFinal:");
  printf("\nP_D:            %s",convertToBitString64(P_D));
  printf("\nP_S_or_P_D:     %s",convertToBitString64(P_S_or_P_D));
  
  //find edit distance score
  //time proportional to number of ones in bit strings
  P_DErase = P_D;
  countOneBits = 0 ;
  while (P_DErase)//stop when zero
  {
    countOneBits++;
    P_DErase &= (P_DErase - 1); //removes last one bit
  }
  PDScore=countOneBits; 
  
  
  P_S_or_P_DErase = P_S_or_P_D;
  countOneBits = 0 ;
  while (P_S_or_P_DErase)//stop when zero
  {
    countOneBits++;
    P_S_or_P_DErase &= (P_S_or_P_DErase - 1); //removes last one bit
  }
  PSorPDScore=countOneBits;//
  
  //debug
  printf("\nPDScore: %d",PDScore);
  printf("\nPSorPDScore: %d",PSorPDScore);
  
  //add +2 to cancel effect of initial 1s from PD and PS||PD due to global alignment
  //the initial 1s shouldn't be counted
  editDistanceScore = m + n - PSorPDScore - PDScore + 2;
  
  
  //printf("\nSingle Word Unit Cost Edit Distance = %d",editDistanceScore);
  
  return(editDistanceScore);
  
  
}

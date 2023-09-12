/*
 *  bitwise edit distance alignment multiple word clean.c
 *  bitwise edit distance alignment
 *
 *  Created by Gary Benson on 7/15/12.
 *  Copyright 2012 Boston University. All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include "bitwise edit distance alignment multiple word clean.h"
#include "bitwise edit distance alignment single word.h"
#include "convert to bitstring64.h"

#define wordSize 64
#define wordSizeMinusOne 63

#define freeMultipleWords \
free(matchA);\
free(matchC);\
free(matchG);\
free(matchT);\
free(matchN);\
free(P_D);\
free(P_S_or_P_D);


int Edit_Distance_multiple_word_clean(char *stringN, char *stringM, int N, int M){
  //the procedure computes fast bitwise edit distance
  //n is the length of string1
  //m is the length of string2
  //n and m are *not* determined by strlen because they
  //may be part of longer strings
  //this procedure does *not* switch strings to make the longer string as string1
  //best if longer string is string1
  //returns -1 if error
  
  int i,j;
  unsigned long long int *matchA;
  unsigned long long int *matchC;
  unsigned long long int *matchG;
  unsigned long long int *matchT;
  unsigned long long int *matchN;
  unsigned long long int *matchVector; //to hold one of matchA, matchC, etc.  no memory allocated
  unsigned long long int bitmask;
  unsigned long long int matchString;
  unsigned long long int *P_D;
  unsigned long long int *P_S_or_P_D;
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
  unsigned long long int highBitMask64 = 0x8000000000000000;
  unsigned long long int carryBitSum;
  unsigned long long int carryBitVC_plus_1_shift;
  unsigned long long int carryBitVC_0_shift;
  unsigned long long int oldCarryBitVC_plus_1_shift;
  unsigned long long int oldCarryBitVC_0_shift;
  
  
  unsigned long long int P_DErase;
  unsigned long long int P_S_or_P_DErase;
  int PDScore;
  int PSorPDScore;
  int editDistanceScore;
  int countOneBits;
  int NWords;
  int afterWords;
  int junkBits;
  unsigned long long int junkBitsMask;
  
  
  //compute number of wordSize-bit words required (NWords).
  //all computation is done in wordSize-bit words.
  //first bit word can only hold wordSize-1 string positions
  //so there is space for the zero column
  
  afterWords = N/wordSize;
  NWords = 1 + afterWords;
  
  //junkbits is used to mask incorrect bits past the end of the string
  //produced by the calculations
  //this is needed to correctly compute the score in time proportional to the
  //number of 1 bits in the results (see score computation below)
  junkBits = wordSizeMinusOne + afterWords*wordSize - N;
  junkBitsMask = 0xFFFFFFFFFFFFFFFF >> junkBits;
  
  //debug
  printf("\n\njunkBits: %d",junkBits);
  printf("\n\njunkBitsMask            %s",convertToBitString64(junkBitsMask));
  printf("\n");
  for (i=0; i<strlen("junkBitsMask            ")+wordSize-junkBits; i++) {
    printf(" ");
  }
  printf("^junkbits start here\n");
  
  
  //storage allocation
  matchA = (unsigned long long int *)calloc(NWords,sizeof(unsigned long long int));
  matchC = (unsigned long long int *)calloc(NWords,sizeof(unsigned long long int));
  matchG = (unsigned long long int *)calloc(NWords,sizeof(unsigned long long int));
  matchT = (unsigned long long int *)calloc(NWords,sizeof(unsigned long long int));
  matchN = (unsigned long long int *)calloc(NWords,sizeof(unsigned long long int));
  P_D = (unsigned long long int *)calloc(NWords,sizeof(unsigned long long int));
  P_S_or_P_D = (unsigned long long int *)calloc(NWords,sizeof(unsigned long long int));
  
  //encode match strings A C G T N for string1
  //loop through stringN and store bits in matchA, matchC, etc.
  //column zero in the score matrix is ignored
  //so we start with i = 0 and bitmask = 1 in the first nWord only
  
  for(j=0;j<NWords;j++)
  {
    bitmask=0x0000000000000001;
    
    matchA[j]=0x0000000000000000;
    matchC[j]=0x0000000000000000;
    matchG[j]=0x0000000000000000;
    matchT[j]=0x0000000000000000;
    matchN[j]=0x0000000000000000;
    for(i=j*wordSize;i<(j+1)*wordSize;i++)
    {
      if (i<=N)
      {
        if(i||j)
        {
          //if both i and j are zero, it means we are at the zero position in the row
          //don't process because it doesn't correspond to a character in the string
          
          switch (stringN[i-1]) {
            case 'A':
              matchA[j] |= bitmask;
              break;
            case 'C':
              matchC[j] |= bitmask;
              break;
            case 'G':
              matchG[j] |= bitmask;
              break;
            case 'T':
              matchT[j] |= bitmask;
              break;
            case 'N':
              matchN[j] |= bitmask;
              break;
            default:
              printf("\nError, non-ACGTN character read at position %d in string:%c",i,stringN[i]);
              freeMultipleWords;
              return(-1);
              break;
          }
        }
        bitmask <<=1; //bitmask = bitmask<<1, moves set bit one position higher
      }
    }
  }
  
  //debug
  /*
   printf("\nmatchA:  ");for(j=0;j<NWords;j++){printf("%s ",convertToBitString64(matchA[j]));}
   printf("\nmatchC:  ");for(j=0;j<NWords;j++){printf("%s ",convertToBitString64(matchC[j]));}
   printf("\nmatchG:  ");for(j=0;j<NWords;j++){printf("%s ",convertToBitString64(matchG[j]));}
   printf("\nmatchT:  ");for(j=0;j<NWords;j++){printf("%s ",convertToBitString64(matchT[j]));}
   printf("\nmatchN:  ");for(j=0;j<NWords;j++){printf("%s ",convertToBitString64(matchN[j]));}
   */
  
  //initialize PD and PS||PD for row zero
  //all of row zero is increase, except zero position
  //which must be S or D so that an initial first R_IS run is computed correctly
  //here it is set to D for global alignment so that the vertical class value
  //below it will be set to +1 (first column increasing)
  //otherwise, set to S for semi-global alignment so no penalty
  //for starting row (first column all zeros)
  
  P_D[0] = 1;        //puts 1 in zero position, as required above for global alignment
  P_S_or_P_D[0] = 1; //required for consistency
  
  //debug
  //printf("\n\nRow 0");
  //printf("\nP_D         ");for(j=0;j<NWords;j++){printf("%s ",convertToBitString64(P_D[j]));}
  //printf("\nP_S_or_P_D  ");for(j=0;j<NWords;j++){printf("%s ",convertToBitString64(P_S_or_P_D[j]));}
  
  //loop for each letter in string2
  for(i=0;i<M;i++)
  {
    carryBitSum = 0x0000000000000000;
    carryBitVC_plus_1_shift = 0x0000000000000001; //set up for adding 1 to first position after shift
    carryBitVC_0_shift = 0x0000000000000000;
    
    //process each row
    //select match vector array for row
    switch (stringM[i]) {
      case 'A':
        matchVector = matchA;
        break;
      case 'C':
        matchVector = matchC;
        break;
      case 'G':
        matchVector = matchG;
        break;
      case 'T':
        matchVector = matchT;
        break;
      case 'N':
        matchVector = matchN;
        break;
      default:
        printf("\nError, non-ACGTN character read in string:%c",stringM[i]);
        freeMultipleWords;
        return(-1);
        break;
    }
    
    //process each nWord
    for(j=0;j<NWords;j++)
    {
      
      //debug
      //printf("\n\nRow %d, Word %d",i,j);
      //printf("\nP_D                     %s",convertToBitString64(P_D[j]));
      //printf("\nP_S_or_P_D              %s",convertToBitString64(P_S_or_P_D[j]));
      
      
      oldCarryBitVC_plus_1_shift = carryBitVC_plus_1_shift;
      oldCarryBitVC_0_shift = carryBitVC_0_shift;
      
      matchString = matchVector[j];
      
      //compute bitstrings
      //match subsets
      M_m = matchString|P_D[j];
      R_IS = ~M_m;
      not_M_I = R_IS|P_S_or_P_D[j];
      not_M_I_xor_P_D = not_M_I^P_D[j];
      
      //the add (with carry)
      //and get carryBitSum for next round
      //saw this online. if a+b causes overflow (carry), a+b<a, a+b<b.
      //bug fix Beth Becker 4/30/14
      //old: sum = not_M_I+P_S_or_P_D[j]+carryBitSum;
      //old: if(sum<not_M_I) carryBitSum = 1;
      sum = not_M_I+carryBitSum;
      carryBitSum = (sum<not_M_I)?1:0;
      sum += P_S_or_P_D[j];
      carryBitSum |= (sum<P_S_or_P_D[j])?1:0;
      
      sum_and_R_IS = sum&R_IS;
      
      //the new vertical classes
      VC_0 = sum_and_R_IS^not_M_I_xor_P_D;
      VC_plus_1 = P_D[j]|(sum_and_R_IS&P_S_or_P_D[j]);
      
      //get shift carry bits
      //mask and shift
      carryBitVC_0_shift = (VC_0&highBitMask64)>>wordSizeMinusOne;
      carryBitVC_plus_1_shift = (VC_plus_1&highBitMask64)>>wordSizeMinusOne;
      
      //the vertical class shifts (with carry)
      VC_0_shift=(VC_0<<1)+oldCarryBitVC_0_shift;
      VC_plus_1_shift=(VC_plus_1<<1)+oldCarryBitVC_plus_1_shift;
      
      //new pair classes
      P_D[j] = M_m&VC_plus_1_shift;
      P_S_or_P_D[j] = VC_plus_1_shift|(M_m&VC_0_shift);

      //debug
      /*
      printf("\n\nRow %d, Word %d",i+1,j);
      printf("\nmatchString             %s",convertToBitString64(matchString));
      printf("\nM_m                     %s",convertToBitString64(M_m));
      printf("\nR_IS                    %s",convertToBitString64(R_IS));
      printf("\nnot_M_I                 %s",convertToBitString64(not_M_I));
      printf("\nnot_M_I_xor_P_D         %s",convertToBitString64(not_M_I_xor_P_D));
      printf("\nsum                     %s",convertToBitString64(sum));
      printf("\nsum_and_R_IS            %s",convertToBitString64(sum_and_R_IS));
      printf("\nVC_0                    %s",convertToBitString64(VC_0));
      printf("\nVC_0_shift              %s",convertToBitString64(VC_0_shift));
      printf("\nVC_plus_1               %s",convertToBitString64(VC_plus_1));
      printf("\nVC_plus_1_shift         %s",convertToBitString64(VC_plus_1_shift));
      
      printf("\n\ncarryBitSum             %s",convertToBitString64(carryBitSum));
      printf("\ncarryBitVC_0_shift      %s",convertToBitString64(carryBitVC_0_shift));
      printf("\ncarryBitVC_plus_1_shift %s",convertToBitString64(carryBitVC_plus_1_shift));
      printf("\nP_D                     %s",convertToBitString64(P_D[j]));
      printf("\nP_S_or_P_D              %s",convertToBitString64(P_S_or_P_D[j]));
      */
    }
    
    //debug
    //printf("\n\nrow %d:",i+1);
    //printf("\nP_D:         ");for(j=0;j<NWords;j++){printf("%s ",convertToBitString64(P_D[j]));}
    //printf("\nP_S_or_P_D:  ");for(j=0;j<NWords;j++){printf("%s ",convertToBitString64(P_S_or_P_D[j]));}
  }
  
  //debug
  
  //debug
  printf("\n\nFinal:");
  printf("\nP_D:         ");
  for(j=0;j<NWords;j++){
    if (j==NWords-1) printf("\n%s ",convertToBitString64(P_D[j]&junkBitsMask));
    else printf("\n%s ",convertToBitString64(P_D[j]));
  }
  printf("\n");
  for (i=0; i<wordSize-junkBits; i++) {
    printf(" ");
  }
  printf("^junkbits start here\n");
  
  //printf("\nP_I = ~P_S_or_P_D:  ");
  //for(j=0;j<NWords;j++){
  //  if (j==NWords-1) printf("\n%s ",convertToBitString64(~P_S_or_P_D[j]&junkBitsMask));
  //  else printf("\n%s ",convertToBitString64(~P_S_or_P_D[j]));
  //}
  printf("\nP_S_or_P_D:");
  for(j=0;j<NWords;j++){
    if (j==NWords-1) printf("\n%s ",convertToBitString64((P_S_or_P_D[j]&junkBitsMask)));
    else printf("\n%s ",convertToBitString64(P_S_or_P_D[j]));
  }
  printf("\n");
  for (i=0; i<wordSize-junkBits; i++) {
    printf(" ");
  }
  printf("^junkbits start here\n");
  
  
  //find edit distance score
  //time proportional to number of ones in bit strings
  //first mask junk bits past end of sequence
  
  PDScore=0;
  PSorPDScore=0;
  for(j=0;j<NWords;j++){
    P_DErase = P_D[j];
    if (j==NWords-1) {
      P_DErase &= junkBitsMask;
    }
    countOneBits = 0 ;
    while (P_DErase)//stop when zero
    {
      countOneBits++;
      P_DErase &= (P_DErase - 1); //removes last one bit
    }
    PDScore+=countOneBits;
    
    
    P_S_or_P_DErase = P_S_or_P_D[j];
    if (j==NWords-1) {
      P_S_or_P_DErase &= junkBitsMask;
    }
    countOneBits = 0 ;
    while (P_S_or_P_DErase)//stop when zero
    {
      countOneBits++;
      P_S_or_P_DErase &= (P_S_or_P_DErase - 1); //removes last one bit
    }
    PSorPDScore+=countOneBits;//
  }
  
  //debug
  printf("\nPDScore: %d",PDScore);
  printf("\nPSorPDScore: %d",PSorPDScore);
  
  //add +2 to cancel effect of initial 1 from PD and PS||PD due to global alignment
  //the initial 1 shouldn't be counted
  editDistanceScore = M + N - PSorPDScore - PDScore + 2;
  
  
  freeMultipleWords;
  return(editDistanceScore);
  
  
  
}


/*   File:                /home/edeveaud/Work/apoptose/src/params.h
 *   Author:              Eric Deveaud edeveaud@pasteur.fr
 *   Date first visited:  "Tue Apr 23 2002"
 *   Time-stamp:          "Thu May 16 2002"
 *   Development Stage :  Under construction
 */


#ifndef __PARAMS_H_
#define __PARAMS_H_

/* how many occurences max we expect */
#define NB_OCCURENCES 10
#define TRUE 1
#define FALSE 0
#define BUFFLEN 100


typedef struct param_S {
  char *pattern;
  char *outfile;
  FILE *OUT;            /* where to display the results */
  int overlap;		/* shall the motif be overlaping */
  int reverse;
}param_t;



#endif /* __PARAMS_H_ */

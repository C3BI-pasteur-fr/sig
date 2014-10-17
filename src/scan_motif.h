/*   File:                /home/edeveaud/Work/apoptose/src/scan_motif.h
 *   Author:              Eric Deveaud edeveaud@pasteur.fr
 *   Date first visited:  "Thu Apr 25 2002"
 *   Time-stamp:          "Tue May 14 2002"
 *   Development Stage :  Under construction / distrib / final
 */


#ifndef __SCAN_MOTIF_H_
#define __SCAN_MOTIF_H_

#include "regex.h"
#include "seq-reader.h"

typedef struct mot_S {
  char *motif;	/* the motif in human expression form */
  regex_t *reg; /* the same in regular expression form */
  int pos;
  int min;
  int max;
}mot_t;

typedef struct motif_S {
  char *motif;
  int start;
  int stop;
  int len;
}motif_t;



//int search_motifs (FILE *OUT, seq_t *seq, mot_t *mot, int n) ;
int search_motifs (param_t params, seq_t *seq, mot_t *mot, int n) ;

#endif  /* __SCAN_MOTIF_H_ */

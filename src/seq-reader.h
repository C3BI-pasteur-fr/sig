/*   File:                /home/edeveaud/Work/toppred/src/seq-reader.h
 *   Author:              Eric Deveaud edeveaud@pasteur.fr
 *   Date first visited:  "Tue Jun 26 2001"
 *   Time-stamp:          "Tue Jul 24 2001"
 *   Development Stage :  Under construction
 */


#ifndef __SEQ_READER_H_
#define __SEQ_READER_H_

#include "params.h"



typedef struct seq_S {
  char *id;
  char *comment;
  char *seq;
  int size;
} seq_t;

#define IDLEN 8
#define COMLEN 100

/* retreive sequence in fasta format from file */
int read_seq(FILE *IN, seq_t  *seq_holder) ;

/* free the sequence holder */
void free_seq(seq_t *seq);

#endif /* __SEQ_READER_H_ */

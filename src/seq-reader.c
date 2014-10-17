/*   File:                /home/edeveaud/Work/toppred/src/seq-reader.c
 *   Author:              Eric Deveaud edeveaud@pasteur.fr
 *   Date first visited:  "Tue Jun 26 2001"
 *   Time-stamp:          "Tue Jun 04 2002"
 *   Development Stage :  Under construction
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>

#ifdef STDC_HEADERS
#include <stdlib.h>
#include <string.h>
#endif

#include <ctype.h>

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#include "seq-reader.h"
#include "error.h"

#define COMMENT 0
#define SEQ 1

int read_seq(FILE *IN, seq_t  *seq_holder) {

  char *BUFF, *id, *comment, *seq, *p, *q;
  int i, state, n, n_tot, len;
  n_tot = BUFFLEN;
  n = 0;
  len = BUFFLEN;

  if((BUFF = malloc(sizeof(char) * (BUFFLEN + 1))) == NULL){
    error_fatal("memory", NULL);
  }

  p = NULL;

  /* procces id and comment */
  state = COMMENT;
  while ((fgets(BUFF, len, IN)) != NULL) {

    /* comment line is not complety read */
    if (strrchr(BUFF, '\n') == NULL) {
      len += BUFFLEN;
      if((BUFF = realloc(BUFF, sizeof(char)*len)) == NULL) {
	error_fatal("memory realloc", NULL);
      }
      if(fseek(IN, -1 *(long)strlen(BUFF), SEEK_CUR)){
	error_fatal("fseek", "fseek");
      }
      continue;
    }

    /* procces id and comment */
    else{

      p = BUFF;
      if (*p != '>') {
	error_fatal("file", "not in fasta format");
      }
      p++;
      while(*p && !isspace((int)*p)) {
	p++;
	n++;
      }

      /* check if sequence is named or not */
      if((n == 0) || (n == 1 && ispunct((int)BUFF[1]))) {
	if((id = malloc(sizeof(char)*(9+1))) == NULL){
 	  error_fatal("memory", NULL);
 	}
 	error_warn("anonymous sequence",
 		   "name will be forced to \"anonymous\"");
 	snprintf(id,10, "anonymous");
      }

      else {
	if((id = malloc(sizeof(char)*(n+1))) == NULL){
	  error_fatal("memory", NULL);
	}
	q = id;
	p = BUFF;
	p++;

	for(i = 0; i<n; i++){
	  *q++ = *p++;
	}
	*q = '\0';
      }

      while(*p && (isspace((int)*p)  || ispunct((int)*p)))
	p++;

      if((comment = malloc((strlen(BUFF) - n)*sizeof(char))) == NULL){
 	error_fatal("memory", "COMLEN");
      }

      q = comment;
      while(*p && *p != '\n') *q++ = *p++;
      *q = '\0';
      state = SEQ;
      break;
    }
  }

  n = 0;
  /*process sequence */
  if(state == SEQ){

    if((seq = malloc(sizeof(char)* (BUFFLEN+1))) == NULL) {
      error_fatal("memory", NULL);
    }
    q =seq;
    while((fgets(BUFF, BUFFLEN, IN)) != NULL) {

      p = BUFF;
      /* allow multiple sequence files to be proceced */
      if (state == SEQ && *p && *p == '>') {
	/* replace the buffer on input*/
	if(fseek(IN, -1 * strlen(p), SEEK_CUR) == -1) {
	  error_fatal("fseek", NULL);
	}
	break;
      }
      /* check if seq buffer is big enough */
      if ((n + strlen(p)) > n_tot) {
	n_tot += BUFFLEN;
	if ((seq = (char *)realloc(seq, n_tot*sizeof(char))) == NULL) {
	  error_fatal("Reallocating seq", NULL);
	}
	q = seq + n;
      }

      while(*p && *p != '\n') {
	if(!isascii((int)*p)) {
	  error_fatal("sequence", "contain non ascii characters");
	}
	if(isspace((int)*p)) {
	  p++;
	  continue;
	}
	n++;
	*q++ = *p++;
      }
    }
    *q = '\0';

    if(n == 0){
      error_fatal(id, "empty sequence");
    }

    /* correct id so that it does not contain \/ */
    if((q = strrchr(id, '/')) != NULL) {
      p = id;
      while(*q) {
	*p++ = *++q;
      }
      *p = '\0';
    }

    /* easiest way to unhide some files, dot is switched to _ */
    if(*id == '.') {
      *id = '_';
    }
    seq_holder->id = id;
    seq_holder->comment = comment;
    seq_holder->seq = seq;
    seq_holder->size = n;
  }
  else {
    seq_holder->id = NULL;
    seq_holder->comment = NULL;
    seq_holder->seq = NULL;
    seq_holder->size = 0;
  }

  free(BUFF);
  return state;
}


void free_seq(seq_t *seq) {
  free(seq->id);
  free(seq->comment);
  free(seq->seq);
}

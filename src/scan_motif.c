/*   File:                /home/edeveaud/Work/apoptose/src/scan_motif.c
 *   Author:              Eric Deveaud <edeveaud@pasteur.fr>
 *   Date first visited:  "Thu Apr 25 2002"
 *   Time-stamp:          "Tue Jun 11 2002"
 *   Development Stage :  Under construction / release
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include <stdio.h>

#ifdef STDC_HEADERS
#include <stdlib.h>
#include <string.h>
#endif

#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>
#endif

#ifdef HAVE_REGCOMP
#include <regex.h>
#endif

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#include "error.h"
#include "scan_motif.h"


/* private struct to hold the matching motif and position */
typedef struct match_S {
  char *motif;
  int pos;
} match_t;


static int new_match (regex_t *preg, char *seq, match_t *ret) ;


int search_motifs (param_t params, seq_t *seq_holder, mot_t *mot, int n) {

  FILE *OUT;
  match_t  ret;			/* the temporary holder for a found motif */
  motif_t *res_holder, *res;   	/* store the motif search results */
  int i, m, l;
  int pos, start, stop, len;
  int min, max, old_stop;
  int overlap;
  int display;
  char *seq, *s;
  int descente;
  char *tmp_seq;

  seq = seq_holder->seq;
  overlap = params.overlap;
  pos = 0;
  display = 0;
  m = 0;
  l = seq_holder->size;

  /* initialisation of reg-exp search*/
  ret.motif=NULL;
  ret.pos=0;

  OUT = params.OUT;

  /* allocating a result buffer for the search */
  if ((res_holder = malloc(sizeof(motif_t)*(n))) == NULL) {
    error_fatal("memory: res" , NULL);
  }

   /* allocating  tmp_seq for included motifs search */
  if ((tmp_seq = malloc(sizeof(char) * (l +1))) == NULL) {
    error_fatal("memory: tmp_seq" , NULL);
  }

  res = res_holder;
  for (i = 0; i < n; i++) {
    res[i].motif = NULL;
    res[i].start = 0;
    res[i].stop = 0;
    res[i].len = 0;
  }

  s = seq;
  descente = 0;
  while(m != -1) {
    if (descente == 1) {
      tmp_seq = strncpy(tmp_seq, seq+res[m].start, res[m].stop-res[m].start);
      *(tmp_seq+res[m].stop-res[m].start) = '\0';
    }
    i = new_match(mot[m].reg, s, &ret);
    if (descente == 1) {
      i = new_match(mot[m].reg, tmp_seq, &ret);
    }
    if (i != -1){	        /* motif found */
      if(descente == 1){
	start = res[m].start;
      }
      else {
	start = i + s - seq;	/* 0 indexed */
      }
      len = strlen(ret.motif);
      stop = start + len - 1;   /* 0 indexed */

      min = mot[m].min;
      max = mot[m].max;

      if(m == 0) old_stop =  0 ; /* 0 indexed */
      else old_stop = res[m-1].stop;

      /* check the minimal distance constraints */
      if(min != 0) {

	if(start - old_stop <= min) {
	  /* shift sequence from 1 */
	  s = seq + start + 1;
	  continue;
	}
      }
      /* check the maximal distance constraints */
      if(max != 0) {
	if(start - old_stop >= max) {
	  /* shift sequence from 1 */
	  s = seq + start + 1;
	  continue;
	}
      }

      /* update the last position where motif was found */
      mot[m].pos = start;

      /* fill results holder */
      if ( res[m].motif != NULL) free(res[m].motif);

      if((res[m].motif = malloc(sizeof(char)*(len+1+1))) == NULL) {
	error_fatal("memory: res.motif", NULL);
      }
      res[m].motif = strncpy(res[m].motif, ret.motif, len+1);
      res[m].start = start;
      res[m].stop = stop;
      res[m].len = len;

      /* shift sequence to check next motif */
      if (overlap == 1) {
	s = seq + start + 1;
      }
      else {
	s = seq + stop + 1;
      }
      /* check if there is some more motifs to find */
      if (m < n-1) {
	/* get next motif */
	descente  = 0;
	m++;
      }
      else {
	(void)fprintf(OUT, "%s ", seq_holder->id);
	for(i = 0; i < n-1; i++){
	  //(void)fprintf(OUT, "/%s/ %s %i %i ",
	  //	params.pattern,
	  //	res[i].motif, res[i].start + 1, res[i].stop + 1);

	  (void)fprintf(OUT, "%s %i %i ",
	                res[i].motif, res[i].start + 1, res[i].stop + 1);
	}
	(void)fprintf(OUT, "%s %i %i\n",
		      res[i].motif, res[i].start + 1, res[i].stop + 1);
	display ++;
	/* for each motif check if there are some included motifs
	 * as regexp return the longuest match, we should cut it from
	 * the end
	 * example motif == [RK]-X(0,1)-V-X-[FW]
	 * seq = pvlvKVVFFasm
	 * KVVFF match [RK]-X(n=1)-V-X-[FW]
	 * but include KVVF that match [RK]-X(n=0)-V-X-[FW
	 */

      }
    }
    else { /* motif not found */
      if(descente == 1) {
	m--;
	descente = 0;
	s = seq + mot[m].pos + 1 ;
	pos =  mot[m].pos +1;

      }
      else {
	descente = 1;
      }
    }

  } /* end while */

  /* some memory cleanning */
  free(ret.motif);
  for(i=0; i<n; i++) {
    free(res_holder[i].motif);
  }
  free(res_holder);
  free(tmp_seq);

  return display;
}


static int new_match (regex_t *preg, char *seq, match_t *ret) {

  char *start;
  char *p;
  int  count;
  int i;
  regmatch_t pmatch;

  start = seq;

  count = regexec(preg, start, 1, &pmatch, 0);
  if( count == REG_NOMATCH ) {
    return -1;
  }

  i = (pmatch.rm_eo - pmatch.rm_so);
  if(ret->motif != NULL) free(ret->motif);
  if((ret->motif = malloc(sizeof(char)*(i+1))) == NULL) {
    error_fatal("memory: new_match", NULL);
  }

  p = ret->motif;
  for (i = pmatch.rm_so; i < pmatch.rm_eo; i++){
    *p++ = seq[i];
  }
  *p++ = '\0';

  ret->pos = pmatch.rm_so;

  return pmatch.rm_so;
}

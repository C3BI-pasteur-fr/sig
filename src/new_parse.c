/*   File:                /home/edeveaud/Work/apoptose/src/new_parse.c
 *   Author:              Eric Deveaud <edeveaud@pasteur.fr>
 *   Date first visited:  "Tue May 07 2002"
 *   Time-stamp:          "Mon Jun 03 2002"
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

#ifdef HAVE_CTYPE_H
#include <ctype.h>
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

#define MOTIF		0
#define SEP		1
#define PAR_OPEN	2
#define PAR_CLOSE	3

static int stock_motif(char *pat, mot_t *mot, int n) ;
static int stock_dist(char *pat, mot_t *mot, int n) ;


static int stock_motif(char *pat, mot_t *mot, int n) {

  int i;
  char *p, *q, *res;
  char message[128];

  p = pat;
  i = strlen(pat);

  if(i == 0){
    error_fatal("pattern", "should contain at least one character");
  }

  /* allocating result buffer: due to prosite syntax, it could not be
   * larger than the original pattern */
  if((res = malloc(sizeof(char)*(i + 1))) == NULL) {
    error_fatal("memory: stock_pattern", NULL);
  }

  q = res;

  /* rewritte pattern skipping the - separator, and replacing prosite
   * location specification to regular expression ones */
  while (*p){

    switch(*p) {
    case '-':      break;

    case '{':
    case '[':      *q++ = '[';      break;
    case '}':
    case ']':      *q++ = ']';      break;
    case '(':	   *q++ = '{'; 	    break;
    case ')':      *q++ = '}';      break;
    case '>':      *q++ = '$';      break;
    case '<':      *q++ = '^';      break;
    case 'x':      *q++ = '.';      break;
    case 'X':	   *q++ = '.';	    break;
    default:      *q++ = *p;
    }
    p++;
  }
  *q++ = '\0';

  /* store the human readable motif for output reasons */
  if((mot[n].motif = malloc(sizeof(char)*(i + 1))) == NULL) {
    error_fatal("memory: stock_motif1", NULL);
  }
  mot[n].motif = strncpy(mot[n].motif, pat, i+1);

  /* compile and store the associated regex pattern */
  if((mot[n].reg = malloc(sizeof(regex_t))) == NULL) {
    error_fatal("memory: stock_motif2", NULL);
  }

  if((i = regcomp(mot[n].reg, res, REG_ICASE | REG_EXTENDED)) != 0) {
    if (regerror(i, mot[n].reg, message, 128)) {
      error_fatal("Translating RE:", message);
    }
    else {
      error_fatal("Translating RE:", res);
    }
  }

  free(res);
  return i;
}

static int stock_dist(char *dist, mot_t *mot, int n) {

  char *p, *q, *BUFF;
  int min, max;

  min = max = 0;

  if(strlen(dist) != 0) {
    if((BUFF = malloc(sizeof(char)*(BUFFLEN+1))) == NULL) {
      error_fatal("memory: stock_dist", NULL);
    }
    q = BUFF;
    p = dist;
    while (*p) {

      switch(*p) {
      case '(':
      case '{': break;

      case ',':
      case '.':  min = atoi(BUFF); q = BUFF; break;

      case '}':
      case ')':  max = atoi(BUFF); break;

      default : *q++ = *p;

      }
      p++;
    }
    free (BUFF);
  }

  mot[n].min = min;
  mot[n].max = max;
  return 0;
}

void free_mot_t (mot_t *mot, int n) {

  int i;

  for (i = 0; i <n; i++) {
    free(mot[i].motif);
    /*freeing regex_t internals */
    regfree(mot[i].reg);
    /* freeing allocated regex_t */
    free(mot[i].reg);
  }
  free(mot);
}


int check_pattern (char *pattern, mot_t **mot) {

  char *BUFF, *p, *q;
  int i, j, len, n;
  int state;
  mot_t *res;

  i = j = 0;
  n = 0;
  len = strlen(pattern);

  if(len == 0){
    error_fatal("pattern", "should contain at least one character");
  }
  /* allocate the motif holder for the first motif*/
  if ((res = malloc(sizeof(mot_t))) == NULL) {
    error_fatal("memory: check_pattern", NULL);
  }

  /* allocate temporary working buffer */
  if((BUFF = malloc(sizeof(char)*(len+1))) == NULL) {
    error_fatal("memory: tmp_pattern", NULL);
  }

  state = MOTIF;
  i = 0;
  p = pattern;
  q = BUFF;

  while (*p) {

    switch (state) {

    case MOTIF:
      if( *p != ' ') {
	*q++ = *p;
      }
      else {
	*q++ = '\0';
	if(stock_motif(BUFF, res, i) != 0) {
	  error_fatal(BUFF, "motif incorrect");
	}
	state = SEP;
	i++;
	q = BUFF;
      }
      break;

    case SEP:
      if(*p != '(') {
	error_fatal(pattern, "pattern syntax error");
      }
      else {
	*q++ = *p;
	state = PAR_OPEN;
      }
      break;

    case PAR_OPEN:
      if(*p != ')') {
	if (isdigit((int)*p) != 0 || *p == ',')  {
	  *q++ = *p;
	}
	else {
	  error_fatal(pattern, "pattern syntax error");
	}
      }
      else {
	*q++ = *p;
	*q = '\0';
	/* as there is a distance, they will be a next motif */
	/* reallocate the motif holder for this one */
	if ((res = realloc(res, sizeof(mot_t)*(i+1) )) == NULL) {
	  error_fatal("memory: res reallocation", NULL);
	}
	stock_dist(BUFF, res, i);
	state = PAR_CLOSE;
      }
      break;

    case PAR_CLOSE:
      if(*p != ' ') {
	error_fatal(pattern, "pattern syntax error");
      }
      else {
	q = BUFF;
	state = MOTIF;
      }
      break;

    default:
      error_fatal(pattern, "pattern syntax error");
      }

      p++;
  }
  if(state!= MOTIF) {
    error_fatal(pattern, "pattern syntax error");
  }
  else {
    *q = '\0';
    if(stock_motif(BUFF, res, i) != 0) {
        error_fatal(BUFF, "motif incorrect");
    }
    i++;
  }

  free(BUFF);

  *mot = res;
  return i;
}





char *read_pattern(FILE *IN) {

  int i,l;
  int len, l_res;
  char *BUFF, *res, *p, *q;

  len = l_res = BUFFLEN;

  if((BUFF = malloc(sizeof(char) * (len + 1))) == NULL){
    error_fatal("memory: read_pattern1", NULL);
  }
  if((res = malloc(sizeof(char) * (l_res + 1))) == NULL){
    error_fatal("memory: read_pattern1", NULL);
  }
  *res='\0'; 
  l = 0;
  p = res;
  while ((fgets(BUFF, len, IN)) != NULL) {
    
    /* line is not compety read */
    if ((strrchr(BUFF, '\n') == NULL) && (feof(IN) !=1)) {
      len += BUFFLEN;
      if((BUFF = realloc(BUFF, sizeof(char)*(len + 1))) == NULL) {
	error_fatal("memory realloc", NULL);
      }
      if(fseek(IN, -1 *(long)strlen(BUFF), SEEK_CUR)){
	error_fatal("fseek", "fseek");
      }
      continue;
    }
    
    /*skip comment lines and empty lines*/
    if(*BUFF == '#'|| *BUFF == '\n') continue;
    
    /* concatenate read lines */
    
    i = strlen(BUFF);
    l = strlen(res);

    if(i + l >= l_res) {
      l_res += i;
      if((res = realloc(res, sizeof(char)*(l_res + 1))) == NULL) {
	error_fatal("memory realloc", NULL);
      }

    }
    p = res + l;
    q = BUFF;
    while(*q && *q != '\n') {
      *p++ = *q++;
      l++;
    }

    *p = ' ';

  }
  *p = '\0';

  free(BUFF);

  if(l == 0) {
    free(res);
    return NULL;
  }
  return res;
}

/*   File:                /home/edeveaud/Work/apoptose/src/new_parse.h
 *   Author:              Eric Deveaud edeveaud@pasteur.fr
 *   Date first visited:  "Tue May 07 2002"
 *   Time-stamp:          "Fri May 17 2002"
 *   Development Stage :  Under construction / distrib / final
 */


#ifndef __NEW_PARSE_H_
#define __NEW_PARSE_H_

#include "scan_motif.h"

#define BUFFLEN 100

int check_pattern (char *pattern, mot_t **mot) ;

void free_mot_t (mot_t *mot, int n) ;

char *read_pattern(FILE *IN) ;

#endif /* __NEW_PARSE_H_ */

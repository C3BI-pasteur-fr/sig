/*   File:                /home/edeveaud/Work/apoptose/src/main.c
 *   Author:              Eric Deveaud <edeveaud@pasteur.fr>
 *   Date first visited:  "Tue Apr 23 2002"
 *   Time-stamp:          "Thu Jul 10 2003"
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

#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>
#endif

#ifdef HAVE_REGCOMP
#include <regex.h>
#endif

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#include "params.h"
#include "error.h"
#include "scan_motif.h"
#include "seq-reader.h"
#include "new_parse.h"


static char *prog;

static void usage(char *prog);

static void process(FILE *IN, param_t params, mot_t *mot, mot_t * rev, int n) ;

int main(int argc, char **argv) {

  /*variables and initialisation */
  int i, n;
  mot_t *mot;	/* store the motifs definition to find */
  mot_t *rev;	/* only used if -r option is given */
  char *p;
  char  *outfile, *pat_file;
  FILE *IN;
  param_t params;

  /* default initialisation */
  outfile = NULL;
  IN = stdin;
  params.OUT = stdout;
  params.pattern = NULL;
  params.outfile = NULL;
  params.overlap = 0;
  params.reverse = 0;
  pat_file = NULL;
  rev = NULL;

  /* get progname */
  prog = argv[0];
  if((p = strrchr(prog, '/')) != NULL)
    prog = ++p;

  /* check syntax option on command line */
  if (argc == 1) {
    usage(prog);
  }

  i = 0;
  while((i = getopt(argc, argv, "f:hio:p:rv")) != -1) {

    switch(i) {

    case 'f':
      pat_file = optarg;
      break;

    case 'h':
      usage(prog);
      return EXIT_SUCCESS;

    case 'i':
      params.overlap = 1;
      break;

    case 'o':
      outfile = optarg;
      break;

    case 'p':
      params.pattern = optarg;
      break;

    case 'r':
      params.reverse = 1;
      break;

    case 'v':
      (void)fprintf(stdout, "%s (%s %s)\n", prog, PACKAGE, VERSION);
      return EXIT_SUCCESS;

    default:
      usage(prog);
      return EXIT_FAILURE;

    }
  }

  /* checking syntax arguments validity */
  if (params.pattern == NULL && pat_file == NULL)
    error_fatal(prog, "pattern is mandatory");

  if (params.pattern != NULL && pat_file != NULL)
    error_fatal(prog, "options -f and -p are mutualy exclusive");

  /* check for pattern definition file, and read pattern from */
  if(pat_file != NULL ){
    if((IN = fopen(pat_file, "r")) == NULL) {
      error_fatal(pat_file, NULL);
    }
    if((params.pattern = read_pattern(IN)) == NULL) {
      error_fatal("pattern", "should contain at least one character");
    }
  }
  if(fclose(IN) != 0) {
    error_fatal(pat_file, NULL);
  }

  /* check output file availability */
  if (outfile != NULL && (params.OUT = fopen(outfile, "w")) == NULL) {
    error_fatal(outfile, NULL);
  }
  params.outfile = outfile;

  n = check_pattern(params.pattern, &mot);
  /* should we reverse the search order*/
  if(params.reverse == 1 ) {
    if ( n > 1) {
      if((rev = malloc(sizeof(mot_t)*n)) == NULL) {
	error_fatal("memory: allocating reverse", NULL);
      }
      /* motif should be interverted */
      for (i = 0; i<n; i++) {
	rev[i].motif = mot[n-i-1].motif;
	rev[i].reg = mot[n-i-1].reg;
      }
      /* but not the distance constraints */
      for (i = 0; i<n; i++) {
	rev[i].min = mot[i].min;
	rev[i].max = mot[i].max;
      }
    }
    /* not disable reverse search, when there's just ONE motif */
    /* we will duplicate the results */
    else {
      params.reverse = 0;
      error_warn("are you stupid ?", "your pattern contains just ONE motif");
    }
  }

  /* getting the sequence */
  for (i = optind; i <argc; i++) {
    if (*argv[i] != '-' && (IN = fopen(argv[i], "r")) == NULL) {
      error_fatal(argv[i], NULL);
    }
    process(IN, params, mot, rev, n);

    if(IN != stdin && fclose(IN) != 0) {
      error_fatal(argv[i], NULL);
    }
  }

  /* memory cleaning */
  free_mot_t(mot, n);
  if (rev != NULL) free(rev);
  if (pat_file  != NULL) free(params.pattern);


  /* closing files */
  if (params.OUT != stdout && fclose(params.OUT) != 0) {
    error_fatal(params.outfile, NULL);
  }

  return 0;
}

static void process(FILE *IN, param_t params, mot_t *mot, mot_t * rev, int n) {

  seq_t seq_holder;
  int j;

  while (read_seq(IN, &seq_holder) != 0) {

    if(rev == NULL) {
      j = search_motifs(params, &seq_holder, mot, n);
    }
    else {
      j = search_motifs(params, &seq_holder, mot, n);
      j = search_motifs(params, &seq_holder, rev, n);
    }
  free_seq(&seq_holder);
  }
}

static void usage(char *prog_name) {

   FILE *PERR = stderr;

   (void)fprintf(PERR, "usage: %s [options] <file>\n", prog_name);
   (void)fprintf(PERR, "  -f <file>      ... Read pattern from file <file>\n");
   (void)fprintf(PERR, "  -h             ... Print this message and exit.\n");
   (void)fprintf(PERR, "  -i             ... Allow overlaping motifs\n");
   (void)fprintf(PERR, "  -p <pattern>   ... Specify the pattern to search\n");
   (void)fprintf(PERR, "  -o <file>      ... Place the output into <file>.\n");
   (void)fprintf(PERR, "  -r             ... Perform a second pass with reverse motif order\n");
   (void)fprintf(PERR, "  -v             ... Print version number and exit.\n");
}

#include <stdio.h>
#include <ctype.h>
#include <stdarg.h>
#include <string.h>

#define SCREAD_MAXB 4096

/*!
 * \file scread.c
 *
 * \brief Codi de scread()
 */

/*!
 * \brief (Skip-comment-and-read) Llegeix un stream a base de
 * <tt>fscanf(fp,"%s",...)</tt>, i intenta llegir cada paraula a
 * partir de la llista de formats-adreces que se li passen als
 * arguments de després de \c fp i \c nitems. Se salta comentaris,
 * considerant comentari el que hi ha entre un # i final de línia.
 */
int scread (FILE *fp, int nitems, ...) {
   va_list ap; char *fmt; void *where;
   char buf[SCREAD_MAXB];
   int c, r=0, i;
   va_start(ap, nitems);
   while (r<nitems) {
   /* Skip blanks */
      do {
         if ((c=fgetc(fp))==EOF) return(r);
         if (c=='!' || c=='#')
            do { if ((c=fgetc(fp))==EOF) return(r); } while (c!='\n');
      } while (isspace(c));
   /* Store token in buf */
      i=0;
      do {
         if (i==SCREAD_MAXB-1)
            fprintf(stdout,
                  "scread()::warning : buf[] capacity (%d) exceded!!\n",
               SCREAD_MAXB);
         if (i<SCREAD_MAXB-1) buf[i++]=c;
         c=fgetc(fp);
         if (c=='!' || c=='#') do { c=fgetc(fp); } while (c!='\n' && c!=EOF);
      } while (!isspace(c));
      buf[i]='\0';
      if (c=='!' || c=='#') do { c=fgetc(fp); } while (c!='\n' && c!=EOF);
   /* Read token */
      fmt=va_arg(ap, char *); where=va_arg(ap, void *);
   /* FORTRAN support: If format is "%lf", switch D/d into E/e */
      if (strcmp(fmt,"%lf")==0) {
         for (i=0; i<strlen(buf); i++) {
            if (buf[i]=='d') buf[i]='e';
            if (buf[i]=='D') buf[i]='E';
         }
      }
      if (sscanf(buf, fmt, where)!=1) return(r);
      r++;
   }
/* All done! */
   return(r);
}

/*!
 * \brief Com scread, però torna el nombre de newlines que ha
 * travessat dins nnl;
 */
int screadcnl (int *nnl, FILE *fp, int nitems, ...) {
   va_list ap; char *fmt; void *where;
   static char buf[SCREAD_MAXB];
   int c, r=0, i;
   va_start(ap, nitems);
   while (r<nitems) {
   /* Skip blanks */
      do {
         if ((c=fgetc(fp))==EOF) return(r);
         if (c=='!' || c=='#')
            do { if ((c=fgetc(fp))==EOF) return(r); } while (c!='\n');
	 if (c=='\n') (*nnl)++;
      } while (isspace(c));
   /* Store token in buf */
      i=0;
      do {
         if (i==SCREAD_MAXB-1)
            fprintf(stdout,
                  "scread()::warning : buf[] capacity (%d) exceded!!\n",
               SCREAD_MAXB);
         if (i<SCREAD_MAXB-1) buf[i++]=c;
         c=fgetc(fp);
         if (c=='!' || c=='#') do { c=fgetc(fp); } while (c!='\n' && c!=EOF);
      } while (!isspace(c));
      if (c=='\n') (*nnl)++;
      buf[i]='\0';
      if (c=='!' || c=='#') do { c=fgetc(fp); } while (c!='\n' && c!=EOF);
   /* Read token */
      fmt=va_arg(ap, char *); where=va_arg(ap, void *);
   /* FORTRAN support: If format is "%lf", switch D/d into E/e */
      if (strcmp(fmt,"%lf")==0) {
         for (i=0; i<strlen(buf); i++) {
            if (buf[i]=='d') buf[i]='e';
            if (buf[i]=='D') buf[i]='E';
         }
      }
      if (sscanf(buf, fmt, where)!=1) return(r);
      r++;
   }
/* All done! */
   return(r);
}

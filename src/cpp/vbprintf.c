#include <stdio.h>
#include <stdarg.h>
// #ifdef __GNUC__
// #include "printf.h"
// #endif /* __GNUC__ */

FILE *vbprintf_log2=NULL;

/*!
 * \file
 *    vbprintf.c
 *
 * \brief
 *    Escriptura de missatges controlats per un índex de verbositat.
 */

/*!
 * \brief
 *    Funciona com <tt>fprintf(stderr, ...)</tt>, però escriu o no
 *    segons un índex de verbositat.
 *
 * El missatge nomes s'escriu si l'índex de verbositat es mes gran o
 * igual que el nivell de verbositat de missatge. Quan mes gran es el
 * nivell de verbositat del missatge, mes costa d'escriure'l.
 *
 * Arguments:
 * \arg
 *    \c vbl (e): nivell de verbositat del missatge,
 * \arg
 *    \c ivb (e): índex de verbositat
 * \arg
 *    resta d'arguments: igual que \c printf.
 *
 * El valor de retorn és com a \c printf.
 */

#define MAXARGTYPES 50

int vbprintf (int vbl, int ivb, const char *fmt, ...) {
   if (ivb>=vbl) {
      va_list ap; int iret;
      va_start(ap, fmt);
      if (vbprintf_log2!=NULL) {
         iret=vfprintf(vbprintf_log2,fmt,ap);
         vfprintf(stderr,fmt,ap);
      } else
         iret=vfprintf(stderr,fmt,ap);
// #ifdef __GNUC__
//       if (vbprintf_log!=NULL) {
//          int argtypes[MAXARGTYPES], n, i;
//          va_list copy;
//          va_copy(copy, ap);
//          n=parse_printf_format(fmt, MAXARGTYPES, argtypes);
//          if (n>MAXARGTYPES) n=MAXARGTYPES;
//       /* Clear strings that begin with an escape code */
//          for (i=0; i<n; i++) {
//             char *s=va_arg(ap,char *);
//             if ( (argtypes[i] & ~PA_FLAG_MASK) == PA_STRING)
//                if (s[0]=='\033') s[0]='\0';
//          }
//       /* Print in log file */
//          vfprintf(vbprintf_log,fmt,copy);
//       }
// #endif /* __GNUC__ */
      va_end(ap);
      return(iret);
   } else return(0);
}

#undef MAXARGTYPES

#ifdef TEST

int main (int argc, char *argv[]) {
   int ivb, vbl;
   if (argc!=4
	 || sscanf(argv[1], "%d", &ivb)!=1
	 || sscanf(argv[2], "%d", &vbl)!=1
      ) {
      fprintf(stderr, "vbfprintf_t.exe ivb vbl msg\n");
      return(-1);
   }
   vbprintf(vbl,ivb,argv[3]);
   printf("\n");
   return(0);
}

#endif /* TEST */

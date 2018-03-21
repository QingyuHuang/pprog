#ifndef HAVE_KOMPLEX_H /* this is necessary for multiple includes */

struct komplex {double re; double im;};
typedef struct komplex komplex;

void    komplex_print    (char* s, komplex z);   /* prints string s and then komplex z */
void    komplex_set      (komplex* z, double x, double y);   /* z is set to x+i*y */
komplex komplex_new      (double x, double y);   /* returns x+i*y */
komplex komplex_add      (komplex a, komplex b); /* returns a+b */
komplex komplex_sub      (komplex a, komplex b); /* returns a-b */
int     komplex_equal    (komplex a, komplex b); /* returns 1 if equal, 0 otherwise */
komplex komplex_mul      (komplex a, komplex b); /* returns a*b */
komplex komplex_div      (komplex a, komplex b); /* returns a/b */
komplex komplex_conjugate(komplex z);            /* returns complex conjugate */
komplex komplex_abs      (komplex z);         /* returns absolute value */
komplex komplex_exp      (komplex z);         /* returns e^z */
komplex komplex_sin      (komplex z);         /* returns sin(z) */
komplex komplex_cos      (komplex z);         /* returns cos(z) */
komplex komplex_sqrt     (komplex z);         /* returns sqrt(z) */

#define HAVE_KOMPLEX_H
#endif

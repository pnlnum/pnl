/*
 * If BASE is undefined we use function names like pnl_name()
 *  and assume that we are using doubles.
 *
 * If BASE is defined we used function names like pnl_BASE_name()
 * and use BASE as the base datatype
*/

#if defined BASE_DOUBLE
#define ORDERED true
#define BASE double
#define SHORT
#define LSHORT
#define ATOMIC double
#define BASE_TYPE DOUBLE
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%lf"
#define OUT_FORMAT "%f"
#define IN_PUT_FORMAT(a) (a)
#define OUT_PUT_FORMAT(a) a
#define ATOMIC_IO ATOMIC
#define ZERO 0.0
#define ONE 1.0
#define BASE_EPSILON PNL_DBL_EPSILON
#define EQ(a,b) (a)==(b)
#define NEQ(a,b) (a)!=(b)
#define PLUS(a,b) (a)+(b)
#define MINUS(a,b) (a)-(b)
#define MULT(a,b) (a)*(b)
#define DIV(a,b) (a)/(b)
#define PLUSEQ(a,b) (a)+=(b)
#define MINUSEQ(a,b) (a)-=(b)
#define MULTEQ(a,b) (a)*=(b)
#define DIVEQ(a,b) (a)/=(b)
#define INV(a) 1/(a)
#define SQUARE_NORM(a) (a)*(a)
#define NORMONE(a) fabs(a)
#define PNL_C2F(f) C2F(d##f)
#define MALLOC_BASE MALLOC_DOUBLE


#elif defined BASE_PNL_COMPLEX
#define BASE dcomplex
#define SHORT Complex
#define LSHORT _complex
#define SHORT_REAL
#define ATOMIC double
#define BASE_TYPE COMPLEX
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%lf + %lfi"
#define OUT_FORMAT "%f + %fi"
#define IN_PUT_FORMAT(a) &((*(a)).r),&((*(a)).i)
#define OUT_PUT_FORMAT(a) (a).r,(a).i
#define ATOMIC_IO ATOMIC
#define ZERO CZERO
#define ONE CONE
#define BASE_EPSILON PNL_DBL_EPSILON
#define EQ(a,b) ((a).r==(b).r)&&((a).i==(b).i)
#define NEQ(a,b) ((a).r!=(b).r)||((a).i!=(b).i)
#define PLUS(a,b) Cadd(a,b)
#define MINUS(a,b) Csub(a,b)
#define MULT(a,b) Cmul(a,b)
#define DIV(a,b) Cdiv(a,b)
#define PLUSEQ(a,b) a=Cadd(a,b)
#define MINUSEQ(a,b) a=Csub(a,b)
#define MULTEQ(a,b) a=Cmul(a,b)
#define DIVEQ(a,b) a=Cdiv(a,b)
#define INV(a) Cinv(a)
#define SQUARE_NORM(a) Csqr_norm(a)
#define NORMONE(a) Cabs(a)
#define PNL_C2F(f) C2F(z##f)
#define MALLOC_BASE MALLOC_COMPLEX

#elif defined BASE_UINT
#define ORDERED true
#define BASE uint
#define SHORT Uint
#define LSHORT _uint
#define ATOMIC unsigned int
#define BASE_TYPE UINT
#define MULTIPLICITY 1
#define IN_FORMAT "%u"
#define OUT_FORMAT "%u"
#define IN_PUT_FORMAT(a) (a)
#define OUT_PUT_FORMAT(a) (a)
#define ATOMIC_IO ATOMIC
#define ZERO 0U
#define ONE 1U
#define EQ(a,b) (a)==(b)
#define NEQ(a,b) (a)!=(b)
#define PLUS(a,b) (a)+(b)
#define MINUS(a,b) (a)>(b)?(a)-(b):0
#define MULT(a,b) (a)*(b)
#define DIV(a,b) (a)/(b)
#define PLUSEQ(a,b) (a)+=(b)
#define MINUSEQ(a,b) (a)-=(b)
#define MULTEQ(a,b) (a)*=(b)
#define DIVEQ(a,b) (a)/=(b)
#define INV(a) 1
#define SQUARE_NORM(a) (a)*(a)
#define NORMONE(a) (a)
#define UNSIGNED 1

#elif defined BASE_INT
#define ORDERED true
#define BASE int
#define SHORT Int
#define LSHORT _int
#define ATOMIC int
#define BASE_TYPE INT
#define MULTIPLICITY 1
#define IN_FORMAT "%d"
#define OUT_FORMAT "%d"
#define IN_PUT_FORMAT(a) (a)
#define OUT_PUT_FORMAT(a) (a)
#define ATOMIC_IO ATOMIC
#define ZERO 0
#define ONE 1
#define EQ(a,b) (a)==(b)
#define NEQ(a,b) (a)!=(b)
#define PLUS(a,b) (a)+(b)
#define MINUS(a,b) (a)-(b)
#define MULT(a,b) (a)*(b)
#define DIV(a,b) (a)/(b)
#define PLUSEQ(a,b) (a)+=(b)
#define MINUSEQ(a,b) (a)-=(b)
#define MULTEQ(a,b) (a)*=(b)
#define DIVEQ(a,b) (a)/=(b)
#define INV(a) 1
#define SQUARE_NORM(a) (a)*(a)
#define NORMONE(a) abs(a)

#elif defined MY_NEW_BASE

#else
#error unknown BASE_ directive
#endif
 
#define CONCAT2x(a,b) a ## b
#define CONCAT2(a,b) CONCAT2x(a,b)
#define CONCAT3x(a,b,c) a ## b ## _ ## c
#define CONCAT3(a,b,c) CONCAT3x(a,b,c)
#define CONCAT4x(a,b,c,d) a ## _ ## b ## _ ## c ## _ ## d
#define CONCAT4(a,b,c,d) CONCAT4x(a,b,c,d)

#define FUNCTION(dir,name) CONCAT3(dir,LSHORT,name)

#define TYPE(dir) CONCAT2(dir,SHORT) /*dir */

#define STRING(x) #x
#define EXPAND(x) STRING(x)
#define NAME(x) EXPAND(TYPE(x))

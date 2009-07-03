/*
 * The idea of using such template files is owed to the GSL and has been
 * adapted to match the need of PNL
 */


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
#define MULTIPLICITY 1
#define FP 1
#define IN_FORMAT "%lf"
#define OUT_FORMAT "%f"
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
#define INV(a) 1/(a)
#define SQUARE_NORM(a) (a)*(a)
#define NORMONE(a) fabs(a)


#elif defined BASE_PNL_COMPLEX
#define BASE fcomplex
#define SHORT Complex
#define LSHORT _complex
#define SHORT_REAL
#define ATOMIC double
#define MULTIPLICITY 2
#define FP 1
#define IN_FORMAT "%lf"
#define OUT_FORMAT "%7.4f + i * %7.4f"
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
#define INV(a) Cinv(a)
#define SQUARE_NORM(a) Csqr_norm(a)
#define NORMONE(a) Cabs(a)

#elif defined BASE_UINT
#define ORDERED true
#define BASE uint
#define SHORT Uint
#define LSHORT _uint
#define ATOMIC unsigned int
#define MULTIPLICITY 1
#define IN_FORMAT "%u"
#define OUT_FORMAT "%u"
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
#define MULTIPLICITY 1
#define IN_FORMAT "%d"
#define OUT_FORMAT "%d"
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
#define INV(a) 1
#define SQUARE_NORM(a) (a)*(a)
#define NORMONE(a) abs(a)

#elif defined BASE_CONTAIN
#define BASE 
#define KEY int
#define VALUE double
#define CONTAIN PnlContains
#define NODE PnlNode
#define NODE_SHORT node
#define CONTAIN_SHORT contains 
#define SHORT  
#define LSHORT sort_list
#define IN_FORMAT "%lf"
#define OUT_FORMAT "%f"
#define OUT_PUT_FORMAT(a) a

#elif defined BASE_SPARSE_POINT
#define BASE 
#define KEY PnlVectUint *
#define VALUE int
#define CONTAIN PnlSparsePoint
#define NODE PnlNodeSparsePoint
#define NODE_SHORT node_sparse_point
#define CONTAIN_SHORT sparse_point
#define SHORT SparsePoint  
#define LSHORT sort_list_sparse_point
#define IN_FORMAT "%lf"
#define OUT_FORMAT "%f"
#define OUT_PUT_FORMAT(a) a

#else
#error unknown BASE_ directive in source.h
#endif
 
#define CONCAT2x(a,b) a ## b
#define CONCAT2(a,b) CONCAT2x(a,b)
#define CONCAT3x(a,b,c) a ## b ## _ ## c
#define CONCAT3(a,b,c) CONCAT3x(a,b,c)
#define CONCAT4x(a,b,c,d) a ## _ ## b ## _ ## c ## _ ## d
#define CONCAT4(a,b,c,d) CONCAT4x(a,b,c,d)


#define USE_QUALIFIER = true
#ifndef USE_QUALIFIER
#define QUALIFIER
#endif

#define FUNCTION(dir,name) CONCAT3(dir,LSHORT,name)
#define FUNCTION_NODE(dir,name) CONCAT3(dir,NODE_SHORT,name)
#define FUNCTION_CONTAIN(dir,name) CONCAT3(dir,CONTAIN_SHORT,name)

#define TYPE(dir) CONCAT2(dir,SHORT) /*dir */
#define VIEW(dir,name) CONCAT2(dir,name)
#define QUALIFIED_TYPE(dir) QUALIFIER dir
#define QUALIFIED_VIEW(dir,name) CONCAT3(dir,QUALIFIER,name)

#define STRING(x) #x
#define EXPAND(x) STRING(x)
#define NAME(x) EXPAND(TYPE(x))

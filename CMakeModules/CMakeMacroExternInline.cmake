
MACRO(CMAKE_CHECK_EXTERN_INLINE)
    #
    # Checks whether the following construction is valid
    #
    # extern int foo (int);
    # extern inline int foo (int x) { return x+1;}
    # int foo (int x) { return x+1;}
    set (HAVE_INLINE false)

    file(WRITE "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/test_extern_inline.c"
        "extern int foo (int);
        extern inline int foo (int x) { return x+1;}
        int foo (int x) { return x+1;}
        int main (int arg, char *argv[]) { foo(1); return 0;}
        ")

    try_compile(C_HAS_INLINE 
        "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/"
        "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/test_extern_inline.c")

    if(C_HAS_INLINE)
        set(HAVE_INLINE true cache internal "have 'inline' definition")
        message (STATUS "Check extern inline...yes")
    else ()
        message (STATUS "Check extern inline...no")
    endif(C_HAS_INLINE)

ENDMACRO(CMAKE_CHECK_EXTERN_INLINE)

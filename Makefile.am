AM_CFLAGS = @GSL_CFLAGS@
## Makefile.am -- Process this file with automake to produce Makefile.in

bin_PROGRAMS = kemp.bin basu.bin pvalue_ab P-value2d pvalue_graph loginorm gamma_test
kemp_bin_SOURCES = kempbasu.c fbst-poisson.c integrate.c P-value.c pvalue.c matrix.c alpha.c matrix.h kempbasu.h params.h
kemp_bin_LDADD = @GSL_LIBS@ -lJudy @PTHREAD_LIBS@ -lm
kemp_bin_CFLAGS = -DKEMP @GSL_CFLAGS@ @PTHREAD_CFLAGS@

basu_bin_SOURCES = kempbasu.c fbst-poisson.c fbst-loginorm.c integrate.c pvalue.c alpha.c  matrix.c matrix.h kempbasu.h params.h
basu_bin_LDADD = @GSL_LIBS@ -lJudy @PTHREAD_LIBS@ -lm
basu_bin_CFLAGS = -DBASU @GSL_CFLAGS@ @PTHREAD_CFLAGS@

P_value2d_SOURCES = P-value2d.c fbst-poisson.c integrate.c P-value.c matrix.c matrix.h kempbasu.h params.h
P_value2d_LDADD = @GSL_LIBS@ -lJudy @PTHREAD_LIBS@ -lm

pvalue_ab_SOURCES = pvalue_ab.c alpha.c
pvalue_ab_LDADD = @GSL_LIBS@ -lm

pvalue_graph_SOURCES = pvalue_graph.c pvalue.c alpha.c
pvalue_graph_LDADD = @GSL_LIBS@ -lm

loginorm_SOURCES = loginorm.c matrix.c
loginorm_LDADD = @GSL_LIBS@ -lm

gamma_test_SOURCES = gamma_test.c
gamma_test_LDADD = @GSL_LIBS@ -lm

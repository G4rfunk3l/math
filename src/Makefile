CC=gcc
ifeq ($(shell uname), Linux)
	CHKFLAGS=-pthread -lcheck_pic -pthread -lrt -lm -lsubunit
else
	CHKFLAGS=-lcheck
endif

all: my_math.a

clean:
	@rm -rf *.o *.gcda *.gcno *.html *.css *.out *.a test gcov_report my_test.c

my_math.o: my_math.c
	gcc -o my_math.o $^ -c

my_math.a: my_math.o
	ar crs $@ $^

test:
	checkmk clean_mode=1 my_test.check > my_test.c
	$(CC) my_math.c -c -o my_math_test.o
	$(CC) my_test.c -c -o test_test.o
	$(CC) test_test.o my_math_test.o -o $@ $(CHKFLAGS)
	./$@

gcov_report:
	checkmk clean_mode=1 my_test.check > my_test.c
	$(CC) --coverage my_math.c -c -o my_math_coverage.o
	$(CC) my_test.c -c -o test_coverage.o
	$(CC) my_math_coverage.o test_coverage.o -o $@ $(CHKFLAGS) -lgcov
	./$@
	gcovr --html --html-details -o report.html
	
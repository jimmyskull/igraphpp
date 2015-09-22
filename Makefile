
.PHONY: test

all:

test:
	mkdir -p build
	cd build && cmake -DCMAKE_BUILD_TYPE=Debug .. && make --no-print-directory -j9 && valgrind --leak-check=full ./test/test_suite -d yes

format:
	clang-format --style=llvm -i test/*.cpp
	clang-format --style=llvm -i examples/*.cpp
	clang-format --style=llvm -i igraphpp/*.hpp

clean:
	$(RM) -r ./build/


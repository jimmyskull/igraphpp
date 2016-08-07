
.PHONY: test

all:

test:
	mkdir -p build
	cd build && cmake -DCMAKE_BUILD_TYPE=Debug -DBuildTest:Bool=On .. && $(MAKE) --no-print-directory
	valgrind ./build/test/test_suite -d yes

format:
	clang-format --style=llvm -i test/*.cpp
	clang-format --style=llvm -i examples/*.cpp
	clang-format --style=llvm -i igraphpp/*.hpp

clean:
	$(RM) -r ./build/


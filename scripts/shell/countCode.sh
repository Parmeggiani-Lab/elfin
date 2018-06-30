#!/bin/bash

find . -type f -name '*.c' -o -name '*.h' -o -name '*.cpp' -o -name '*.hpp' -o -name '*.sh' | grep -v 'venv\|json.hpp' | xargs wc -l

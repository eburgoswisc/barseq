#!/usr/bin/env bash
# Make temporary directory

set -e

CREATE_DUMP() {

if mkdir dump_test > /dev/null 2>&1; then
    cd dump_test
else
    cd dump_test
    rm -r *
fi

}

CLEAN_UP() {
    rm -r ../dump_test
}
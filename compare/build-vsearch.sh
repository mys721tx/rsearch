#!/bin/sh
# Builds vsearch from its source tree without autotools, for comparing
# rsearch's output with vsearch's.
#
# usage: compare/build-vsearch.sh [VSEARCH_SOURCE] [OUTPUT_DIR]
# defaults: ../vsearch next to this repository, and compare/vsearch-build
#
# The binary is written to OUTPUT_DIR/vsearch. The architecture and the
# operating system come from uname, and pick vsearch's backends as its
# Makefile.am does: arch/x86_64 (SSE2 and SSSE3), arch/aarch64 (NEON),
# arch/ppc64le (AltiVec) or arch/generic; os/macos, os/freebsd, or
# os/linux for Linux and the other Unix systems, with os/posix.
# VSEARCH_ARCH and VSEARCH_OS override uname -m and uname -s; DRY_RUN=1 prints
# what would be compiled and stops. Tested on macOS on Apple silicon.
set -eu

here=$(cd "$(dirname "$0")" && pwd)
source_dir=$(cd "${1:-$here/../../vsearch}" && pwd)/src
out=$(mkdir -p "${2:-$here/vsearch-build}" && cd "${2:-$here/vsearch-build}" && pwd)
arch=${VSEARCH_ARCH:-$(uname -m)}
os=${VSEARCH_OS:-$(uname -s)}
CXX=${CXX:-c++}
CC=${CC:-cc}
jobs=$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4)

# the architecture's backend: its sources, and the flags of those that
# need their own (x86_64 builds both SSE levels and picks one at run time)
arch_flags=""
case $arch in
  x86_64 | amd64)
    arch_sources="arch/x86_64/cpu_features.cpp"
    simd_sources="arch/x86_64/SSE2/increment_counters.cpp:-msse2 arch/x86_64/SSSE3/increment_counters.cpp:-mssse3" ;;
  arm64 | aarch64)
    arch_sources="arch/aarch64/increment_counters.cpp arch/generic/cpu_features.cpp"
    simd_sources="" ;;
  ppc64le)
    arch_sources="arch/ppc64le/increment_counters.cpp arch/generic/cpu_features.cpp"
    simd_sources=""
    arch_flags="-mcpu=powerpc64le -maltivec" ;;
  *)
    arch_sources="arch/generic/increment_counters.cpp arch/generic/cpu_features.cpp"
    simd_sources="" ;;
esac

# the operating system's backend; os/linux is vsearch's generic Unix one
libs="-lpthread"
case $os in
  Darwin) os_sources="os/posix/dynlib_loader.cc os/posix/system.cc os/macos/system_memory.cc os/macos/system_cpu.cc" ;;
  FreeBSD) os_sources="os/posix/dynlib_loader.cc os/posix/system.cc os/freebsd/system_memory.cc os/freebsd/system_cpu.cc" ;;
  MINGW* | MSYS* | CYGWIN* | Windows*)
    echo "build-vsearch.sh: Windows builds need vsearch's own build" >&2; exit 1 ;;
  *)
    os_sources="os/posix/dynlib_loader.cc os/posix/system.cc os/linux/system_memory.cc os/linux/system_cpu.cc os/linux/cgroup.cc"
    [ "$os" = Linux ] && libs="$libs -ldl" ;;
esac

# what autotools would write: config.h for cityhash, and the package
# name and version vsearch prints
mkdir -p "$out/obj" "$out/inc"
rm -f "$out"/obj/*.o
: > "$out/inc/config.h"
printf '#define PACKAGE "vsearch"\n#define PACKAGE_VERSION "dev"\n' > "$out/inc/package.h"

# zlib and bzip2 headers, as configure finds them: vsearch then opens libz
# and libbz2 at run time (dlopen) to read compressed input
flags="-std=c++20 -O2 -DNDEBUG $arch_flags -include $out/inc/package.h"
for h in zlib bzlib; do
  if echo "#include <$h.h>" | $CXX -E -x c++ - > /dev/null 2>&1; then
    flags="$flags -DHAVE_$(echo $h | tr a-z A-Z)_H=1"
  else
    echo "build-vsearch.sh: no $h.h, so no compressed input" >&2
  fi
done

cd "$source_dir"
# Left out: the throwing fatal() and the library API (the CLI uses
# fatal_exit), every platform backend (the chosen ones are added below),
# and cityhash, built with its own flags.
find . \( -name '*.cc' -o -name '*.cpp' \) \
  | sed 's|^\./||' \
  | grep -v -e fatal_throw -e vsearch_api -e '^arch/' -e '^os/[^/]*/' \
      -e '^vendored/city\.cc$' > "$out/sources"
for s in $arch_sources $os_sources; do echo "$s"; done >> "$out/sources"

if [ "${DRY_RUN:-}" = 1 ]; then
  echo "arch=$arch os=$os"
  echo "flags: $flags"
  echo "libs: $libs"
  grep -e '^arch/' -e '^os/' "$out/sources"
  for s in $simd_sources; do echo "${s%%:*} (${s#*:})"; done
  exit 0
fi

obj() { echo "$out/obj/$(echo "$1" | tr / _).o"; }
export CXX flags out
# shellcheck disable=SC2016
xargs -P "$jobs" -I{} sh -c '$CXX $flags -fno-exceptions -I. -Ivendored -c "$1" -o "$2/obj/$(echo "$1" | tr / _).o"' _ {} "$out" < "$out/sources"
for s in $simd_sources; do
  # shellcheck disable=SC2086
  $CXX $flags ${s#*:} -fno-exceptions -I. -c "${s%%:*}" -o "$(obj "${s%%:*}")"
done
# shellcheck disable=SC2086
$CXX $flags -I"$out/inc" -I. -Ivendored -c vendored/city.cc -o "$out/obj/city.o"
$CC -O2 -c vendored/md5.c -o "$out/obj/md5.o"
$CC -O2 -c vendored/sha1.c -o "$out/obj/sha1.o"

# shellcheck disable=SC2086
$CXX "$out"/obj/*.o $libs -o "$out/vsearch"
"$out/vsearch" --version 2>&1 | head -1

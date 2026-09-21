# XDR serialization

This directory contains the Sun RPC / Oracle XDR implementation used by ALPS's Osiris checkpoint serialization. The imported files carry Oracle America copyright notices from 2010 and 2012 and the accompanying BSD-style [license](LICENSE). The original import revision was not recorded.

The four C source files implement generic XDR operations, arrays, floating-point values, and standard-I/O streams. ALPS builds them into `libalps` on every platform to use one wire format and avoid differences between platform RPC types and ABIs. The ALPS C++ dump classes remain in `src/alps/osiris/`.

This is a locally adapted copy. Retained changes include portable fixed-width RPC types, byte-order conversion, Windows support, and ALPS symbol-export declarations. Source includes now refer to the relocated headers. Keep the original notices when modifying these files.

The headers remain available as `<alps/osiris/xdrcore.h>` and `<alps/osiris/rpc_types.h>` because the public ALPS dump classes expose their types. `ALPS::headers` supplies their build include directory; SDK installation places them under `include/alps/osiris/`. The SDK installs the license under `share/alps/licenses/xdr/`, and Python wheels copy it from the SDK into `pyalps/licenses/xdr/`.

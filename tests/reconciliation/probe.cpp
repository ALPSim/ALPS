// SPDX-License-Identifier: MIT
// The same characterization probe is compiled against one provider at a time.
#ifdef PROBE_ALPSCORE
#include <alps/params.hpp>
#else
#include <alps/ngs/params.hpp>
#endif
#include <alps/hdf5/archive.hpp>
#include <alps/hdf5/complex.hpp>
#include <alps/hdf5/map.hpp>
#include <alps/hdf5/vector.hpp>

#include <complex>
#include <climits>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

// Core's supported native long is wide only on LP64 platforms.
constexpr long wide_value = sizeof(long) >= 8 ? static_cast<long>(1099511627776LL) : 42L;

template<class F> void report(std::string const& key, F action) {
    // Evaluate before writing so an exception cannot leave half a record.
    try {
        auto result = action();
        std::cout << key << '\t' << result << '\n';
    } catch (std::exception const& error) {
        std::cout << key << "\tthrows\n";
        std::cerr << key << ": " << error.what() << '\n';
    }
}

template<class T, class P> T get(P const& p, std::string const& key) {
#ifdef PROBE_ALPSCORE
    return p[key].template as<T>();
#else
    return p[key].template cast<T>();
#endif
}

void semantics() {
    report("native_long_bits", [] { return sizeof(long) * CHAR_BIT; });
    report("missing_lookup_size", [] {
        alps::params p;
        (void)p["missing"];
        return p.size();
    });
    report("iteration", [] {
        alps::params p;
        p["z"] = 1; p["a"] = 2;
        std::string keys;
        for (auto it = p.begin(); it != p.end(); ++it)
            keys += (keys.empty() ? "" : ",") + it->first;
        return keys;
    });
    report("erase_missing", [] { alps::params p; p.erase("absent"); return "ok"; });
    report("string_to_int", [] { alps::params p; p["x"] = "4.5"; return get<int>(p, "x"); });
    report("double_to_int", [] { alps::params p; p["x"] = 4.5; return get<int>(p, "x"); });
    report("int_to_double", [] { alps::params p; p["x"] = 4; return get<double>(p, "x"); });
    report("scalar_to_vector", [] {
        alps::params p; p["x"] = 3;
        return get<std::vector<double>>(p, "x") == std::vector<double>{3} ? "ok" : "mismatch";
    });
    report("wide_integer", [] {
        alps::params p; p["x"] = wide_value;
        return get<long>(p, "x");
    });
}

// Each case verifies a value, not merely whether the read completed.
struct Payload {
    alps::hdf5::archive& ar;
    bool writing;
    template<class T> void operator()(std::string const& path, T const& expected) {
        report(path, [&] {
            if (writing) ar[path] << expected;
            else {
                T actual{};
                ar[path] >> actual;
                if (actual != expected) return "mismatch";
            }
            return "ok";
        });
    }
};

void archive_payload(std::string const& filename, bool writing) {
    alps::hdf5::archive ar(filename, writing ? "w" : "r");
    Payload check{ar, writing};
    check("/int", 42);
    check("/double", 1.25);
    check("/bool", true);
    check("/int8", static_cast<signed char>(-7));
    check("/string", std::string("HDF5 \xCE\xBB"));
    check("/empty_string", std::string());
    check("/vector", std::vector<double>{1.25, -2.5, 0});
    check("/empty_vector", std::vector<double>{});
    check("/strings", std::vector<std::string>{"first", "", "last"});
    check("/bools", std::vector<bool>{true, false, true});
    check("/int8s", std::vector<signed char>{-7, 0, 42});
    check("/complex", std::complex<double>{1.25, -2.5});
    check("/complex_vector", std::vector<std::complex<double>>{{1, 2}, {-3, 4}});
    check("/ragged", std::vector<std::vector<double>>{{1, 2}, {3}});
    check("/map", std::map<std::string, double>{{"a/b", 1}, {"c&d", 2}});
    check("/int/@number", 3);
    check("/int/@string", std::string("attribute"));
    check("/int/@strings", std::vector<std::string>{"a", "", "b"});
    check("/int/@vector", std::vector<double>{1, 2});
    check("/int/@bool", true);
    if (!writing) {
        report("bool_type_marker", [&] { return ar.is_attribute("/bool/@__alps_type__"); });
        report("int8_type_marker", [&] { return ar.is_attribute("/int8/@__alps_type__"); });
    }
}

void archive_semantics(std::string const& filename) {
    alps::hdf5::archive ar(filename, "w");
    report("initial_context", [&] { return ar.get_context() == "/" ? "root" : ar.get_context().empty() ? "empty" : "other"; });
    report("empty_map_exists", [&] {
        ar["/empty_map"] << std::map<std::string, double>{};
        return ar.is_group("/empty_map");
    });
    report("overflow_string_to_int", [&] {
        ar["/overflow"] << std::string("99999999999999999999999999999");
        int value = 0;
        ar["/overflow"] >> value;
        return value;
    });
    report("context_complex", [&] {
        ar.set_context("/context_complex");
        ar[""] << std::complex<double>{1, 2};
        std::complex<double> actual;
        ar[""] >> actual;
        return actual == std::complex<double>{1, 2} ? "ok" : "mismatch";
    });
}

void write_params(std::string const& filename, bool extended) {
    alps::params p;
    p["integer"] = 42;
    p["real"] = 1.25;
    p["text"] = "example";
    p["boolean"] = true;
    p["vector"] = std::vector<double>{1, 2};
    if (extended) {
        p["wide"] = wide_value;
        p["unsigned"] = 42U;
        p["float"] = 1.25F;
    }
    alps::hdf5::archive ar(filename, "w");
    report("save", [&] { ar["/parameters"] << p; return "ok"; });
}

template<class P> void read_parameter_values(P& p, alps::hdf5::archive& ar, bool extended) {
    report("load", [&] { ar["/parameters"] >> p; return "ok"; });
    report("integer", [&] { return get<int>(p, "integer"); });
    report("real", [&] { return get<double>(p, "real"); });
    report("text", [&] { return get<std::string>(p, "text"); });
    report("boolean", [&] { return get<bool>(p, "boolean"); });
    report("vector", [&] { return get<std::vector<double>>(p, "vector") == std::vector<double>{1, 2} ? "ok" : "mismatch"; });
    if (extended) {
        report("wide", [&] { return get<long>(p, "wide"); });
        report("unsigned", [&] { return get<unsigned>(p, "unsigned"); });
        report("float", [&] { return get<float>(p, "float"); });
    }
}

int main(int argc, char** argv) {
    try {
        if (argc == 2 && std::string(argv[1]) == "semantics") semantics();
        else if (argc == 3) {
            std::string command = argv[1], file = argv[2];
            if (command == "write-archive") archive_payload(file, true);
            else if (command == "read-archive") archive_payload(file, false);
            else if (command == "archive-semantics") archive_semantics(file);
            else if (command == "write-params" || command == "write-extended-params")
                write_params(file, command == "write-extended-params");
            else if (command == "read-params" || command == "read-extended-params") {
                alps::hdf5::archive ar(file, "r");
                alps::params p;
                read_parameter_values(p, ar, command == "read-extended-params");
            }
#ifdef PROBE_ALPSCORE
            else if (command == "read-dictionary") {
                alps::hdf5::archive ar(file, "r");
                alps::params_ns::dictionary p;
                read_parameter_values(p, ar, false);
            }
#endif
            else throw std::invalid_argument("unknown command");
        } else throw std::invalid_argument("expected semantics, or <command> <file>");
    } catch (std::exception const& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}

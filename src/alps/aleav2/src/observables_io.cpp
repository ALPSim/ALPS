/*
 * Copyright (C) 2026 ALPS Collaboration. See COPYRIGHT.TXT
 * All rights reserved. Use is subject to license terms. See LICENSE.TXT
 */
#include <alps/alea/observables_io.hpp>

#include "observables_impl.hpp"

#include <alps/alea/hdf5.hpp>
#include <alps/alea/plugin/stream_serializer.hpp>
#include <alps/hdf5/archive.hpp>
#include <alps/hdf5/vector.hpp>
#include <alps/utilities/encode.hpp>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <set>

namespace alps { namespace alea {

namespace {

typedef observables_access access;
typedef std::shared_ptr<observable_result::impl> impl_pointer;

const uint64_t RESULT_SCHEMA = 1;
const uint64_t CHECKPOINT_VERSION = 1;
const double nan = std::numeric_limits<double>::quiet_NaN();

std::string signed_name(const std::string &sign, const std::string &name) { return sign + " * " + name; }

// RESULTS, WRITING

void write_values(hdf5::archive &ar, const std::string &path, const std::vector<double> &values, bool scalar)
{
    if (scalar)
        ar[path] << values.at(0);
    else
        ar[path] << values;
}

/** Attributes can only be written to a group that exists */
void need_group(hdf5::archive &ar, const std::string &group)
{
    if (!ar.is_group(group))
        ar.create_group(group);
}

void write_header(hdf5::archive &ar, const std::string &group, const char *kind)
{
    need_group(ar, group);
    ar[group + "/@alea_schema"] << RESULT_SCHEMA;
    ar[group + "/@kind"] << std::string(kind);
    ar[group + "/@value_type"] << std::string("real");
}

void write_statistics(hdf5::archive &ar, const std::string &group, const observable_result &result)
{
    const bool scalar = result.is_scalar();
    write_values(ar, group + "/mean/value", result.mean_vector(), scalar);
    write_values(ar, group + "/mean/error", result.error_vector(), scalar);
    const std::vector<error_convergence> converged = result.converged_vector();
    const std::vector<int> flags(converged.begin(), converged.end());
    if (scalar)
        ar[group + "/mean/error_convergence"] << flags[0];
    else
        ar[group + "/mean/error_convergence"] << flags;
}

/** An observable that is not signed, or the samples of a signed one */
void write_plain(hdf5::archive &ar, const std::string &group, const observable_result &result)
{
    const observable_result::impl &self = access::get(result);
    if (access::is_histogram(self.declaration)) {
        write_header(ar, group, "histogram");
        ar[group + "/count"] << self.count;
        ar[group + "/histogram"] << self.counts;
        // As in ALPS 2, @max is one past the last value counted.
        ar[group + "/@min"] << access::min(self.declaration);
        ar[group + "/@max"] << static_cast<std::int64_t>(access::min(self.declaration) + self.counts.size());
        ar[group + "/@stepsize"] << std::int64_t(1);
        return;
    }

    const bool scalar = result.is_scalar();
    write_header(ar, group, scalar ? "scalar" : "vector");
    ar[group + "/count"] << self.count;
    if (self.count == 0)
        return;
    write_statistics(ar, group, result);
    write_values(ar, group + "/variance/value", result.variance_vector(), scalar);
    write_values(ar, group + "/tau/value", result.tau_vector(), scalar);
    if (!result.labels().empty())
        ar[group + "/labels"] << result.labels();

    if (self.batch) {
        const std::vector<std::vector<double> > bins = result.bin_means();
        if (scalar) {
            std::vector<double> flat;
            for (const std::vector<double> &bin : bins)
                flat.push_back(bin[0]);
            ar[group + "/timeseries/data"] << flat;
        } else {
            ar[group + "/timeseries/data"] << bins;
        }
        ar[group + "/timeseries/data/@binningtype"] << std::string("linear");
        ar[group + "/timeseries/data/@binsize"] << result.bin_size();
        ar[group + "/timeseries/data/@minbinsize"] << uint64_t(1);
        ar[group + "/timeseries/data/@maxbinnum"] << static_cast<uint64_t>(self.batch->store().num_batches());
    }

    // The binning analysis: one row per bin length, one column per component.
    std::vector<std::vector<double> > table;
    for (size_t i = 0; i != result.size(); ++i) {
        const std::vector<double> errors = result.binning_errors(i);
        if (table.size() < errors.size())
            table.resize(errors.size(), std::vector<double>(result.size()));
        for (size_t level = 0; level != errors.size(); ++level)
            table[level][i] = errors[level];
    }
    if (!table.empty()) {
        std::vector<uint64_t> binsize(table.size());
        for (size_t level = 0; level != table.size(); ++level)
            binsize[level] = uint64_t(1) << level;
        ar[group + "/binning/binsize"] << binsize;
        if (scalar) {
            std::vector<double> flat;
            for (const std::vector<double> &row : table)
                flat.push_back(row[0]);
            ar[group + "/binning/error"] << flat;
        } else {
            ar[group + "/binning/error"] << table;
        }
    }
}

// RESULTS, READING

std::vector<double> read_values(hdf5::archive &ar, const std::string &path, size_t size)
{
    if (!ar.is_data(path))
        return std::vector<double>(size, nan);
    std::vector<double> values;
    if (ar.is_scalar(path)) {
        values.resize(1);
        ar[path] >> values[0];
    } else {
        ar[path] >> values;
    }
    return values;
}

impl_pointer read_plain(hdf5::archive &ar, const std::string &group, const std::string &name)
{
    observable declaration(name);
    if (ar.is_attribute(group + "/@sign"))
        ar[group + "/@sign"] >> access::sign(declaration);
    if (ar.is_data(group + "/labels"))
        ar[group + "/labels"] >> access::labels(declaration);
    access::with_bins(declaration) = false;

    if (ar.is_attribute(group + "/@min") && ar.is_attribute(group + "/@max") && ar.is_attribute(group + "/@stepsize")) {
        if (!ar.is_datatype<int>(group + "/@min") && !ar.is_datatype<long>(group + "/@min")
            && !ar.is_datatype<long long>(group + "/@min"))
            throw std::runtime_error("The histogram '" + name + "' is not one of integers; such are not supported");
        access::is_histogram(declaration) = true;
        ar[group + "/@min"] >> access::min(declaration);
        impl_pointer out(new observable_result::impl(declaration));
        ar[group + "/histogram"] >> out->counts;
        access::max(out->declaration) = access::min(out->declaration) + std::int64_t(out->counts.size()) - 1;
        out->size = out->counts.size();
        out->count = std::accumulate(out->counts.begin(), out->counts.end(), uint64_t(0));
        return out;
    }

    impl_pointer out(new observable_result::impl(declaration));
    if (!ar.is_data(group + "/count"))
        throw std::runtime_error("'" + group + "' does not hold a stored observable: it has no count");
    ar[group + "/count"] >> out->count;
    if (out->count == 0 || !ar.is_data(group + "/mean/value")) {
        out->count = 0;
        return out;
    }

    std::shared_ptr<observable_result::impl::stored_values> stored(new observable_result::impl::stored_values);
    stored->mean = read_values(ar, group + "/mean/value", 0);
    out->size = stored->mean.size();
    // A vector with one component is stored as an array, a number is not.
    const bool scalar = ar.is_scalar(group + "/mean/value");
    access::is_vector(out->declaration) = !scalar;
    stored->error = read_values(ar, group + "/mean/error", out->size);
    stored->variance = read_values(ar, group + "/variance/value", out->size);
    stored->tau = read_values(ar, group + "/tau/value", out->size);
    stored->convergence.assign(out->size, MAYBE_CONVERGED);
    if (ar.is_data(group + "/mean/error_convergence")) {
        std::vector<int> flags(1);
        if (ar.is_scalar(group + "/mean/error_convergence"))
            ar[group + "/mean/error_convergence"] >> flags[0];
        else
            ar[group + "/mean/error_convergence"] >> flags;
        for (size_t i = 0; i != out->size && i != flags.size(); ++i)
            stored->convergence[i] = static_cast<error_convergence>(flags[i]);
    }
    if (ar.is_data(group + "/binning/error")) {
        if (scalar) {
            std::vector<double> flat;
            ar[group + "/binning/error"] >> flat;
            for (double error : flat)
                stored->binning_errors.push_back(std::vector<double>(1, error));
        } else {
            ar[group + "/binning/error"] >> stored->binning_errors;
        }
    }
    out->stored = stored;

    // Bins, for functions of results. "timeseries/data" holds bin means in
    // files of this schema, in the evaluated results of ALPS 2 and 3.0 and in
    // those of ALPSCore; the last two mark the group with "@cannotrebin" (or
    // do not record the bin length), as ALEA's mcdata reads them. The
    // observables a 2.x or 3.0 simulation wrote while running hold bin sums.
    const std::string data = group + "/timeseries/data";
    if (ar.is_data(data)) {
        const bool has_binsize = ar.is_attribute(data + "/@binsize");
        const bool sums = has_binsize && !ar.is_attribute(group + "/@alea_schema")
                          && !ar.is_attribute(group + "/@cannotrebin");
        uint64_t binsize = 0;
        if (has_binsize)
            ar[data + "/@binsize"] >> binsize;
        std::vector<std::vector<double> > bins;
        if (scalar) {
            std::vector<double> flat;
            ar[data] >> flat;
            for (double mean : flat)
                bins.push_back(std::vector<double>(1, mean));
        } else {
            ar[data] >> bins;
        }
        if (!has_binsize && !bins.empty())
            binsize = out->count / bins.size();
        if (!bins.empty() && binsize != 0) {
            batch_data<double> store(out->size, bins.size());
            store.set_batch_size(binsize);
            for (size_t i = 0; i != bins.size(); ++i) {
                if (bins[i].size() != out->size)
                    throw std::runtime_error("The bins of '" + name + "' do not have the size of its mean");
                for (size_t k = 0; k != out->size; ++k)
                    store.batch()(k, i) = sums ? bins[i][k] : bins[i][k] * double(binsize);
                store.count()(i) = binsize;
            }
            out->batch.reset(new batch_result<double>(store));
            out->parts.push_back(out->batch);
            access::with_bins(out->declaration) = true;
        }
    }
    return out;
}

// TRANSPORT

/** Bytes in, bytes out, for the stream serializer and for the rest */
class byte_stream
{
public:
    byte_stream() : position_(0) { }
    explicit byte_stream(const std::string &bytes) : bytes_(bytes), position_(0) { }

    const std::string &bytes() const { return bytes_; }

    template <typename T>
    byte_stream &operator<<(const T &value)
    {
        bytes_.append(reinterpret_cast<const char *>(&value), sizeof(T));
        return *this;
    }

    template <typename T>
    byte_stream &operator>>(T &value)
    {
        if (position_ + sizeof(T) > bytes_.size())
            throw std::runtime_error("Packed results end early");
        std::memcpy(&value, bytes_.data() + position_, sizeof(T));
        position_ += sizeof(T);
        return *this;
    }

    byte_stream &operator<<(const std::string &value)
    {
        *this << static_cast<uint64_t>(value.size());
        bytes_.append(value);
        return *this;
    }

    byte_stream &operator>>(std::string &value)
    {
        uint64_t size;
        *this >> size;
        if (position_ + size > bytes_.size())
            throw std::runtime_error("Packed results end early");
        value.assign(bytes_, position_, size);
        position_ += size;
        return *this;
    }

private:
    std::string bytes_;
    size_t position_;
};

const uint64_t PACK_MAGIC = 0x414c454132763031ull;   // "ALEA2v01"

void pack_declaration(byte_stream &out, const observable &declaration_)
{
    observable declaration(declaration_);
    out << declaration.name() << access::sign(declaration) << static_cast<uint64_t>(access::labels(declaration).size());
    for (const std::string &label : access::labels(declaration))
        out << label;
    out << static_cast<uint64_t>(access::num_bins(declaration)) << static_cast<uint64_t>(access::components(declaration))
        << uint64_t(access::with_bins(declaration)) << uint64_t(access::is_vector(declaration))
        << uint64_t(access::is_histogram(declaration)) << access::min(declaration) << access::max(declaration);
}

observable unpack_declaration(byte_stream &in)
{
    std::string name;
    in >> name;
    observable declaration(name);
    uint64_t labels, num_bins, components, with_bins, is_vector, is_histogram;
    in >> access::sign(declaration) >> labels;
    access::labels(declaration).resize(labels);
    for (std::string &label : access::labels(declaration))
        in >> label;
    in >> num_bins >> components >> with_bins >> is_vector >> is_histogram >> access::min(declaration) >> access::max(declaration);
    access::num_bins(declaration) = num_bins;
    access::components(declaration) = components;
    access::with_bins(declaration) = with_bins != 0;
    access::is_vector(declaration) = is_vector != 0;
    access::is_histogram(declaration) = is_histogram != 0;
    return declaration;
}

}

void save_results(hdf5::archive &ar, const std::string &path, const result_set &results)
{
    for (const auto &named : access::results(results)) {
        const observable_result::impl &self = access::get(named.second);
        const std::string group = path + "/" + hdf5_name_encode(named.first);
        if (access::sign(self.declaration).empty()) {
            write_plain(ar, group, named.second);
            continue;
        }

        // A signed observable: the ratio under its name, and the samples as
        // measured under "<sign> * <name>", as ALPS 2 does.
        impl_pointer raw(new observable_result::impl(self.raw ? *self.raw : self));
        if (!self.raw) {
            access::sign(raw->declaration).clear();
            raw->parts.clear();
            if (raw->batch)
                raw->parts.push_back(raw->batch);
            raw->function = nullptr;
            raw->convergence = CONVERGED;
            raw->evaluated.reset();
        }
        write_plain(ar, path + "/" + hdf5_name_encode(signed_name(access::sign(self.declaration), named.first)),
                    access::make_result(raw));

        write_header(ar, group, named.second.is_scalar() ? "scalar" : "vector");
        ar[group + "/@sign"] << access::sign(self.declaration);
        ar[group + "/count"] << self.count;
        if (self.count != 0)
            write_statistics(ar, group, named.second);
        if (!named.second.labels().empty())
            ar[group + "/labels"] << named.second.labels();
    }
}

void save_result(hdf5::archive &ar, const std::string &group, const observable_result &result)
{
    const observable_result::impl &self = access::get(result);
    if (access::sign(self.declaration).empty()) {
        write_plain(ar, group, result);
        return;
    }
    write_header(ar, group, result.is_scalar() ? "scalar" : "vector");
    ar[group + "/count"] << self.count;
    if (self.count != 0)
        write_statistics(ar, group, result);
    if (!result.labels().empty())
        ar[group + "/labels"] << result.labels();
}

observable_result stored_result(const std::string &name, const stored_numbers &numbers)
{
    observable declaration(name);
    access::sign(declaration) = numbers.sign;
    access::labels(declaration) = numbers.labels;
    access::with_bins(declaration) = false;
    impl_pointer out(new observable_result::impl(declaration));
    out->count = numbers.count;
    const size_t size = numbers.mean.size();
    if (numbers.count == 0 || size == 0) {
        out->count = 0;
        return access::make_result(out);
    }
    auto check = [&](size_t other, bool optional, const char *what) {
        if (other != size && !(optional && other == 0)) {
            throw std::invalid_argument("The stored " + std::string(what) + " of '" + name
                                        + "' does not have the size of its mean");
        }
    };
    check(numbers.error.size(), false, "error");
    check(numbers.variance.size(), true, "variance");
    check(numbers.tau.size(), true, "autocorrelation time");
    check(numbers.converged.size(), true, "convergence");
    check(numbers.binning_errors.size(), true, "binning analysis");
    if (numbers.scalar && size != 1)
        throw std::invalid_argument("The stored result '" + name + "' is a number with several components");

    out->size = size;
    access::is_vector(out->declaration) = !numbers.scalar;
    std::shared_ptr<observable_result::impl::stored_values> stored(new observable_result::impl::stored_values);
    stored->mean = numbers.mean;
    stored->error = numbers.error;
    stored->variance = numbers.variance.empty() ? std::vector<double>(size, nan) : numbers.variance;
    stored->tau = numbers.tau.empty() ? std::vector<double>(size, nan) : numbers.tau;
    stored->convergence = numbers.converged.empty()
        ? std::vector<error_convergence>(size, MAYBE_CONVERGED) : numbers.converged;
    // The levels that every component has, as error per level
    if (!numbers.binning_errors.empty()) {
        size_t levels = numbers.binning_errors[0].size();
        for (const std::vector<double> &component : numbers.binning_errors)
            levels = std::min(levels, component.size());
        for (size_t level = 0; level != levels; ++level) {
            std::vector<double> errors(size);
            for (size_t i = 0; i != size; ++i)
                errors[i] = numbers.binning_errors[i][level];
            stored->binning_errors.push_back(errors);
        }
    }
    out->stored = stored;

    return access::make_result(out);
}

observable_result stored_histogram(const std::string &name, std::int64_t min, const std::vector<std::uint64_t> &counts)
{
    if (counts.empty())
        throw std::invalid_argument("The stored histogram '" + name + "' has no entries");
    observable declaration(name);
    access::with_bins(declaration) = false;
    access::is_histogram(declaration) = true;
    access::min(declaration) = min;
    access::max(declaration) = min + std::int64_t(counts.size()) - 1;
    impl_pointer out(new observable_result::impl(declaration));
    out->counts = counts;
    out->size = counts.size();
    out->count = std::accumulate(counts.begin(), counts.end(), uint64_t(0));
    return access::make_result(out);
}

observable_result load_result(hdf5::archive &ar, const std::string &group)
{
    const std::string name = hdf5_name_decode(group.substr(group.find_last_of('/') + 1));
    return access::make_result(read_plain(ar, group, name));
}

result_set load_results(hdf5::archive &ar, const std::string &path)
{
    std::map<std::string, impl_pointer> loaded;
    for (const std::string &child : ar.list_children(path)) {
        if (ar.is_group(path + "/" + child))
            loaded[hdf5_name_decode(child)] = read_plain(ar, path + "/" + child, hdf5_name_decode(child));
    }

    // The samples of signed observables are not results of their own.
    result_set out;
    std::set<std::string> samples;
    for (auto &named : loaded) {
        observable_result::impl &self = *named.second;
        const std::string &sign = access::sign(self.declaration);
        if (sign.empty())
            continue;
        const std::string raw_name = signed_name(sign, named.first);
        samples.insert(raw_name);
        if (!loaded.count(raw_name))
            continue;
        self.raw = loaded[raw_name];
        if (self.raw->batch && loaded.count(sign) && loaded[sign]->batch && loaded[sign]->size == 1) {
            self.parts.assign(1, self.raw->batch);
            self.parts.push_back(loaded[sign]->batch);
            const size_t size = self.raw->size;
            self.function = [size](const std::vector<double> &means) {
                std::vector<double> ratio(means.begin(), means.begin() + size);
                for (double &value : ratio)
                    value /= means[size];
                return ratio;
            };
        }
    }
    for (const auto &named : loaded) {
        if (!samples.count(named.first))
            access::results(out)[named.first] = access::make_result(named.second);
    }
    return out;
}

void save_checkpoint(hdf5::archive &ar, const std::string &path, const observable_set &set)
{
    need_group(ar, path);
    ar[path + "/@alea_state"] << CHECKPOINT_VERSION;
    for (const auto &named : access::observables(set)) {
        const accumulator::impl &self = access::get(named.second);
        observable declaration(self.declaration);
        const std::string group = path + "/" + hdf5_name_encode(named.first);
        need_group(ar, group);
        ar[group + "/@kind"] << std::string(access::is_histogram(self.declaration) ? "histogram"
                                             : access::with_bins(self.declaration) ? "binned" : "simple");
        ar[group + "/@sign"] << access::sign(declaration);
        ar[group + "/@bins"] << static_cast<uint64_t>(access::num_bins(declaration));
        ar[group + "/@components"] << static_cast<uint64_t>(access::components(declaration));
        ar[group + "/@vector"] << uint64_t(access::is_vector(declaration));
        ar[group + "/@size"] << static_cast<uint64_t>(self.size);
        ar[group + "/@min"] << access::min(declaration);
        ar[group + "/@max"] << access::max(declaration);
        if (!access::labels(declaration).empty())
            ar[group + "/labels"] << access::labels(declaration);
        if (access::is_histogram(self.declaration))
            ar[group + "/counts"] << self.counts;

        hdf5_serializer state(ar, group);
        if (self.simple)
            self.simple->save_state(state, "simple");
        if (self.autocorr) {
            self.autocorr->save_state(state, "autocorr");
            self.batch->save_state(state, "batch");
        }
    }
}

void load_checkpoint(hdf5::archive &ar, const std::string &path, observable_set &set)
{
    uint64_t version = 0;
    if (ar.is_attribute(path + "/@alea_state"))
        ar[path + "/@alea_state"] >> version;
    if (version != CHECKPOINT_VERSION) {
        throw std::runtime_error(
            "'" + path + "' does not hold measurements that this version of ALPS can continue"
            + (version == 0 ? std::string(": they were written by an earlier version, or are results.")
                            : ": they have state version " + std::to_string(version) + ".")
            + " Results of a finished run can still be read; to add to them, start a new run.");
    }
    for (const std::string &child : ar.list_children(path)) {
        const std::string group = path + "/" + child;
        if (!ar.is_group(group))
            continue;
        observable declaration(hdf5_name_decode(child));
        std::string kind;
        uint64_t num_bins, size, components;
        ar[group + "/@kind"] >> kind;
        ar[group + "/@sign"] >> access::sign(declaration);
        ar[group + "/@bins"] >> num_bins;
        ar[group + "/@components"] >> components;
        access::components(declaration) = components;
        if (ar.is_attribute(group + "/@vector")) {
            uint64_t is_vector;
            ar[group + "/@vector"] >> is_vector;
            access::is_vector(declaration) = is_vector != 0;
        }
        ar[group + "/@size"] >> size;
        ar[group + "/@min"] >> access::min(declaration);
        ar[group + "/@max"] >> access::max(declaration);
        if (ar.is_data(group + "/labels"))
            ar[group + "/labels"] >> access::labels(declaration);
        access::num_bins(declaration) = num_bins;
        access::with_bins(declaration) = kind == "binned";
        access::is_histogram(declaration) = kind == "histogram";

        accumulator loaded = access::make_accumulator(declaration);
        accumulator::impl &self = access::get(loaded);
        if (kind == "histogram")
            ar[group + "/counts"] >> self.counts;
        if (size != 0) {
            self.set_size(size);
            hdf5_serializer state(ar, group);
            if (self.simple)
                self.simple->load_state(state, "simple");
            if (self.autocorr) {
                self.autocorr->load_state(state, "autocorr");
                self.batch->load_state(state, "batch");
            }
        }
        access::observables(set).erase(declaration.name());
        access::observables(set).insert(std::make_pair(declaration.name(), loaded));
    }
}

std::string pack(const result_set &results)
{
    byte_stream out;
    stream_serializer<byte_stream> serializer(out);
    out << PACK_MAGIC << static_cast<uint64_t>(access::results(results).size());
    for (const auto &named : access::results(results)) {
        const observable_result::impl &self = access::get(named.second);
        if (self.stored)
            throw std::invalid_argument("'" + named.first + "' was read from a result file and cannot be packed");
        observable declaration(self.declaration);
        access::name(declaration) = named.first;
        pack_declaration(out, declaration);
        out << static_cast<uint64_t>(self.size) << self.count << static_cast<uint64_t>(self.counts.size());
        for (uint64_t count : self.counts)
            out << count;
        out << uint64_t(bool(self.simple)) << uint64_t(bool(self.autocorr));
        if (self.simple)
            serialize(serializer, "", *self.simple);
        if (self.autocorr) {
            serialize(serializer, "", *self.autocorr);
            serialize(serializer, "", *self.batch);
        }
    }
    return out.bytes();
}

result_set unpack(const std::string &bytes)
{
    byte_stream in(bytes);
    stream_deserializer<byte_stream> deserializer(in);
    uint64_t magic, number;
    in >> magic >> number;
    if (magic != PACK_MAGIC)
        throw std::runtime_error("These bytes are not packed results of this version of ALPS");
    std::map<std::string, impl_pointer> parts;
    for (uint64_t n = 0; n != number; ++n) {
        impl_pointer self(new observable_result::impl(unpack_declaration(in)));
        uint64_t size, counts, has_simple, has_autocorr;
        in >> size >> self->count >> counts;
        self->size = size;
        self->counts.resize(counts);
        for (uint64_t &count : self->counts)
            in >> count;
        in >> has_simple >> has_autocorr;
        if (has_simple) {
            var_result<double> simple;
            deserialize(deserializer, "", simple);
            self->simple.reset(new var_result<double>(simple));
        }
        if (has_autocorr) {
            autocorr_result<double> autocorr;
            batch_result<double> batch;
            deserialize(deserializer, "", autocorr);
            deserialize(deserializer, "", batch);
            self->autocorr.reset(new autocorr_result<double>(autocorr));
            self->batch.reset(new batch_result<double>(batch));
            self->parts.push_back(self->batch);
        }
        parts[self->declaration.name()] = self;
    }
    result_set out;
    access::link(out, parts);
    return out;
}

void observable_set::save(hdf5::archive &ar) const { save_checkpoint(ar, ar.get_context(), *this); }

void observable_set::load(hdf5::archive &ar) { load_checkpoint(ar, ar.get_context(), *this); }

void result_set::save(hdf5::archive &ar) const { save_results(ar, ar.get_context(), *this); }

void result_set::load(hdf5::archive &ar) { *this = load_results(ar, ar.get_context()); }

}} /* namespace alps::alea */

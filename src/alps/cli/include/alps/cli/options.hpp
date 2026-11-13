// SPDX-License-Identifier: MIT
#ifndef ALPS_CLI_OPTIONS_HPP
#define ALPS_CLI_OPTIONS_HPP

#include <alps/cli_export.h>

#include <charconv>
#include <cstddef>
#include <cstdlib>
#include <functional>
#include <iosfwd>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace alps {
namespace cli {

/// A command line that does not fit the declared options.
class error : public std::logic_error {
public:
  using std::logic_error::logic_error;
};

namespace detail {

inline bool convert(const std::string& text, std::string& value) { value = text; return true; }

inline bool convert(const std::string& text, double& value) {
  // strtod would also skip leading blanks and read hexadecimal.
  if (text.empty() || text.find_first_of(" \t\n\v\f\rxX") != std::string::npos) return false;
  char* end = nullptr;
  value = std::strtod(text.c_str(), &end);
  return end == text.c_str() + text.size();
}

/// Whole decimal integers that fit T. A leading '+' is accepted; a negative
/// value for an unsigned T is not.
template<class T>
typename std::enable_if<std::is_integral<T>::value, bool>::type
convert(const std::string& text, T& value) {
  const char* first = text.c_str() + (text.size() > 1 && text[0] == '+' && text[1] != '-' ? 1 : 0);
  const char* last = text.c_str() + text.size();
  auto result = std::from_chars(first, last, value);
  return result.ec == std::errc() && result.ptr == last;
}

} // namespace detail

/// The options of one program, and the values a command line gave them.
///
/// The grammar is the one ALPS programs have always had:
///   --name value, --name=value, -n value, -nvalue, grouped flags (-cv, -cT5);
///   a long name may be shortened to any prefix that names one option;
///   "--" ends the options, and everything after it is positional;
///   an option may be given once, unless it was declared with list().
/// Names are case-sensitive. Any violation throws alps::cli::error.
class ALPS_CLI_DECL options {
public:
  explicit options(const std::string& caption = "Allowed options");

  /// Start a group of options with its own caption in the help text.
  options& section(const std::string& caption);

  /// Declare an option. names is "long" or "long,s", s being the short letter.
  options& flag(const std::string& names, const std::string& help);

  /// An option with one value of type T. Read it back with count() and get<T>().
  template<class T>
  options& value(const std::string& names, const std::string& help) {
    return add(names, help, [](const std::string& text) { T checked; return detail::convert(text, checked); });
  }
  /// As above, and the value is also written to target, which parse() must
  /// find alive.
  template<class T>
  options& value(const std::string& names, T& target, const std::string& help) {
    return add(names, help, [&target](const std::string& text) { return detail::convert(text, target); });
  }
  /// As above; target starts as fallback, which the help text shows.
  template<class T>
  options& value(const std::string& names, T& target, const typename std::common_type<T>::type& fallback,
                 const std::string& help) {
    target = fallback;
    std::ostringstream shown;
    shown << fallback;
    return value(names, target, help).fallback_text(shown.str());
  }
  /// An option that may be repeated; each use adds one string to values().
  options& list(const std::string& names, const std::string& help);

  /// Leave the option declared last out of the help text.
  options& hidden();

  /// Arguments that are not options go, in the order of these calls, to the
  /// named options: the next max_count of them, or all that remain if negative.
  options& positional(const std::string& name, int max_count = 1);

  /// Read argv[1] .. argv[argc - 1].
  void parse(int argc, const char* const* argv);

  /// How often the option was given, by its long name; 0 if it is not declared.
  std::size_t count(const std::string& name) const;
  /// The strings given to the option, in order.
  const std::vector<std::string>& values(const std::string& name) const;
  /// The value of an option that was given; check count() first.
  template<class T>
  T get(const std::string& name) const {
    T result;
    const std::vector<std::string>& given = values(name);
    if (given.empty() || !detail::convert(given.front(), result))
      throw error("option '--" + name + "' has no value");
    return result;
  }

  /// Write the help text, one entry per option.
  void print(std::ostream& os) const;

private:
  struct option {
    std::string name;      // long name, without dashes
    char letter;           // short name, or 0
    std::string help, caption, fallback;
    bool takes_value, repeatable, hidden;
    std::function<bool(const std::string&)> assign;
    std::vector<std::string> given;
    std::size_t seen;
  };
  options& add(const std::string& names, const std::string& help, std::function<bool(const std::string&)> assign);
  options& fallback_text(const std::string& text);
  option& find(const std::string& name);
  const option& find(const std::string& name) const;
  option& match_long(const std::string& name);
  option* match_short(char letter);
  void store(option& opt, const std::string* text);

  std::string caption_;
  std::vector<option> options_;
  std::vector<std::pair<std::string, int> > positional_;
};

inline std::ostream& operator<<(std::ostream& os, const options& opts) {
  opts.print(os);
  return os;
}

} // namespace cli
} // namespace alps

#endif
